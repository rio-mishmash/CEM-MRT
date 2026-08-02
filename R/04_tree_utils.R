# ==============================================================================
# R/04_tree_utils.R: Helper Functions for Tree-Based Stratification and MRT
# ==============================================================================

#' Prepare Tree Metadata from a Ranger Forest
#' @param ranger_model A trained ranger object
#' @param tree_indices Vector of tree indices to process
#' @return A list of tree metadata (parent maps, descendant lists, depths)
prepare_tree_info <- function(ranger_model, tree_indices) {
  lapply(tree_indices, function(k) {
    t_info <- ranger::treeInfo(ranger_model, k)
    max_node <- max(t_info$nodeID, na.rm = TRUE)
    
    pmap <- rep(NA, max_node + 1)
    for (j in 1:nrow(t_info)) {
      if (!is.na(t_info$leftChild[j]))  pmap[t_info$leftChild[j]  + 1] <- t_info$nodeID[j]
      if (!is.na(t_info$rightChild[j])) pmap[t_info$rightChild[j] + 1] <- t_info$nodeID[j]
    }
    
    node_depth <- rep(0, max_node + 1)
    for (n in 0:max_node) {
      curr <- n
      d <- 0
      while (!is.na(curr) && curr >= 0 && curr <= max_node && !is.na(pmap[curr + 1])) {
        curr <- pmap[curr + 1]
        d <- d + 1
        if (d > 100) break 
      }
      node_depth[n + 1] <- d
    }
    
    d_list <- vector("list", max_node + 1)
    sorted_rows <- order(t_info$nodeID, decreasing = TRUE)
    for (j in sorted_rows) {
      n <- t_info$nodeID[j]
      if (t_info$terminal[j]) {
        d_list[[n + 1]] <- numeric(0)
      } else {
        lc <- t_info$leftChild[j]
        rc <- t_info$rightChild[j]
        lc_desc <- if (!is.na(lc) && lc <= max_node) c(lc, d_list[[lc + 1]]) else numeric(0)
        rc_desc <- if (!is.na(rc) && rc <= max_node) c(rc, d_list[[rc + 1]]) else numeric(0)
        d_list[[n + 1]] <- unique(c(lc_desc, rc_desc))
      }
    }
    
    list(parent_map = pmap, descendants_list = d_list, node_depth = node_depth, tree_df = t_info)
  })
}

#' Calculate Variable Importance for Pruned Forest
#' @param ranger_model The original ranger model
#' @param tree_indices Indices of trees being used
#' @param tree_info_list Output from prepare_tree_info
#' @param final_node_matrix Matrix of terminal node IDs per observation
calculate_pruned_importance <- function(ranger_model, tree_indices, tree_info_list, final_node_matrix) {
  all_vars  <- ranger_model$forest$independent.variable.names
  importance_scores <- setNames(numeric(length(all_vars)), all_vars)
  
  K <- length(tree_indices)
  for (i in 1:K) {
    k <- tree_indices[i]
    t_info <- ranger::treeInfo(ranger_model, k)
    
    col_idx <- if (ncol(final_node_matrix) == K) i else k
    active_leaves <- unique(final_node_matrix[, col_idx])
    active_leaves <- active_leaves[!is.na(active_leaves)]
    if (length(active_leaves) == 0) next
    
    max_node_id <- max(t_info$nodeID)
    parent_map <- integer(max_node_id + 1)
    valid_left <- t_info$leftChild[!is.na(t_info$leftChild)] + 1
    valid_right <- t_info$rightChild[!is.na(t_info$rightChild)] + 1
    parent_map[valid_left] <- which(!is.na(t_info$leftChild))
    parent_map[valid_right] <- which(!is.na(t_info$rightChild))
    
    active_split_nodes <- c()
    for (leaf_id in active_leaves) {
      curr_id <- leaf_id
      while (TRUE) {
        parent_row <- parent_map[curr_id + 1]
        if (parent_row == 0) break
        p_id <- t_info$nodeID[parent_row]
        active_split_nodes <- c(active_split_nodes, p_id)
        curr_id <- p_id
      }
    }
    active_split_nodes <- unique(active_split_nodes)
    
    if (length(active_split_nodes) > 0) {
      node_row_map <- integer(max_node_id + 1)
      node_row_map[t_info$nodeID + 1] <- 1:nrow(t_info)
      for (node_id in active_split_nodes) {
        node_row <- node_row_map[node_id + 1]
        v_name <- as.character(t_info$splitvarName[node_row])
        stat <- t_info$splitStat[node_row]
        if (!is.na(stat)) {
          importance_scores[v_name] <- importance_scores[v_name] + stat
        }
      }
    }
  }
  
  total_imp <- sum(importance_scores)
  if (total_imp > 0) importance_scores <- importance_scores / total_imp
  return(importance_scores)
}

#' Prune Forest Structure to Reach Target Sample Size
#' @param initial_nodes Matrix of terminal nodes for each observation
#' @param tree_info_list Output from prepare_tree_info
#' @param df Data frame with covariates and treatment Z
#' @param target_n The desired matched sample size
#' @param tolerance Allowable shortfall from target_n to prevent over-pruning
prune_forest_to_target <- function(initial_nodes, tree_info_list, df, target_n, tolerance = 10) {
  K <- length(tree_info_list)
  current_nodes <- as.matrix(initial_nodes)
  idx_z1 <- which(df$Z == 1)
  idx_z0 <- which(df$Z == 0)
  X_cols <- grep("^X", names(df), value = TRUE)
  
  Sigma <- cov(df[, X_cols, drop = FALSE])
  Sigma_inv <- MASS::ginv(Sigma)
  
  # 1. Pre-calculate branch scores (Mahalanobis Distance) for pruning
  branch_list <- list()
  for (i in 1:K) {
    tree_info <- tree_info_list[[i]]
    parent_nodes <- unique(tree_info$parent_map[!is.na(tree_info$parent_map)])
    
    for (p in parent_nodes) {
      # Find children IDs for node p
      children <- which(tree_info$parent_map == p) - 1
      
      # Ensure there are at least two valid split branches
      if (length(children) >= 2) {
        # Extract left and right child IDs explicitly
        child_1 <- children[1]
        child_2 <- children[2]
        
        # Correctly define two separate branches to compare
        # Pass a single integer into [[ ]] to avoid out-of-bounds hierarchical indexing
        desc_1 <- c(child_1, tree_info$descendants_list[[child_1 + 1]])
        desc_2 <- c(child_2, tree_info$descendants_list[[child_2 + 1]])
        
        idx_1 <- which(initial_nodes[, i] %in% desc_1)
        idx_2 <- which(initial_nodes[, i] %in% desc_2)
        
        if (length(idx_1) > 0 && length(idx_2) > 0) {
          mean_1 <- colMeans(df[idx_1, X_cols, drop = FALSE])
          mean_2 <- colMeans(df[idx_2, X_cols, drop = FALSE])
          diff_vec <- as.numeric(mean_1 - mean_2)
          
          # Calculate Mahalanobis distance between the two split branches
          md_score <- as.numeric(t(diff_vec) %*% Sigma_inv %*% diff_vec)
          
          branch_list[[length(branch_list) + 1]] <- list(
            tree_idx = i, parent = p, children = children, score = md_score
          )
        }
      }
    }
  }
  
  scores <- sapply(branch_list, function(x) x$score)
  branch_ranking <- branch_list[order(scores)]
  
  # 2. Pruning Loop
  strata_chars <- do.call(paste, c(as.data.frame(current_nodes), sep = "_"))
  
  count_actual_n <- function(s_chars) {
    s_int <- as.integer(factor(s_chars))
    u_strata <- unique(s_int)
    c1 <- tabulate(s_int[idx_z1], nbins = length(u_strata))
    c0 <- tabulate(s_int[idx_z0], nbins = length(u_strata))
    return(sum((c1 >= 2 & c0 >= 2)[s_int]))
  }
  
  actual_n <- count_actual_n(strata_chars)
  best_diff <- abs(actual_n - target_n)
  best_snapshot <- list(strata = as.integer(factor(strata_chars)), actual_n = actual_n, final_nodes = current_nodes)
  
  # Early return if actual_n is already within acceptable range (target_n - tolerance)
  if (actual_n >= (target_n - tolerance)) return(best_snapshot)
  
  for (step in 1:length(branch_ranking)) {
    task <- branch_ranking[[step]]
    
    # Extract left and right child IDs explicitly to avoid vector indexing inside [[ ]]
    child_1 <- task$children[1]
    child_2 <- task$children[2]
    
    # Shortcut to the target tree's descendants list
    target_desc_list <- tree_info_list[[task$tree_idx]]$descendants_list
    
    # Combine both child nodes and all their respective descendants safely
    descendants <- c(
      task$children, 
      target_desc_list[[child_1 + 1]], 
      target_desc_list[[child_2 + 1]]
    )
    
    update_idx <- which(current_nodes[, task$tree_idx] %in% descendants)
    if (length(update_idx) > 0) {
      current_nodes[update_idx, task$tree_idx] <- task$parent
      strata_chars[update_idx] <- do.call(paste, c(as.data.frame(current_nodes[update_idx, , drop = FALSE]), sep = "_"))
      actual_n <- count_actual_n(strata_chars)
      current_diff <- abs(actual_n - target_n)
      
      if (current_diff < best_diff) {
        best_diff <- current_diff
        best_snapshot <- list(strata = as.integer(factor(strata_chars)), actual_n = actual_n, final_nodes = current_nodes)
      }
      
      # Early return if the sample size reaches target_n - tolerance
      if (actual_n >= (target_n - tolerance)) return(best_snapshot)
    }
  }
  return(best_snapshot)
}


#' Select the Most Representative Tree (MRT) using WSV Distance
#' @param ranger_model A trained ranger random forest object
#' @param data Data frame used for training
#' @return A list containing the MRT index, distance scores, and structure details
get_representative_tree <- function(ranger_model, data) {
  num_trees <- ranger_model$num.trees
  all_vars  <- ranger_model$forest$independent.variable.names
  num_vars  <- length(all_vars)
  
  # Initialize the U matrix (Trees x Variables) to store variable split weights
  U_matrix <- matrix(0, nrow = num_trees, ncol = num_vars, dimnames = list(NULL, all_vars))
  
  for (i in 1:num_trees) {
    t_info <- ranger::treeInfo(ranger_model, i)
    
    # 1. Identify valid split nodes (where splitvarName is not NA)
    split_nodes <- t_info[!is.na(t_info$splitvarName), ]
    
    if (nrow(split_nodes) > 0) {
      max_node <- max(t_info$nodeID)
      
      # 2. Build parent node tracking map
      # pmap[child_nodeID + 1] stores its parent_nodeID
      pmap <- rep(NA, max_node + 1)
      for (j in 1:nrow(t_info)) {
        if (!is.na(t_info$leftChild[j])) {
          pmap[t_info$leftChild[j] + 1] <- t_info$nodeID[j]
        }
        if (!is.na(t_info$rightChild[j])) {
          pmap[t_info$rightChild[j] + 1] <- t_info$nodeID[j]
        }
      }
      
      # 3. Calculate node depths securely using a complete while-loop traceback
      depth_map <- numeric(max_node + 1)
      for (nid in t_info$nodeID) {
        d <- 1
        curr <- nid
        # Trace up to the root node (root has no parent, so pmap will be NA)
        while (!is.na(pmap[curr + 1])) {
          curr <- pmap[curr + 1]
          d <- d + 1
        }
        depth_map[nid + 1] <- d
      }
      
      # 4. Extract max depth among valid split nodes
      max_depth <- max(depth_map[split_nodes$nodeID + 1])
      
      # 5. Populate weights into the U matrix
      for (k in 1:nrow(split_nodes)) {
        v_name <- as.character(split_nodes$splitvarName[k])
        if (v_name %in% all_vars) {
          nid   <- split_nodes$nodeID[k]
          level <- depth_map[nid + 1]
          U_matrix[i, v_name] <- U_matrix[i, v_name] + (1 / (2^(level - 1)))
        }
      }
      
      # 6. Normalize the row vector by max depth
      if (max_depth > 0) {
        U_matrix[i, ] <- U_matrix[i, ] / max_depth
      }
    }
  }
  
  # Compute pairwise Euclidean distances between tree vectors (WSV space)
  dist_matrix <- as.matrix(dist(U_matrix, method = "euclidean"))
  total_distances <- rowSums(dist_matrix)
  
  # Identify the centroid tree (Representative Tree) minimizing global distance
  rep_tree_idx <- which.min(total_distances)
  
  return(list(
    index = as.integer(rep_tree_idx),
    scores = total_distances,
    U_matrix = U_matrix,
    tree_structure = ranger::treeInfo(ranger_model, tree = as.integer(rep_tree_idx))
  ))
}