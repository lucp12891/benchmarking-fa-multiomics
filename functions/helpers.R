extract_id_from_feature <- function(df, feature_col = "feature") {
  df$ID <- as.numeric(sub(".*omic[12]_feature_(\\d+)", "\\1", df[[feature_col]]))
  df
}

extract_id_from_sample <- function(df, sample_col = "sample") {
  df$ID <- as.numeric(sub(".*sample_(\\d+)", "\\1", df[[sample_col]]))
  df
}

run_with_seed <- function(expr, seed = NULL) {
  old_seed <- .Random.seed
  set.seed(seed)
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  force(expr)
}

align_loading_columns_by_truth <- function(df, true_vector) {
  # Identify loading columns
  loading_cols <- names(df)[grepl("loading", names(df))]
  
  # Extract only loading columns
  loading_matrix <- df[, loading_cols, drop = FALSE]
  
  # Check dimensions
  if (nrow(loading_matrix) != length(true_vector)) {
    stop("Number of rows in loading columns must match length of true_vector.")
  }
  
  # Compute correlation with true vector
  cor_vals <- sapply(loading_matrix, function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order columns by absolute correlation (descending)
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- loading_cols[ordered_idx]
  
  # Extract method name (e.g., "FABIA" from "loading_FABIA1")
  method_name <- sub(".*loading_([A-Za-z]+).*", "\\1", ordered_cols[1])
  
  # Rename columns based on match order
  new_names <- paste0("loading_", method_name, seq_along(ordered_cols))
  
  # Replace names in original df
  df_renamed <- df %>%
    dplyr::select(-all_of(loading_cols)) %>%
    dplyr::bind_cols(setNames(loading_matrix[, ordered_idx], new_names))
  
  return(df_renamed)
}

align_loading_columns_by_truth_overall <- function(weights_df, true_df, method_name) {
  # Merge by feature
  df_merged <- merge(true_df, weights_df, by = "feature")
  
  # Pull true values
  true_vector <- df_merged$true_loading_F1
  
  # Extract only numeric loading columns for this method
  loading_cols <- grep(paste0("^loading_", method_name, "[0-9]+$"), names(df_merged), value = TRUE)
  
  if (length(loading_cols) < 2) {
    stop(paste("Expected at least 2 loading columns for method", method_name))
  }
  
  # Ensure loading columns are numeric
  loading_matrix <- df_merged#[, loading_cols, drop = FALSE]
  for (col in loading_cols) {
    if (!is.numeric(loading_matrix[[col]])) {
      loading_matrix[[col]] <- as.numeric(as.character(loading_matrix[[col]]))
    }
  }
  
  # Compute correlations
  cor_vals <- sapply(loading_matrix[, loading_cols], function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order by absolute correlation
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- loading_cols[ordered_idx]
  
  new_names <- paste0("matched_", method_name, seq_along(ordered_cols))
  old_names <- character(length(ordered_idx))
  for (j in seq_along(ordered_idx)) {
    old_names[j] <- paste0("loading_", method_name, ordered_idx[j])
  }
  
  # Rename the columns in df_merged directly
  names(df_merged)[match(old_names, names(df_merged))] <- new_names
  
  return(df_merged)
}

align_scores_columns_by_truth <- function(df, true_vector) {
  # Identify score columns
  score_cols <- names(df)[grepl("score", names(df))]
  
  # Extract only score columns
  score_matrix <- df[, score_cols, drop = FALSE]
  
  # Check dimensions
  if (nrow(score_matrix) != length(true_vector)) {
    stop("Number of rows in score columns must match length of true_vector.")
  }
  
  # Compute correlation with true vector
  cor_vals <- sapply(score_matrix, function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order columns by absolute correlation (descending)
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- score_cols[ordered_idx]
  
  # Extract method name (e.g., "FABIA" from "loading_FABIA1")
  method_name <- sub(".*score_([A-Za-z]+).*", "\\1", ordered_cols[1])
  
  # Rename columns based on match order
  new_names <- paste0("score_", method_name, seq_along(ordered_cols))
  
  # Replace names in original df
  df_renamed <- df %>%
    dplyr::select(-all_of(score_cols)) %>%
    dplyr::bind_cols(setNames(score_matrix[, ordered_idx], new_names))
  
  return(df_renamed)
}


align_score_columns_by_truth_overall <- function(scores_df, true_df, method_name) {
  # Merge by feature
  df_merged <- merge(true_df, scores_df, by = "sample")
  
  # Pull true values
  true_vector <- df_merged$true_score_F1
  
  # Extract only numeric score columns for this method
  score_cols <- grep(paste0("^score_", method_name, "[0-9]+$"), names(df_merged), value = TRUE)
  
  if (length(score_cols) < 2) {
    stop(paste("Expected at least 2 score columns for method", method_name))
  }
  
  # Ensure score columns are numeric
  score_matrix <- df_merged#[, loading_cols, drop = FALSE]
  for (col in score_cols) {
    if (!is.numeric(score_matrix[[col]])) {
      score_matrix[[col]] <- as.numeric(as.character(score_matrix[[col]]))
    }
  }
  
  # Compute correlations
  cor_vals <- sapply(score_matrix[, score_cols], function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order by absolute correlation
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- score_cols[ordered_idx]
  
  new_names <- paste0("matched_", method_name, seq_along(ordered_cols))
  old_names <- character(length(ordered_idx))
  for (j in seq_along(ordered_idx)) {
    old_names[j] <- paste0("score_", method_name, ordered_idx[j])
  }
  
  # Rename the columns in df_merged directly
  names(df_merged)[match(old_names, names(df_merged))] <- new_names
  
  return(df_merged)
}

# ----------------- Align Loadings -----------------
align_loadings <- function(weights_df, truth, method_name = NULL, truth_col = NULL) {
  # Case 1: truth is a vector
  if (is.vector(truth)) {
    true_vector <- truth
    loading_cols <- grep("loading", names(weights_df), value = TRUE)
    if (length(loading_cols) < 1) stop("No loading columns found.")
    if (nrow(weights_df) != length(true_vector)) {
      stop("Length of truth vector must match number of rows in weights_df.")
    }
    loading_matrix <- weights_df[, loading_cols, drop = FALSE]
    
    # correlations
    cor_vals <- sapply(loading_matrix, function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
    ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
    ordered_cols <- loading_cols[ordered_idx]
    method <- if (is.null(method_name)) sub(".*loading_([A-Za-z]+).*", "\\1", ordered_cols[1]) else method_name
    new_names <- paste0("loading_", method, seq_along(ordered_cols))
    
    df_out <- dplyr::bind_cols(
      weights_df %>% dplyr::select(-all_of(loading_cols)),
      setNames(loading_matrix[, ordered_idx], new_names)
    )
    return(df_out)
  }
  
  # Case 2: truth is a dataframe
  if (is.data.frame(truth)) {
    if (is.null(truth_col)) stop("Please provide truth_col when truth is a dataframe.")
    df_merged <- merge(truth, weights_df, by = "feature")
    true_vector <- df_merged[[truth_col]]
    loading_cols <- grep(paste0("^loading_", method_name, "[0-9]+$"), names(df_merged), value = TRUE)
    if (length(loading_cols) < 1) stop("No loading columns found for method.")
    
    # correlations
    cor_vals <- sapply(df_merged[, loading_cols, drop = FALSE], function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
    ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
    ordered_cols <- loading_cols[ordered_idx]
    
    new_names <- paste0("matched_", method_name, seq_along(ordered_cols))
    names(df_merged)[match(ordered_cols, names(df_merged))] <- new_names
    return(df_merged)
  }
  
  stop("truth must be either a vector or a dataframe.")
}

# ----------------- Align Scores -----------------
align_scores <- function(scores_df, truth, method_name = NULL, truth_col = NULL) {
  # Case 1: truth is a vector
  if (is.vector(truth)) {
    true_vector <- truth
    score_cols <- grep("score", names(scores_df), value = TRUE)
    if (length(score_cols) < 1) stop("No score columns found.")
    if (nrow(scores_df) != length(true_vector)) {
      stop("Length of truth vector must match number of rows in scores_df.")
    }
    score_matrix <- scores_df[, score_cols, drop = FALSE]
    
    # correlations
    cor_vals <- sapply(score_matrix, function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
    ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
    ordered_cols <- score_cols[ordered_idx]
    method <- if (is.null(method_name)) sub(".*score_([A-Za-z]+).*", "\\1", ordered_cols[1]) else method_name
    new_names <- paste0("score_", method, seq_along(ordered_cols))
    
    df_out <- dplyr::bind_cols(
      scores_df %>% dplyr::select(-all_of(score_cols)),
      setNames(score_matrix[, ordered_idx], new_names)
    )
    return(df_out)
  }
  
  # Case 2: truth is a dataframe
  if (is.data.frame(truth)) {
    if (is.null(truth_col)) stop("Please provide truth_col when truth is a dataframe.")
    df_merged <- merge(truth, scores_df, by = "sample")
    true_vector <- df_merged[[truth_col]]
    score_cols <- grep(paste0("^score_", method_name, "[0-9]+$"), names(df_merged), value = TRUE)
    if (length(score_cols) < 1) stop("No score columns found for method.")
    
    # correlations
    cor_vals <- sapply(df_merged[, score_cols, drop = FALSE], function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
    ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
    ordered_cols <- score_cols[ordered_idx]
    
    new_names <- paste0("matched_", method_name, seq_along(ordered_cols))
    names(df_merged)[match(ordered_cols, names(df_merged))] <- new_names
    return(df_merged)
  }
  
  stop("truth must be either a vector or a dataframe.")
}

# ------------------------- Apply Cutoffs -------------------------
apply_cutoffs <- function(aligned_list, varphi_function, type = c("loading","score"), sigma = NULL) {
  type <- match.arg(type)
  vf <- get(varphi_function, mode = "function")
  
  out <- list()
  
  for (method in names(aligned_list)) {
    df <- aligned_list[[method]]
    
    # pick only matched cols
    cols <- grep(paste0("^matched_", method), names(df), value = TRUE)
    
    for (j in seq_along(cols)) {
      colname <- cols[j]
      newname <- paste0(tolower(method), "_index_", letters[j])
      
      if (identical(varphi_function, "varphi.five")) {
        # pass context arguments
        df[[newname]] <- ifelse(
          abs(df[[colname]]) >= vf(
            loading_or_score = df[[colname]],
            factor = j,
            type = type
          ), 1, 0
        )
      } else {
        df[[newname]] <- ifelse(abs(df[[colname]]) >= vf(df[[colname]]), 1, 0)
      }
    }
    
    out[[method]] <- df
  }
  
  # merge all methods back
  Reduce(function(x,y) merge(x,y, by = intersect(names(x),names(y)), all=TRUE), out)
}

# ------------------------- Compare Jaccard -------------------------
CompareJaccard <- function(FeatureData, SampleData) {
  methods <- c("fabia","mofa","mfa","gfa")
  
  ji_features <- lapply(methods, function(m) {
    jaccard.index.sim(FeatureData[[paste0(m,"_index_a")]], FeatureData$signal_index_a)
  })
  names(ji_features) <- paste0(methods, "_features_a")
  
  ji_samples <- lapply(methods, function(m) {
    jaccard.index.sim(SampleData[[paste0(m,"_index_a")]], SampleData$signal_index_a)
  })
  names(ji_samples) <- paste0(methods, "_samples_a")
  
  return(list(features = ji_features, samples = ji_samples))
}

# ------------------------- Compare Pairs -------------------------
ComparePairs <- function(FeatureData, SampleData) {
  methods <- c("fabia","mofa","mfa","gfa")
  
  pairwise <- function(df, type = "features") {
    combn(methods, 2, function(pair) {
      m1 <- paste0(pair[1], "_index_a")
      m2 <- paste0(pair[2], "_index_a")
      jaccard.index.sim(df[[m1]], df[[m2]])
    })
  }
  
  ji_feat_pairs <- pairwise(FeatureData, "features")
  ji_smp_pairs  <- pairwise(SampleData, "samples")
  
  names(ji_feat_pairs) <- combn(methods, 2, paste, collapse = "_vs_")
  names(ji_smp_pairs)  <- combn(methods, 2, paste, collapse = "_vs_")
  
  return(list(feature_pairs = ji_feat_pairs, sample_pairs = ji_smp_pairs))
}

