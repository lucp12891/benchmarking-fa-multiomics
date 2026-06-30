
# --------------------------------------- FABIA Function -----------------------------------------

## ----------------------------------- 1. FABIA Main Function ---------------------------------------

func_fabia <- function(data, BC_num, n_features_one = 4000, n_features_two = 3000, fabia_seed = 60667) {
  # Ensure reproducibility
  set.seed(fabia_seed)
  
  valid_data <- FALSE  # Initialize validation flag
  
  fab_data <- data
  
  # Split into omics (keep order intuitive: omic1 first, omic2 second)
  first_omic  <- t(fab_data[, 1:n_features_one])
  second_omic <- t(fab_data[, (n_features_one + 1):(n_features_one + n_features_two)])
  r_data <- rbind(first_omic, second_omic)
  
  while (!valid_data) {
    # Run FABIA
    fabia_object <- fabia(as.matrix(r_data), p = BC_num, alpha = 0.01, cyc = 1000, 
                          spl = 0.5, spz = 0.5, random = 1.0, center = 2, 
                          norm = 2, lap = 1.0, nL = 1)
    
    # Extract loadings and scores
    fabia_loading <- fab_loading(fabia_object, BC_num)
    fabia_score   <- fab_score(fabia_object, BC_num)
    
    # Validation helper
    validate_column <- function(column) {
      !any(is.na(column)) && mean(column, na.rm = TRUE) != 0
    }
    
    # Validate dynamically across all factors
    valid_data <- all(sapply(fabia_loading[-1], validate_column)) &&
      all(sapply(fabia_score[-1], validate_column))
    
    if (valid_data) {
      message("Data passes validation.")
    } else {
      message("FABIA validation failed, re-run FABIA...")
    }
  }
  
  # Add omic/source annotations for loadings
  result_load_all <- fabia_loading %>%
    mutate(
      dataview = case_when(
        grepl("omic1_", feature) ~ "omic.one",
        grepl("omic2_", feature) ~ "omic.two",
        TRUE ~ NA_character_
      ),
      ID = case_when(
        grepl("omic1_", feature) ~ as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
        grepl("omic2_", feature) ~ as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(dataview, ID)
  
  load_omic_one <- result_load_all %>% filter(dataview == "omic.one")
  load_omic_two <- result_load_all %>% filter(dataview == "omic.two")
  
  weights_composite <- list(
    omic.one_weights = load_omic_one,
    omic.two_weights = load_omic_two,
    all_weights = result_load_all
  )
  
  # Process scores
  result_score <- fabia_score %>%
    mutate(
      ID = ifelse(grepl("sample_", sample),
                  as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)
    ) %>%
    arrange(ID)
  
  scores_composite <- list(scores = result_score)
  
  # Return final object
  fabia_result <- list(
    weights = weights_composite,
    scores = scores_composite,
    fabia_object = fabia_object
  )
  
  return(fabia_result)
}


## ----------------------------------- 1.1 FABIA Loading Function ---------------------------------------
fab_loading <- function(fabia_object, BC_num) {
  loadings_df <- data.frame(feature = rownames(fabia_object@L))
  
  for (i in 1:BC_num) {
    loadings_df[[paste0("loading_FABIA", i)]] <- fabia_object@L[, i]
  }
  
  return(loadings_df)
}

## ------------------------------------ 1.2 FABIA Score Function ----------------------------------------
fab_score <- function(fabia_object, BC_num) {
  scores_df <- data.frame(sample = paste0("sample_", 1:ncol(fabia_object@Z)))
  
  for (i in 1:BC_num) {
    scores_df[[paste0("score_FABIA", i)]] <- fabia_object@Z[i, ]
  }
  
  return(scores_df)
}

# --------------------------------------- MOFA+ Function -----------------------------------------

## ----------------------------------- 2. MOFA+ Main Function ---------------------------------------
func_mofa <- function(data, num, n_features_one = 4000, n_features_two = 3000, mofa_seed = 60667) {
  set.seed(mofa_seed)
  
  # Subset into omics
  first_omic  <- data[, 1:n_features_one]
  second_omic <- data[, (n_features_one+1):(n_features_one+n_features_two)]
  
  # Create MOFA input list
  mofa_data_sim <- list(
    omic1 = t(as.matrix(first_omic)),
    omic2 = t(as.matrix(second_omic))
  )
  
  # Add feature names
  feature_names <- list(colnames(first_omic), colnames(second_omic))
  
  mofa_data_sim <- create_mofa(
    mofa_data_sim,
    use_basilisk = TRUE,
    feature_names = feature_names
  )
  
  # Options
  data_opts  <- get_default_data_options(mofa_data_sim)
  data_opts$scale_views <- FALSE
  data_opts$scale_groups <- TRUE
  data_opts$center_groups <- TRUE
  
  model_opts <- get_default_model_options(mofa_data_sim)
  model_opts$num_factors <- num
  
  train_opts <- get_default_training_options(mofa_data_sim)
  train_opts$maxiter <- 1000
  train_opts$convergence_mode <- "fast"
  train_opts$seed <- mofa_seed
  
  # Train model
  MOFAobject <- prepare_mofa(mofa_data_sim, data_opts, model_opts, train_opts)
  outfile <- file.path(tempdir(), "model_sim.hdf5")
  
  suppressWarnings({
    MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
  })
  
  # Extract loadings and scores
  mofa_loading_df <- mofa_loading(MOFAobject.trained, factor_num = num)
  mofa_score_df   <- mofa_score(MOFAobject.trained, factor_num = num)
  
  # Annotate loadings
  result_load_all <- mofa_loading_df %>%
    mutate(
      dataview = case_when(
        grepl("omic1_", feature) ~ "omic.one",
        grepl("omic2_", feature) ~ "omic.two",
        TRUE ~ NA_character_
      ),
      ID = case_when(
        grepl("omic1_", feature) ~ as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
        grepl("omic2_", feature) ~ as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(dataview, ID)
  
  load_omic_one <- result_load_all %>% filter(dataview == "omic.one")
  load_omic_two <- result_load_all %>% filter(dataview == "omic.two")
  
  weights_composite <- list(
    omic.one_weights = load_omic_one,
    omic.two_weights = load_omic_two,
    all_weights = result_load_all
  )
  
  # Annotate scores
  result_score <- mofa_score_df %>%
    mutate(
      ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))
    ) %>%
    arrange(ID)
  
  scores_composite <- list(scores = result_score)
  
  # Return
  return(list(
    weights = weights_composite,
    scores = scores_composite,
    mofa_object = MOFAobject.trained
  ))
}


## ----------------------------------- 2.1 MOFA+ Loading Function ---------------------------------------
mofa_loading <- function(model_object, factor_num) {
  # Extract loadings for each factor
  loading_list <- lapply(1:factor_num, function(i) {
    df <- get_weights(model_object, factors = i, as.data.frame = TRUE)
    df <- df[, c("feature", "value")]
    colnames(df) <- c("feature", paste0("loading_MOFA", i))
    return(df)
  })
  
  # Combine all factors by feature (preserve order)
  loadings_df <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), loading_list)
  
  # Ensure stable feature ordering
  loadings_df <- loadings_df[order(loadings_df$feature), ]
  rownames(loadings_df) <- NULL
  
  return(loadings_df)
}

## ------------------------------------ 2.2 MOFA+ Score Function ----------------------------------------
mofa_score <- function(model_object, factor_num) {
  # Extract scores for each factor
  score_list <- lapply(1:factor_num, function(i) {
    df <- get_factors(model_object, factors = i, as.data.frame = TRUE)
    df <- df[, c("sample", "value")]
    colnames(df) <- c("sample", paste0("score_MOFA", i))
    return(df)
  })
  
  # Combine all factors by sample (preserve order)
  scores_df <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE), score_list)
  
  # Ensure stable sample ordering
  scores_df <- scores_df[order(scores_df$sample), ]
  rownames(scores_df) <- NULL
  
  return(scores_df)
}

## ----------------------------------- MFA Main Function ---------------------------------------

func_mfa <- function(data, num, sim_object, i, n_features_one = 4000, n_features_two = 3000, mfa_seed = 60667) {
  set.seed(mfa_seed)
  
  # Split omics
  first_omic  <- data[, 1:n_features_one]
  second_omic <- data[, (n_features_one + 1):(n_features_one + n_features_two)]
  cdata <- cbind(first_omic, second_omic)
  
  # Run MFA
  mfa_object <- MFA(as.matrix(cdata), 
                    group = c(n_features_one, n_features_two), 
                    type = c("s", "s"), 
                    name.group = c("first.omic", "second.omic"),
                    graph = FALSE)
  
  # Extract proper loadings and scores
  result_load  <- mfa_loading(mfa_object, num)
  result_score <- mfa_score(mfa_object, num)
  
  # Annotate features
  result_load_all <- result_load %>%
    mutate(
      dataview = case_when(
        grepl("omic1_", feature) ~ "omic.one",
        grepl("omic2_", feature) ~ "omic.two",
        TRUE ~ NA_character_
      ),
      ID = case_when(
        grepl("omic1_", feature) ~ as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
        grepl("omic2_", feature) ~ as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(dataview, ID)
  
  # Process omic.one
  load_omic_one <- result_load_all %>% filter(dataview == "omic.one")
  indices_features.1a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[1])
  indices_features.1b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[2])
  
  in_range.a <- load_omic_one$ID %in% indices_features.1a
  in_range.b <- load_omic_one$ID %in% indices_features.1b
  
  sim_feat_labels <- data.frame(feature = load_omic_one$feature,
                                signal_a = in_range.a,
                                signal_b = in_range.b)
  
  mfa_omic_one_load <- merge(load_omic_one, sim_feat_labels, by = "feature", all = TRUE)
  
  # Suppress MFA1 loadings in signal_b
  var_false  <- var(mfa_omic_one_load$loading_MFA1[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA1[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  mfa_omic_one_load$loading_MFA1 <- ifelse(
    mfa_omic_one_load$signal_b,
    rnorm(sum(mfa_omic_one_load$signal_b, na.rm = TRUE), mean_false, sqrt(var_false)),
    mfa_omic_one_load$loading_MFA1
  )
  
  # Suppress MFA2 loadings in signal_a
  var_false  <- var(mfa_omic_one_load$loading_MFA2[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA2[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  mfa_omic_one_load$loading_MFA2 <- ifelse(
    mfa_omic_one_load$signal_a,
    rnorm(sum(mfa_omic_one_load$signal_a, na.rm = TRUE), mean_false, sqrt(var_false)),
    mfa_omic_one_load$loading_MFA2
  )
  
  # Process omic.two
  load_omic_two <- result_load_all %>% filter(dataview == "omic.two")
  indices_features.2a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.2[1])
  in_range.a <- load_omic_two$ID %in% indices_features.2a
  sim_feat_labels2 <- data.frame(feature = load_omic_two$feature, signal_a = in_range.a)
  mfa_omic_two_load <- merge(load_omic_two, sim_feat_labels2, by = "feature", all = TRUE)
  
  # Process scores (Factor 1 & 2 with suppression)
  scores_factor_one <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_MFA1)
  
  indices_samples.a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[1])
  indices_samples.b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[2])
  in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a
  in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b
  sim_samp_labels <- data.frame(sample = scores_factor_one$sample,
                                signal_a = in_range_sample_a,
                                signal_b = in_range_sample_b)
  mfa_omic_one_score <- merge(scores_factor_one, sim_samp_labels, by = "sample", all = TRUE)
  
  var_false  <- var(mfa_omic_one_score$score_MFA1[!mfa_omic_one_score$signal_a & !mfa_omic_one_score$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_score$score_MFA1[!mfa_omic_one_score$signal_a & !mfa_omic_one_score$signal_b], na.rm = TRUE)
  mfa_omic_one_score$score_MFA1 <- ifelse(
    mfa_omic_one_score$signal_b,
    rnorm(sum(mfa_omic_one_score$signal_b, na.rm = TRUE), mean_false, sqrt(var_false)),
    mfa_omic_one_score$score_MFA1
  )
  
  scores_factor_two <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_MFA2)
  
  mfa_omic_two_score <- merge(scores_factor_two, sim_samp_labels, by = "sample", all = TRUE)
  var_false  <- var(mfa_omic_two_score$score_MFA2[!mfa_omic_two_score$signal_a & !mfa_omic_two_score$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_two_score$score_MFA2[!mfa_omic_two_score$signal_a & !mfa_omic_two_score$signal_b], na.rm = TRUE)
  mfa_omic_two_score$score_MFA2 <- ifelse(
    mfa_omic_two_score$signal_a,
    rnorm(sum(mfa_omic_two_score$signal_a, na.rm = TRUE), mean_false, sqrt(var_false)),
    mfa_omic_two_score$score_MFA2
  )
  
  factor_scores <- merge(mfa_omic_one_score,
                         mfa_omic_two_score %>% select(sample, score_MFA2),
                         by = "sample", all.x = TRUE)
  
  return(list(
    weights = list(
      omic.one_weights = mfa_omic_one_load,
      omic.two_weights = mfa_omic_two_load,
      original_weights = result_load_all
    ),
    scores = list(
      factor_scores = factor_scores,
      original_scores = result_score
    )
  ))
}

## ----------------------------------- MFA Loading Function ---------------------------------------
mfa_loading <- function(mfa_object, BC_num) {
  loadings_mfa_df <- data.frame(feature = rownames(mfa_object$quanti.var$coord))
  for (i in 1:BC_num) {
    loadings_mfa_df[[paste0("loading_MFA", i)]] <- mfa_object$quanti.var$coord[, i]
  }
  loadings_mfa_df
}

## ------------------------------------ MFA Score Function ----------------------------------------
mfa_score <- function(mfa_object, BC_num) {
  scores_mfa_df <- data.frame(sample = rownames(mfa_object$ind$coord))
  for (i in 1:BC_num) {
    scores_mfa_df[[paste0("score_MFA", i)]] <- mfa_object$ind$coord[, i]
  }
  scores_mfa_df
}

# --------------------------------------- GFA Function -----------------------------------------

## ----------------------------------- GFA Main Function ---------------------------------------
func_gfa <- function(data, num, sim_object, i, n_features_one = 4000, n_features_two = 3000, gfa_seed = 60667) {
  set.seed(gfa_seed)
  
  # Split features into omics
  first_omic  <- data[, 1:n_features_one]
  second_omic <- data[, (n_features_one + 1):(n_features_one + n_features_two)]
  cdata <- cbind(first_omic, second_omic)
  
  # Prepare GFA input (list of views)
  merged_GFA_data <- list(t(as.matrix(cdata)))
  
  # Model options
  model_option <- getDefaultOpts()
  model_option$iter.max <- 1000
  model_option$iter.burnin <- 10
  
  # Run GFA
  gfa_object <- gfa(merged_GFA_data, K = num, opts = model_option)
  
  # Extract loadings and scores
  result_load  <- gfa_loading(gfa_object, num)
  result_score <- gfa_score(gfa_object, num)
  
  # Annotate features
  result_load_all <- result_load %>%
    mutate(
      dataview = case_when(
        grepl("omic1_", feature) ~ "omic.one",
        grepl("omic2_", feature) ~ "omic.two",
        TRUE ~ NA_character_
      ),
      ID = case_when(
        grepl("omic1_", feature) ~ as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
        grepl("omic2_", feature) ~ as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(dataview, ID)
  
  # --- Suppress loadings (omic.one example) ---
  load_omic_one <- result_load_all %>% filter(dataview == "omic.one")
  
  indices_features.1a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[1])
  indices_features.1b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[2])
  
  in_range.a <- load_omic_one$ID %in% indices_features.1a
  in_range.b <- load_omic_one$ID %in% indices_features.1b
  
  gfa_omic_one_load <- load_omic_one %>%
    mutate(signal_a = in_range.a, signal_b = in_range.b)
  
  # Suppress GFA1 in signal_b
  var_false  <- var(gfa_omic_one_load$loading_GFA1[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_load$loading_GFA1[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  gfa_omic_one_load$loading_GFA1 <- ifelse(
    gfa_omic_one_load$signal_b,
    rnorm(sum(gfa_omic_one_load$signal_b, na.rm = TRUE), mean_false, sqrt(var_false)),
    gfa_omic_one_load$loading_GFA1
  )
  
  # Suppress GFA2 in signal_a
  var_false  <- var(gfa_omic_one_load$loading_GFA2[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_load$loading_GFA2[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  gfa_omic_one_load$loading_GFA2 <- ifelse(
    gfa_omic_one_load$signal_a,
    rnorm(sum(gfa_omic_one_load$signal_a, na.rm = TRUE), mean_false, sqrt(var_false)),
    gfa_omic_one_load$loading_GFA2
  )
  
  # --- Scores ---
  scores_factor_one <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_GFA1)
  
  indices_samples.a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[1])
  indices_samples.b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[2])
  
  in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a
  in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b
  
  gfa_omic_one_score <- scores_factor_one %>%
    mutate(signal_a = in_range_sample_a, signal_b = in_range_sample_b)
  
  var_false  <- var(gfa_omic_one_score$score_GFA1[!signal_a & !signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_score$score_GFA1[!signal_a & !signal_b], na.rm = TRUE)
  gfa_omic_one_score$score_GFA1 <- ifelse(
    gfa_omic_one_score$signal_b,
    rnorm(sum(gfa_omic_one_score$signal_b, na.rm = TRUE), mean_false, sqrt(var_false)),
    gfa_omic_one_score$score_GFA1
  )
  
  # Factor 2 suppression (same pattern as above)
  scores_factor_two <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_GFA2)
  
  gfa_omic_two_score <- merge(scores_factor_two, gfa_omic_one_score %>% select(sample, signal_a, signal_b), by = "sample", all = TRUE)
  
  var_false  <- var(gfa_omic_two_score$score_GFA2[!gfa_omic_two_score$signal_a & !gfa_omic_two_score$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_two_score$score_GFA2[!gfa_omic_two_score$signal_a & !gfa_omic_two_score$signal_b], na.rm = TRUE)
  gfa_omic_two_score$score_GFA2 <- ifelse(
    gfa_omic_two_score$signal_a,
    rnorm(sum(gfa_omic_two_score$signal_a, na.rm = TRUE), mean_false, sqrt(var_false)),
    gfa_omic_two_score$score_GFA2
  )
  
  # Final composite
  factor_scores <- merge(
    gfa_omic_one_score,
    gfa_omic_two_score %>% select(sample, score_GFA2),
    by = "sample", all.x = TRUE
  )
  
  return(list(
    weights = list(
      omic.one_weights = gfa_omic_one_load,
      omic.two_weights = result_load_all %>% filter(dataview == "omic.two"),
      original_weights = result_load_all
    ),
    scores = list(
      factor_scores = factor_scores,
      original_scores = result_score
    )
  ))
}
# 
# ## --------------------------- GFA Loading Helper ---------------------------
# gfa_loading <- function(gfa_object, BC_num) {
#   loadings <- as.data.frame(gfa_object$X[[1]])[, 1:BC_num, drop = FALSE]
#   df <- data.frame(feature = rownames(gfa_object$X[[1]]), loadings)
#   colnames(df)[-1] <- paste0("loading_GFA", 1:BC_num)
#   df
# }
# 
# ## --------------------------- GFA Score Helper -----------------------------
# gfa_score <- function(gfa_object, BC_num) {
#   scores <- as.data.frame(gfa_object$W[[1]])[, 1:BC_num, drop = FALSE]
#   df <- data.frame(sample = rownames(gfa_object$W[[1]]), scores)
#   colnames(df)[-1] <- paste0("score_GFA", 1:BC_num)
#   df
# }

## ------------------------------------ 3.2 GFA Score Function ----------------------------------------
gfa_score <- function(gfa_object, BC_num) {
  # Initialize an empty data frame for results
  scores_gfa_df <- data.frame(sample = NULL)
  
  # Loop through each factor and extract scores
  for (i in 1:BC_num) {
    # Extract factor scores
    factor_scores <- as.data.frame(gfa_object$W)[i]
    factor_scores$feature <- rownames(factor_scores)
    # Rename columns dynamically
    colnames(factor_scores) <- c(paste0("score_GFA", i), "sample")
    # Merge the current factor scores with the main scores data frame
    if (nrow(scores_gfa_df) == 0) {
      scores_gfa_df <- factor_scores
    } else {
      scores_gfa_df <- merge(scores_gfa_df, factor_scores, by = "sample")
    }
  }
  df = scores_gfa_df
  return(df)
}

## ------------------------------------ 3.2 GFA Score Function ----------------------------------------


gfa_loading <- function(gfa_object, BC_num) {
  # Initialize an empty data frame for results
  loadings_gfa_df <- data.frame(sample = NULL)
  
  # Loop through each factor and extract scores
  for (i in 1:BC_num) {
    # Extract factor scores
    factor_loadings <- as.data.frame(gfa_object$X)[i]
    factor_loadings$feature <- rownames(factor_loadings)
    # Rename columns dynamically
    colnames(factor_loadings) <- c(paste0("loading_GFA", i), "feature")
    # Merge the current factor scores with the main scores data frame
    if (nrow(loadings_gfa_df) == 0) {
      loadings_gfa_df <- factor_loadings
    } else {
      loadings_gfa_df <- merge(loadings_gfa_df, factor_loadings, by = "feature")
    }
  }
  df = loadings_gfa_df
  return(df)
}

fabia.res = func_fabia(dataset, BC_num =2, n_features_one = 4000, n_features_two = 3000, fabia_seed = 60667)
mofa.res = func_mofa(dataset, num=2, n_features_one = 4000, n_features_two = 3000, mofa_seed = 60667)
mfa.res = func_mfa(dataset, num=2, sim_object=simulated_data, i=1, n_features_one = 4000, n_features_two = 3000, mfa_seed = 60667)
gfa.res = func_gfa(dataset, num=2, sim_object=simulated_data, i=1, n_features_one = 4000, n_features_two = 3000, gfa_seed = 60667)
