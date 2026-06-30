# -------------------------- Cutoff Strategies --------------------------

# Cutoff 1: max abs value / sqrt(noise variance)
varphi.one <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    return(max(abs(numeric_values)) / sqrt(sigma))
  }
  return(NA)
}

# Cutoff 2: max abs value / (sigma/2)
varphi.two <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    return(max(abs(numeric_values)) / (sigma/2))
  }
  return(NA)
}

# Cutoff 3: mean of top 100 loadings or top 9 scores / sqrt(sigma)
varphi.three <- function(loading_or_score) {
  vector_name <- deparse(substitute(loading_or_score))
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    x <- abs(numeric_values)
    if (grepl("loading", vector_name)) {
      return(mean(head(sort(x, decreasing = TRUE), 100)) / sqrt(sigma))
    } else if (grepl("score", vector_name)) {
      return(mean(head(sort(x, decreasing = TRUE), 9)) / sqrt(sigma))
    }
  }
  return(NA)
}

# Cutoff 4: mean of top 100 loadings or top 9 scores / (sigma/2)
varphi.four <- function(loading_or_score) {
  vector_name <- deparse(substitute(loading_or_score))
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    x <- abs(numeric_values)
    if (grepl("loading", vector_name)) {
      return(mean(head(sort(x, decreasing = TRUE), 100)) / (sigma/2))
    } else if (grepl("score", vector_name)) {
      return(mean(head(sort(x, decreasing = TRUE), 9)) / (sigma/2))
    }
  }
  return(NA)
}

# Cutoff 5: based on true proportion of signal (quantile thresholding)
varphi.five <- function(loading_or_score, omic = NULL, factor = NULL, type = c("loading", "score")) {
  type <- match.arg(type)
  numeric_values <- as.numeric(unlist(loading_or_score))
  buffer <- 1e-6
  numeric_values <- numeric_values + buffer
  
  indices_features.OMIC1.A <- global_env$indices_features.OMIC1.A 
  indices_features.OMIC1.B <- global_env$indices_features.OMIC1.B  
  indices_features.OMIC2.A <- global_env$indices_features.OMIC2.A
  indices_samples.1A <- global_env$indices_samples.1A  
  indices_samples.2B <- global_env$indices_samples.2B
  n_features_one <- global_env$n_features_one
  n_features_two <- global_env$n_features_two
  n_samples <- global_env$n_samples
  
  if (type == "loading") {
    if (omic == 1 && factor == 1) {
      return(quantile(abs(numeric_values), probs = 1 - (length(indices_features.OMIC1.A) / n_features_one))[[1]])
    } else if (omic == 1 && factor == 2) {
      return(quantile(abs(numeric_values), probs = 1 - (length(indices_features.OMIC1.B) / n_features_one))[[1]])
    } else if (omic == 2 && factor == 1) {
      return(quantile(abs(numeric_values), probs = 1 - (length(indices_features.OMIC2.A) / n_features_two))[[1]])
    } else if (omic == 2 && factor == 2) {
      return(max(abs(numeric_values)))
    }
  }
  
  if (type == "score") {
    if (factor == 1) {
      return(quantile(abs(numeric_values), probs = 1 - (length(indices_samples.1A) / n_samples))[[1]])
    } else if (factor == 2) {
      return(quantile(abs(numeric_values), probs = 1 - (length(indices_samples.2B) / n_samples))[[1]])
    }
  }
  return(NA)
}

# Cutoff 6: half the normalized [0,1] range
normalize_columns <- function(df, columns, clip_percentiles = TRUE) {
  normalize <- function(x, lower_bound = NULL, upper_bound = NULL) {
    if (!is.null(lower_bound) && !is.null(upper_bound)) {
      x <- pmin(pmax(x, lower_bound), upper_bound)
    }
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    if (x_max - x_min == 0) return(rep(0.5, length(x)))
    (x - x_min) / (x_max - x_min)
  }
  
  for (col in columns) {
    if (clip_percentiles) {
      lower_bound <- quantile(df[[col]], 0.01, na.rm = TRUE)
      upper_bound <- quantile(df[[col]], 0.99, na.rm = TRUE)
    } else {
      lower_bound <- upper_bound <- NULL
    }
    df[[col]] <- normalize(df[[col]], lower_bound, upper_bound)
  }
  return(df)
}

varphi.six <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    normalized_values <- normalize_columns(
      data.frame(values = numeric_values), columns = "values"
    )$values
    return(0.5 * (max(normalized_values) - min(normalized_values)))
  }
  return(NA)
}

# Cutoff 7–9: quantiles
varphi.seven  <- function(x) ifelse(all(!is.na(x <- as.numeric(unlist(x)))), quantile(x+1e-6, 0.80)[[1]], NA)
varphi.eight  <- function(x) ifelse(all(!is.na(x <- as.numeric(unlist(x)))), quantile(x+1e-6, 0.85)[[1]], NA)
varphi.nine   <- function(x) ifelse(all(!is.na(x <- as.numeric(unlist(x)))), quantile(x+1e-6, 0.90)[[1]], NA)

# Cutoff 10: k-means clustering (force k = 2 for signal/noise)
varphi.ten <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  numeric_values <- numeric_values[!is.na(numeric_values)]  # drop NA
  
  # Force at least 2 unique values by jittering if necessary
  if (length(unique(numeric_values)) < 2) {
    numeric_values <- jitter(numeric_values, amount = 1e-6)
  }
  
  # Run k-means with centers = 2
  km <- kmeans(numeric_values, centers = 2, nstart = 10)
  
  # Assign the higher-mean cluster as "signal"
  cluster_means <- tapply(numeric_values, km$cluster, mean)
  signal_cluster <- which.max(cluster_means)
  
  # Cutoff = mean of signal cluster
  return(mean(numeric_values[km$cluster == signal_cluster], na.rm = TRUE))
}

# Cutoff 11: Gaussian Mixture Model (force G = 2 for signal/noise)
varphi.eleven <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  numeric_values <- numeric_values[!is.na(numeric_values)]
  
  # Force at least 2 unique values by jittering if necessary
  if (length(unique(numeric_values)) < 2) {
    numeric_values <- jitter(numeric_values, amount = 1e-6)
  }
  
  # Fit Gaussian Mixture Model with G = 2
  gmm <- Mclust(numeric_values, G = 2, verbose = FALSE)
  
  # Pick cluster with higher mean as "signal"
  signal_cluster <- which.max(gmm$parameters$mean)
  
  return(mean(numeric_values[gmm$classification == signal_cluster], na.rm = TRUE))
}

# Cutoff 12: Rolling mean (zoo)
varphi.twelve <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    rolling_mean <- zoo::rollmean(numeric_values, k = 10, fill = NA, align = "center")
    rolling_mean[is.na(rolling_mean)] <- 0
    binary_vector <- ifelse(numeric_values > rolling_mean, 1, 0)
    signal_values <- numeric_values[binary_vector == 1]
    if (length(signal_values) > 0) return(mean(signal_values))
  }
  return(NA)
}

# Collect in a list
varphi_functions <- list(
  varphi.one    = varphi.one,
  varphi.two    = varphi.two,
  varphi.three  = varphi.three,
  varphi.four   = varphi.four,
  varphi.five   = varphi.five,
  varphi.six    = varphi.six,
  varphi.seven  = varphi.seven,
  varphi.eight  = varphi.eight,
  varphi.nine   = varphi.nine,
  varphi.ten    = varphi.ten,
  varphi.eleven = varphi.eleven,
  varphi.twelve = varphi.twelve
)
