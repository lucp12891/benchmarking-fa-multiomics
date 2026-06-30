## ------------------------------- Jaccard similarity --------------------------------
jaccard.index.sim <- function(X, Y, type = c("vector", "matrix")) {
  type <- match.arg(type)
  
  # Helper: compute Jaccard for two binary vectors
  jaccard_vec <- function(a, b) {
    if (!all(a %in% c(0,1)) || !all(b %in% c(0,1))) {
      stop("Inputs must be binary (0/1).")
    }
    similar_count <- sum(a == 1 & b == 1)
    sum_a <- sum(a == 1)
    sum_b <- sum(b == 1)
    denom <- sum_a + sum_b - similar_count
    if (denom == 0) return(1) else return(similar_count / denom)
  }
  
  # --- Case 1: vectors ---
  if (is.vector(X) && is.vector(Y)) {
    if (length(X) != length(Y)) stop("Vectors must have the same length.")
    return(jaccard_vec(X, Y))
  }
  
  # --- Case 2: matrices ---
  if (is.matrix(X) && is.matrix(Y)) {
    if (!all(dim(X) == dim(Y))) stop("Matrices must have the same dimensions.")
    
    p <- ncol(X)
    
    if (type == "vector") {
      out <- numeric(p)
      names(out) <- colnames(X) %||% paste0("Col", 1:p)
      for (i in 1:p) {
        out[i] <- jaccard_vec(X[, i], Y[, i])
      }
      return(out)
    }
    
    if (type == "matrix") {
      out <- matrix(0, nrow = p, ncol = p, dimnames = list(colnames(X), colnames(Y)))
      for (i in 1:p) {
        for (j in 1:p) {
          out[i, j] <- jaccard_vec(X[, i], Y[, j])
        }
      }
      return(out)
    }
  }
  
  stop("Inputs must both be vectors or both be matrices.")
}

## ------------------------------- Cosine similarity --------------------------------
cosine.sim <- function(x, y) {
  if (length(x) != length(y)) stop("Vectors must have the same length.")
  num <- sum(x * y)
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (denom == 0) return(NA)  # undefined if one vector is all zeros
  return(num / denom)
}

## ------------------------------- Metrics --------------------------------
# CompareMetrics <- function(data, signal_index, methods) {
#   
#   # Helper function to calculate metrics
#   calculate_metrics <- function(actual, predicted) {
#     # Build full confusion matrix (always 2x2, ensures no missing rows/cols)
#     cm <- table(
#       Actual = factor(actual, levels = c(0,1)),
#       Predicted = factor(predicted, levels = c(0,1))
#     )
#     
#     TP <- cm["1","1"]
#     TN <- cm["0","0"]
#     FP <- cm["0","1"]
#     FN <- cm["1","0"]
#     
#     # Core metrics
#     sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))   # recall
#     specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
#     accuracy    <- (TP + TN) / sum(cm)
#     precision   <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
#     f1_score    <- ifelse(is.na(precision) || is.na(sensitivity) || 
#                             (precision + sensitivity) == 0, NA,
#                           2 * (precision * sensitivity) / (precision + sensitivity))
#     
#     # MCC
#     denom <- sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
#     mcc   <- ifelse(denom == 0, NA, ((TP * TN) - (FP * FN)) / denom)
#     
#     # AUC (only if predictions are not all binary)
#     auc_value <- NA
#     if (length(unique(actual)) > 1) {
#       if (length(unique(predicted)) > 2) {
#         auc_value <- pROC::auc(pROC::roc(actual, predicted))
#       } else {
#         auc_value <- 0.5  # binary predictions → chance level
#       }
#     }
#     
#     # Return as data frame row
#     return(data.frame(
#       TP = TP, FP = FP, TN = TN, FN = FN,
#       Sensitivity = sensitivity,
#       Specificity = specificity,
#       Accuracy = accuracy,
#       Precision = precision,
#       F1 = f1_score,
#       MCC = mcc,
#       AUC = auc_value
#     ))
#   }
#   
#   # --- Compare methods against ground truth ---
#   result_true_list <- lapply(methods, function(m) {
#     calculate_metrics(data[[signal_index]], data[[m]])
#   })
#   names(result_true_list) <- methods
#   
#   # --- Pairwise comparisons between methods ---
#   result_comparison_list <- list()
#   for (i in 1:(length(methods) - 1)) {
#     for (j in (i+1):length(methods)) {
#       name <- paste(methods[i], "vs", methods[j], sep="_")
#       result_comparison_list[[name]] <- calculate_metrics(data[[methods[i]]], data[[methods[j]]])
#     }
#   }
#   
#   # Combine into data frames
#   per_measures_true <- do.call(rbind, lapply(names(result_true_list), function(x) {
#     cbind(Method = x, result_true_list[[x]])
#   }))
#   
#   per_measures_comparison <- do.call(rbind, lapply(names(result_comparison_list), function(x) {
#     cbind(Method_Comparison = x, result_comparison_list[[x]])
#   }))
#   
#   return(list(
#     per_measures_true = per_measures_true,
#     per_measures_comparison = per_measures_comparison
#   ))
# }
# ------------------------- Compare Metrics -------------------------
CompareMetrics <- function(data, signal_index, methods) {
  results <- lapply(methods, function(method) {
    pred <- data[[method]]
    truth <- data[[signal_index]]
    
    TP <- sum(pred == 1 & truth == 1, na.rm = TRUE)
    TN <- sum(pred == 0 & truth == 0, na.rm = TRUE)
    FP <- sum(pred == 1 & truth == 0, na.rm = TRUE)
    FN <- sum(pred == 0 & truth == 1, na.rm = TRUE)
    
    acc  <- (TP + TN) / (TP + TN + FP + FN)
    prec <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    rec  <- ifelse((TP + FN) > 0, TP / (TP + FN), NA) # sensitivity
    spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    f1   <- ifelse(!is.na(prec) & !is.na(rec) & (prec + rec) > 0,
                   2 * prec * rec / (prec + rec), NA)
    mcc  <- ifelse((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) > 0,
                   (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),
                   NA)
    
    return(data.frame(
      method = method,
      Accuracy = acc,
      Precision = prec,
      Recall = rec,
      Specificity = spec,
      F1 = f1,
      MCC = mcc
    ))
  })
  
  do.call(rbind, results)
}
