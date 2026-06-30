varphi.five_multisim <- function(loading_or_score, omic = NULL, factor = NULL, type = c("loading", "score")) {
  type <- match.arg(type)
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  if (type == "loading") {
    if (method == 'multiple.factor'){
      if (omic == 1 && factor == 1) {
        x <- (1 - (length(indices_features.OMIC1.A) / n_features_one))
        varphi <- quantile(abs(numeric_values), probs = x)[[1]]
      }else if (omic == 1 && factor == 2) {
        x <- (1 - (length(indices_features.OMIC1.B) / n_features_one))
        varphi <- quantile(abs(numeric_values), probs = x)[[1]]
      }else if (omic == 2 && factor == 1) {
        y <- (1 - (length(indices_features.OMIC2.A) / n_features_two))
        varphi <- quantile(abs(numeric_values), probs = y)[[1]]
      } 
    } else {
      stop("Invalid method specified.")
    }
    return(varphi)
  }
  
  if (type == "score") {
    if (method =='multiple.factor') {
      if(factor == 1){
        y <- (1 - (length(indices_samples.1A) / n_samples))
        varphi <- quantile(abs(numeric_values), probs = y)[[1]]
      } else if(factor == 2){
        y <- (1 - (length(indices_samples.2B) / n_samples))
        varphi <- quantile(abs(numeric_values), probs = y)[[1]]
      }
    } else {
      stop("Invalid method specified.")
    }
    return(varphi)
  }
  
  stop("Invalid type: must be either 'loading' or 'score'")
}
