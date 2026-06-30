# Load required packages
library(ggplot2)
library(reshape2)

# Define the simulation function (as discussed earlier)
simulate_multi_omics_pmf <- function(n_samples = 100, 
                                     omic_features = c(50, 60, 40),
                                     r_shared = 5, 
                                     r_specific = 3, 
                                     noise_sd = 0.1, 
                                     missing_prob = 0.2, 
                                     seed = 123) {
  set.seed(seed)
  K <- length(omic_features)
  U_shared <- matrix(rnorm(n_samples * r_shared), n_samples, r_shared)
  omics <- list()
  
  for (k in 1:K) {
    p_k <- omic_features[k]
    V_k <- matrix(rnorm(p_k * r_shared), p_k, r_shared)
    W_k <- matrix(rnorm(p_k * r_specific), p_k, r_specific)
    U_k <- matrix(rnorm(n_samples * r_specific), n_samples, r_specific)
    
    X_k <- U_shared %*% t(V_k) + U_k %*% t(W_k)
    X_k <- X_k + matrix(rnorm(n_samples * p_k, sd = noise_sd), n_samples, p_k)
    
    mask <- matrix(runif(n_samples * p_k) > missing_prob, n_samples, p_k)
    X_k[!mask] <- NA
    
    omics[[k]] <- X_k
  }
  
  names(omics) <- paste0("omic", seq_len(K))
  return(list(omics = omics, U_shared = U_shared, V_k=V_k))
}

# 🔄 Run the simulation
sim_data <- simulate_multi_omics_pmf(n_samples = 80, 
                                     omic_features = c(30, 40), 
                                     r_shared = 4, 
                                     r_specific = 2, 
                                     noise_sd = 0.05, 
                                     missing_prob = 0.25)

# 📦 View the structure
str(sim_data)

# 👁️ Preview omic1
head(sim_data$omics$omic1)

# 📊 Visualize the missing pattern for omic1
X1 <- sim_data$omics$omic1
X1_df <- melt(X1)
names(X1_df) <- c("Sample", "Feature", "Value")
X1_df$Missing <- is.na(X1_df$Value)

ggplot(X1_df, aes(x = Feature, y = Sample, fill = Missing)) +
  geom_tile() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white")) +
  labs(title = "Missing Value Pattern in Omic1", fill = "Missing") +
  theme_minimal()
