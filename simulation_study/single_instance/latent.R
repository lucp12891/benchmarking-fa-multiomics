set.seed(123)  # For reproducibility

# Load required package
library(MASS)  # For generating multivariate normal data

# Parameters
n_samples <- 100   # Number of samples
n_features_omics1 <- 50  # Features in first omics dataset (e.g., transcriptomics)
n_features_omics2 <- 50  # Features in second omics dataset (e.g., metabolomics)
n_latent <- 5   # Number of shared latent factors

# 1. Generate shared latent factors (representing biological signals)
latent_factors <- matrix(rnorm(n_samples * n_latent, mean=0, sd=1), nrow=n_samples, ncol=n_latent)

# 2. Create loading matrices (effect of latent factors on features)
loadings_omics1 <- matrix(rnorm(n_features_omics1 * n_latent, mean=0.5, sd=0.2), nrow=n_features_omics1, ncol=n_latent)
loadings_omics2 <- matrix(rnorm(n_features_omics2 * n_latent, mean=0.5, sd=0.2), nrow=n_features_omics2, ncol=n_latent)

# 3. Generate omics datasets by combining latent factors and loadings
omics1_shared <- latent_factors %*% t(loadings_omics1)
omics2_shared <- latent_factors %*% t(loadings_omics2)

# 4. Add dataset-specific variation (e.g., technical noise, unique effects)
omics1_specific <- matrix(rnorm(n_samples * n_features_omics1, mean=0, sd=0.3), n_samples, n_features_omics1)
omics2_specific <- matrix(rnorm(n_samples * n_features_omics2, mean=0, sd=0.3), n_samples, n_features_omics2)

# 5. Final simulated datasets
omics1 <- omics1_shared + omics1_specific  # Simulated transcriptomics
omics2 <- omics2_shared + omics2_specific  # Simulated metabolomics

# 6. Store as data frames
omics1_df <- as.data.frame(omics1)
omics2_df <- as.data.frame(omics2)

# Check structure
dim(omics1_df)  # 100 samples x 50 features
dim(omics2_df)  # 100 samples x 40 features

dataset <- as.matrix(t(omics1_df))#dataset <- omic.two[[1]]
image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")
