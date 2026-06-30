simulate_nonoverlapping_signals_varying <- function(n_samples = 100,
                                                    n_features = c(7000, 3000),
                                                    n_blocks = 2,
                                                    sample_range = c(5, 20),
                                                    feature_range = c(100, 500),
                                                    signal_strength = 2,
                                                    noise_sd = 0.2,
                                                    seed = 123) {
  set.seed(seed)
  K <- length(n_features)
  omics <- list()
  signal_info <- list()
  
  for (k in 1:K) {
    n_samp <- n_samples
    n_feat <- n_features[k]
    mat <- matrix(rnorm(n_samp * n_feat, mean = 0, sd = noise_sd), nrow = n_samp, ncol = n_feat)
    
    used_sample_idx <- c()
    used_feature_idx <- c()
    blocks <- list()
    
    for (b in 1:n_blocks) {
      repeat {
        # Random block size
        samp_block <- sample(seq(sample_range[1], sample_range[2]), 1)
        feat_block <- sample(seq(feature_range[1], feature_range[2]), 1)
        
        samp_start <- sample(1:(n_samp - samp_block), 1)
        feat_start <- sample(1:(n_feat - feat_block), 1)
        
        samp_idx <- samp_start:(samp_start + samp_block - 1)
        feat_idx <- feat_start:(feat_start + feat_block - 1)
        
        # Ensure no overlap
        if (!any(samp_idx %in% used_sample_idx) && !any(feat_idx %in% used_feature_idx)) {
          break
        }
      }
      
      mat[samp_idx, feat_idx] <- mat[samp_idx, feat_idx] + signal_strength
      
      used_sample_idx <- c(used_sample_idx, samp_idx)
      used_feature_idx <- c(used_feature_idx, feat_idx)
      
      blocks[[b]] <- list(samples = samp_idx, features = feat_idx,
                          sample_size = samp_block, feature_size = feat_block)
    }
    
    omics[[k]] <- mat
    signal_info[[k]] <- blocks
  }
  
  names(omics) <- paste0("omic", seq_len(K))
  return(list(omics = omics, signals = signal_info))
}


sim <- simulate_nonoverlapping_signals_varying(
  n_samples = 100,
  n_features = c(7000, 3000),
  n_blocks = 3,
  sample_range = c(5, 20),
  feature_range = c(200, 500),
  signal_strength = 3
)

# Visualize OMIC2 with signal patterns
library(ggplot2)
library(reshape2)

omic2 <- sim$omics$omic2
omic2_melt <- melt(omic2)
colnames(omic2_melt) <- c("Sample", "Feature", "Value")

ggplot(omic2_melt, aes(x = Sample, y = Feature, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "#fddc8a", high = "#622400") +
  theme_minimal() +
  labs(title = "OMIC 2: Varying Non-overlapping Signal Blocks",
       x = "samples", y = "features")

# Merge omic1 and omic2 by stacking features (rows)
dataset <- t(cbind(sim$omics$omic1, sim$omics$omic2))

# Basic heatmap of merged dataset
set.seed(657)
image(
  x = 1:dim(dataset)[2],  # samples on x-axis
  y = 1:dim(dataset)[1],  # features on y-axis
  z = t(dataset[nrow(dataset):1, ]),  # flip vertically for correct orientation
  col = heat.colors(100),
  xlab = "Samples",
  ylab = "Features",
  font.lab = 2
)

# Add a dashed horizontal line to separate OMIC 1 (4000 features)
abline(h = 7000, col = "yellow", lty = "dashed", lwd = 2)
