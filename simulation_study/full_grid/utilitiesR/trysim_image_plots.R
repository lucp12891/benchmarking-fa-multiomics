# Single factor
simulated_data <- single.factor(
  n_features_one = n_features_one,
  n_features_two = n_features_two,
  n_samples = n_samples,
  sigmas = sigma,
  iterations = iter
)

# Shared factor
simulated_data <- shared.factor(
  n_features_one = n_features_one,
  n_features_two = n_features_two,
  n_samples = n_samples,
  sigmas = sigma,
  iterations = iter 
)

# Multiple(single + shared) factors
simulated_data <- multiple_factor(
  n_features_one = 4000, 
  n_features_two = 3000, 
  n_samples = 100, 
  sigmas = 7, 
  n_factors = 2, 
  iterations = 1
)

# Extra cata
dataset <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]

set.seed(657)
image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples", font.lab = 2)
abline(h = 4000, col = "brown", lty = "dashed")  # lty = "dashed" sets the line type
# Add text annotations for the matrices
# Adjust x and y coordinates based on the location of the dark matrices
text(x = 28, y = 2550, labels = "shared factor", col = "black", cex = 0.8, font = 2)
text(x = 60.5, y = 5540, labels = "unique factor", col = "black", cex = 0.80, font = 2)
persp(c(1:dim(t(dataset))[2]), c(1:dim(t(dataset))[1]), dataset, theta = 325, phi = 15, col = "yellow", xlab = "", ylab = "Samples", zlab = " ")
