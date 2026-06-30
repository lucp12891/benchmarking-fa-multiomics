
# Set the working directory
#setwd("C:/Users/Lenovo/Documents/PhD Integromics/Simulation Study/R Code/test_sim")
setwd("C:/Users/bosangir/OneDrive - Studiecentrum voor Kernenergie/Desktop/Ausan - DEL/Output - Paper II")
#=------------------------------------------------------------------------------------

# Generate sigma and corresponding SNR values for OMIC Two, Samples, and OMIC One
sigma <- seq(2, 30, by = 1)
snr_omic_two <- 5.5 / sigma
snr_samples <- 3 / sigma
snr_omic_one <- 5 / sigma

snr_shifted_samples <- 4.5 / sigma

# Set the working directory
setwd("C:/Users/Lenovo/Documents/PhD Integromics/Simulation Study/R Code/test_sim")

# Save plot as EPS
setEPS()
postscript("snr_plot.eps", width = 7, height = 5, paper = "special", horizontal = FALSE)

# Set up the plot
plot(sigma, snr_samples, type = "o", col = "orange", pch = 16, lty = 1,
     ylim = c(0, max(snr_omic_two)), xlab = expression("Noise Level (" * sigma * ")"),
     ylab = "SNR")

# Add the second line
lines(sigma, snr_omic_one, type = "o", col = "green", pch = 16, lty = 1, cex = 1.0)

# Add the third line
lines(sigma, snr_omic_two, type = "o", col = "blue", pch = 16, lty = 1, cex = 1.0)

# Add a legend
legend("topright", legend = c("Samples (μ = 3)", "Omic One (μ = 5)", "Omic Two (μ = 5.5)"),
       col = c("orange", "green", "blue"), pch = 16, lty = 1, bty = "n", cex = 1.0)

# Add formula for SNR
#text(20, max(snr_omic_two) * 0.57, labels = expression(SNR == frac(mu[s], sigma)), cex = 1.1)

# Close the EPS device
dev.off()


# Generate sigma and corresponding SNR values for OMIC Two, Samples, OMIC One, and Shifted Samples
sigma <- seq(2, 30, by = 1)
snr_omic_two <- 5.5 / sigma
snr_samples <- 3 / sigma
snr_omic_one <- 5 / sigma
snr_shifted_samples <- 4.5 / sigma

# Set the working directory
setwd("C:/Users/Lenovo/Documents/PhD Integromics/Simulation Study/R Code/test_sim")

# Save plot as EPS
setEPS()
postscript("snr_plot_with_shift.eps", width = 7, height = 5, paper = "special", horizontal = FALSE)

# Set up the plot
plot(sigma, snr_samples, type = "o", col = "orange", pch = 16, lty = 1,
     ylim = c(0, max(snr_omic_two)), xlab = expression("Noise Level (" * sigma * ")"),
     ylab = "Signal-to-Noise Ratio (SNR)")

# Add the second line for Omic One
lines(sigma, snr_omic_one, type = "o", col = "green", pch = 16, lty = 1, cex = 1.0)

# Add the third line for Omic Two
lines(sigma, snr_omic_two, type = "o", col = "blue", pch = 16, lty = 1, cex = 1.0)

# Add the shifted samples line
lines(sigma, snr_shifted_samples, type = "o", col = "red", pch = 16, lty = 2, cex = 1.0)

# Add arrows demonstrating the shift
arrows(sigma, snr_samples, sigma, snr_shifted_samples, length = 0.1, angle = 15, col = "purple")

# Add a legend
legend("topright", 
       legend = c("Samples (μ = 3)", "Omic One (μ = 5)", "Omic Two (μ = 5.5)", "Shifted Samples (μ = 4.5)"),
       col = c("orange", "green", "blue", "red"), pch = 16, lty = c(1, 1, 1, 2), bty = "n", cex = 0.9)

# Add formula for SNR
text(24, max(snr_omic_two) * 0.75, labels = expression(SNR == frac(mu[s], sigma)), cex = 0.8)

# Close the EPS device
dev.off()


# Set the working directory
setwd("C:/Users/Lenovo/Documents/PhD Integromics/Simulation Study/R Code/test_sim")

# Save plot as EPS
setEPS()
postscript("snr_plot_with_shift_arrows.eps", width = 7, height = 5, paper = "special", horizontal = FALSE)

# Set up the plot
plot(sigma, snr_samples, type = "o", col = "orange", pch = 16, lty = 1,
     ylim = c(0, max(snr_omic_two)), xlab = expression("Noise Level (" * sigma * ")"),
     ylab = "Signal-to-Noise Ratio (SNR)")

# Add the second line for Omic One
lines(sigma, snr_omic_one, type = "o", col = "green", pch = 16, lty = 1, cex = 1.0)

# Add the third line for Omic Two
lines(sigma, snr_omic_two, type = "o", col = "blue", pch = 16, lty = 1, cex = 1.0)

# Add the shifted samples line
lines(sigma, snr_shifted_samples, type = "o", col = "red", pch = 16, lty = 2, cex = 1.0)

# Specify arrow positions
arrow_positions <- c(2, 4, 6, 9, 12, 16, 20, 25, 30)
arrow_indices <- match(arrow_positions, sigma)  # Get indices for specified sigma values

# Add arrows demonstrating the shift
for (i in arrow_indices) {
  arrows(sigma[i], snr_samples[i], sigma[i], snr_shifted_samples[i], 
         length = 0.1, angle = 15, col = "purple")
}
#sigma, snr_samples, sigma, snr_shifted_samples, length = 0.1, angle = 15, col = "purple"
# Add a legend
legend("topright", 
       legend = c("Samples (μ = 3)", "Omic One (μ = 5)", "Omic Two (μ = 5.5)", "Shifted Samples (μ = 4.5)"),
       col = c("orange", "green", "blue", "red"), pch = 16, lty = c(1, 1, 1, 2), bty = "n", cex = 0.8)

# Add formula for SNR
text(24, max(snr_omic_two) * 0.67, labels = expression(SNR == frac(mu, sigma)), cex = 0.8)

# Close the EPS device
dev.off()


# Generate sigma and SNR values
sigma <- seq(2, 30, by = 1)

snr_samples <- 3 / sigma
snr_omic_one <- 5 / sigma
snr_omic_two <- 5.5 / sigma
snr_shifted_samples <- 4.5 / sigma

# Plot directly in RStudio Plot Viewer
plot(
  sigma, snr_samples,
  type = "o",
  col = "orange",
  pch = 16,
  lty = 1,
  ylim = c(0, max(snr_omic_two)),
  xlab = expression("Noise Level (" * sigma^2 * ")"),
  ylab = "Signal-to-Noise Ratio (SNR)"
)

lines(sigma, snr_omic_one, type = "o", col = "green", pch = 16, lty = 1)
lines(sigma, snr_omic_two, type = "o", col = "blue", pch = 16, lty = 1)
lines(sigma, snr_shifted_samples, type = "o", col = "red", pch = 16, lty = 2)

# Add arrows only at selected sigma values
arrow_positions <- c(2, 4, 6, 9, 12, 16, 20, 25, 30)
arrow_indices <- match(arrow_positions, sigma)

for (i in arrow_indices) {
  arrows(
    sigma[i], snr_samples[i],
    sigma[i], snr_shifted_samples[i],
    length = 0.1,
    angle = 15,
    col = "purple"
  )
}

legend(
  "topright",
  legend = c(
    "Samples (mu = 3)",
    "Omic One (mu = 5)",
    "Omic Two (mu = 5.5)",
    "Shifted Samples (mu = 4.5)"
  ),
  col = c("orange", "green", "blue", "red"),
  pch = 16,
  lty = c(1, 1, 1, 2),
  bty = "n",
  cex = 0.8
)

text(
  24,
  max(snr_omic_two) * 0.67,
  labels = expression(SNR == frac(mu, sigma)),
  cex = 0.8
)

# Generate sigma and SNR values
sigma <- seq(2, 30, by = 1)

snr_samples  <- 3 / sigma
snr_omic_one <- 5 / sigma
snr_omic_two <- 5.5 / sigma

# Plot in RStudio viewer
plot(
  sigma, snr_samples,
  type = "o",
  col = "orange",
  pch = 16,
  lty = 1,
  ylim = c(0, max(snr_omic_two)),
  xlab = expression("Noise Level (" * sigma^2 * ")"),
  ylab = "Signal-to-Noise Ratio (SNR)"
)

lines(sigma, snr_omic_one, type = "o", col = "green", pch = 16, lty = 1)
lines(sigma, snr_omic_two, type = "o", col = "blue", pch = 16, lty = 1)

legend(
  "topright",
  legend = c("Samples (mu = 3)", "Omic One (mu = 5)", "Omic Two (mu = 5.5)"),
  col = c("orange", "green", "blue"),
  pch = 16,
  lty = 1,
  bty = "n",
  cex = 1.0
)

# Optional formula annotation
text(
  20,
  max(snr_omic_two) * 0.6,
  labels = expression(SNR == frac(mu, sigma)),
  cex = 0.9
)