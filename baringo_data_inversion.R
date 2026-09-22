# ==========================================================
# PROSPECT inversion to retrieve EWT and LMA
#
# Uses full measured reflectance and transmittance spectra
#
# Outputs:
#   Sample
#   Retrieved_EWT
#   Retrieved_LMA
#
# Author: Victor Korir
# ==========================================================

library(prospect)
library(dplyr)
library(signal)  # For Savitzky-Golay filter

# ==========================================================
# Savitzky-Golay smoothing function
# ==========================================================

smooth_spectrum <- function(spectrum, p = 3, n = 11) {
  # p = polynomial order (typically 2-4)
  # n = window size (must be odd, typically 7-21)
  sgolayfilt(spectrum, p = p, n = n)
}

# ==========================================================
# Load measured spectra
# ==========================================================

reflect <- read.table(
  "Data/Baringo_data/Merged_Data/Reflectance_merged.txt",
  header = TRUE
)

transmit <- read.table(
  "Data/Baringo_data/Merged_data/merged_scaled_Transmittance_W.txt",
  header = TRUE
)

# ==========================================================
# Load metadata/results table
# ==========================================================

results <- read.csv(
  "Data/Baringo_data/estimated_N_results.csv"
)

# ==========================================================
# Prepare spectra
# ==========================================================

reflect_l <- reflect[reflect$wavelength >= 750 & reflect$wavelength <= 820 , ]
transmit_l <- transmit[transmit$wavelength >= 750 & transmit$wavelength <= 820, ]

lambda <- reflect_l$wavelength

refl_samples <- reflect_l[, -1]
refl_samples <- refl_samples %>% select(-'S387')

trans_samples <- transmit_l[, -1]
trans_samples <- trans_samples %>% select(-'S387')

# ==========================================================
# Apply Savitzky-Golay smoothing to all spectra
# ==========================================================

# Smoothing parameters
sg_order <- 3    # polynomial order
sg_window <- 11  # window size (must be odd)

# Smooth reflectance
refl_smoothed <- as.data.frame(
  apply(refl_samples, 2, smooth_spectrum, p = sg_order, n = sg_window)
)

# Smooth transmittance
trans_smoothed <- as.data.frame(
  apply(trans_samples, 2, smooth_spectrum, p = sg_order, n = sg_window)
)

# Ensure values stay in valid range [0, 1]
#refl_smoothed <- pmax(pmin(refl_smoothed, 1), 0)
#trans_smoothed <- pmax(pmin(trans_smoothed, 1), 0)

cat("Smoothing applied: SG filter, order =", sg_order, ", window =", sg_window, "\n")

# Quick visual check (first sample)
# plot(lambda, refl_samples[, 1], type = "l", col = "gray", 
#      main = "Original vs Smoothed", xlab = "Wavelength", ylab = "Reflectance")
# lines(lambda, refl_smoothed[, 1], col = "red", lwd = 2)
# legend("topright", c("Original", "Smoothed"), col = c("gray", "red"), lwd = c(1, 2))


# ==========================================================
# Parameters to estimate
# ==========================================================

parms_to_estimate <- c("lma")

# ==========================================================
# Storage vectors
# ==========================================================
sample_names <- colnames(refl_smoothed)
retrieval_results <- data.frame(
  Sample = sample_names,
  Estimated_N = numeric(length(sample_names)),
  Retrieved_EWT = numeric(length(sample_names)),
  Retrieved_LMA = numeric(length(sample_names)),
  stringsAsFactors = FALSE
)

# ==========================================================
# Loop through samples
# ==========================================================

# Use smoothed spectra



for (i in seq_along(sample_names)) {
  
  sample_name <- sample_names[i]
  match_idx <- which(results$Sample == sample_name)
  
  if (length(match_idx) == 0) {
    warning(paste("No match found for sample:", sample_name))
    next
  }
  
  subset_lop <- fit_spectral_data(
    lambda = lambda,
    refl = refl_smoothed[, sample_name], 
    tran = trans_smoothed[, sample_name]
  )
  
  n_struct_value <- results$Estimated_N[match_idx]
  

  options <- set_options_prospect(fun = 'invert_prospect')
  options$spec_prospect <- subset_lop$spec_prospect
  
  # Fix N 
  options$fix_values <- list(n_struct = n_struct_value)
  
  
  res <- invert_prospect(
    refl = subset_lop$refl, 
    tran = subset_lop$tran, 
    prospect_version = 'D',
    parms_to_estimate = c( "lma", 'ewt'),
    options = options
  )
  
  #print(paste("Sample:", sample_name, "N:", n_struct_value, "LMA:", res$lma))
  print(options$init_values)
  print(options$fix_values)
  
  
 # Retrieved_EWT[i] <- res$ewt
  retrieval_results$Retrieved_LMA[i] <- res$lma * 1000
  retrieval_results$Sample[i] <- sample_name
 
}

# ==========================================================
# Retrieving EWT fixing LMA and N
# ==========================================================

reflect<- reflect[reflect$wavelength >= 400 & reflect$wavelength <= 2200 , ]
transmit <- transmit[transmit$wavelength >= 400 & transmit$wavelength <= 2200, ]

lambda <- reflect$wavelength

refl_samples <- reflect[, -1]
refl_samples <- refl_samples %>% select(-'S387')

trans_samples <- transmit[, -1]
trans_samples <- trans_samples %>% select(-'S387')

# ==========================================================
# Apply Savitzky-Golay smoothing to all spectra
# ==========================================================

# Smoothing parameters
sg_order <- 3    # polynomial order
sg_window <- 11  # window size (must be odd)

# Smooth reflectance
refl_smoothed <- as.data.frame(
  apply(refl_samples, 2, smooth_spectrum, p = sg_order, n = sg_window)
)

# Smooth transmittance
trans_smoothed <- as.data.frame(
  apply(trans_samples, 2, smooth_spectrum, p = sg_order, n = sg_window)
)

# Ensure values stay in valid range [0, 1]
#refl_smoothed <- pmax(pmin(refl_smoothed, 1), 0)
#trans_smoothed <- pmax(pmin(trans_smoothed, 1), 0)

cat("Smoothing applied: SG filter, order =", sg_order, ", window =", sg_window, "\n")

for (i in seq_along(sample_names)) {
  
  sample_name <- sample_names[i]
  match_idx <- which(results$Sample == sample_name)
  match_idx1 <- which(retrieval_results$Sample == sample_name)
  if (length(match_idx) == 0) {
    warning(paste("No match found for sample:", sample_name))
    next
  }
  
  n_struct_value <- results$Estimated_N[match_idx]
  lma = retrieval_results$Retrieved_LMA[match_idx1]
  # Use smoothed spectra
  subset_lop <- fit_spectral_data(
    lambda = lambda,
    refl = refl_smoothed[, i], 
    tran = trans_smoothed[, i]
  )
  
  options <- set_options_prospect(fun = 'invert_prospect_opt')
  options$spec_prospect <- subset_lop$spec_prospect
  
  # Fix N 
  options$fix_values <- list(n_struct = n_struct_value, lma = lma/1000)
  
  
  res <- invert_prospect_opt( lambda = lambda,
    refl = subset_lop$refl, s
    tran = subset_lop$tran, 
    prospect_version = 'D',
    parms_to_estimate = c( "lma", 'ewt'),
    options = options
  )
  
  #print(paste("Sample:", sample_name, "N:", n_struct_value, "LMA:", res$lma))
  print(options$init_values)
  print(options$fix_values)
  
  
  # Retrieved_EWT[i] <- res$ewt
  if(retrieval_results$Sample[i] == sample_name){
    retrieval_results$Retrieved_EWT[i] <- res$ewt * 1000
  }
  
  
}




retrieval_results <- results %>%
  left_join(
    retrieval_results %>% select(Sample,Retrieved_EWT, Retrieved_LMA),
    by = "Sample"
  )

write.csv(
  retrieval_results,
  "Data/Baringo_data/EWT_LMA_inversion_results_LMA.csv",
  row.names = FALSE
)

# ==========================================================
# Summary
# ==========================================================

summary(results$Retrieved_EWT)
summary(results$Retrieved_LMA)

hist(results$Retrieved_EWT, main = "Retrieved EWT", xlab = "EWT")
hist(results$Retrieved_LMA, main = "Retrieved LMA", xlab = "LMA")


# ==========================================================
# Regression plotting function
# ==========================================================

plot_trait_regression <- function(data, measured_col, retrieved_col, 
                                  species_col = "Species",
                                  trait_name = "Trait", units = "") {
  
  library(ggplot2)
  
  # Remove NA values
  plot_data <- data[!is.na(data[[measured_col]]) & !is.na(data[[retrieved_col]]), ]
  
  # Calculate statistics
  r2 <- cor(plot_data[[measured_col]], plot_data[[retrieved_col]])^2
  rmse <- sqrt(mean((plot_data[[measured_col]] - plot_data[[retrieved_col]])^2))
  bias <- mean(plot_data[[retrieved_col]] - plot_data[[measured_col]])
  n_samples <- nrow(plot_data)
  
  # Axis range
  axis_min <- min(c(plot_data[[measured_col]], plot_data[[retrieved_col]]))
  axis_max <- max(c(plot_data[[measured_col]], plot_data[[retrieved_col]]))
  padding <- 0.05 * (axis_max - axis_min)
  axis_range <- c(axis_min - padding, axis_max + padding)
  
  # Stats label
  stats_label <- sprintf("R² = %.3f\nRMSE = %.3f\nBias = %.3f\nn = %d", 
                         r2, rmse, bias, n_samples)
  
  # Plot
  p <- ggplot(plot_data, aes(x = .data[[measured_col]], 
                             y = .data[[retrieved_col]], 
                             color = .data[[species_col]])) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, 
                color = "blue", fill = "lightblue", alpha = 0.2) +
    coord_fixed(ratio = 1, xlim = axis_range, ylim = axis_range) +
    labs(
      title = paste("Retrieved vs Measured", trait_name),
      x = paste("Measured", trait_name, units),
      y = paste("Retrieved", trait_name, units),
      color = "Species"
    ) +
    annotate("text", 
             x = axis_min, 
             y = axis_max,
             label = stats_label, 
             hjust = 0, vjust = 1, 
             size = 4, fontface = "bold") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}

# ==========================================================
# Generate plots
# ==========================================================

# Adjust column names to match your data
p_ewt <- plot_trait_regression(
  results,
  measured_col = "Measured_EWT",
  retrieved_col = "Retrieved_EWT",
  species_col = "Species",
  trait_name = "EWT",
  units = "(mg/cm²)"
)

p_lma <- plot_trait_regression(
  results,
  measured_col = "Measured_LMA",
  retrieved_col = "Retrieved_LMA",
  species_col = "Species",
  trait_name = "LMA",
  units = "(mg/cm²)"
)

print(p_ewt)
print(p_lma)

