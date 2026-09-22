# ==========================================================
# Prepare data for N structural parameter optimization
#
# Includes:
# 1. Savitzky-Golay smoothing (order 2)
# 2. Reflectance metrics
# 3. Transmittance metrics
# 4. Absorptance metrics
# 5. Retrieval of R and T at selected wavelengths
# 6. Matching pigments using sample IDs
#
# Author: Victor Korir
# Refactored for Type-Safety and Speed (No NA Interpolation)
# ==========================================================

library(dplyr)
library(signal)

# ==========================================================
# Load data
# ==========================================================

reflect <- read.table(
  "Data/Baringo_data/Merged_data/Reflectance_merged.txt", 
  sep = '\t', header = TRUE
)

transmit <- read.table(
  "Data/Baringo_data/Merged_data/merged_scaled_Transmittance_T.txt", 
  sep = '\t', header = TRUE
)

pigment <- read.table(
  "Data/Baringo_data/Merged_output_with_species.csv",
  sep = ",", header = TRUE
)

pigment <- pigment %>%
  mutate(
    LMA.mg.cm2.       = as.numeric(LMA.mg.cm2.),
    Cw.mg.cm2.        = as.numeric(Cw.mg.cm2.),
    CHLab_tot_ug_cm2  = as.numeric(CHLab_tot_ug_cm2)
  )

# ==========================================================
# Restrict wavelength range & Extract
# ==========================================================

reflect_750_900 <- reflect[reflect$wavelength >= 750 & reflect$wavelength <= 900, ]
transmit_750_900 <- transmit[transmit$wavelength >= 750 & transmit$wavelength <= 900, ]

wavelengths <- reflect_750_900$wavelength

refl_samples <- reflect_750_900[, -1]
trans_samples <- transmit_750_900[, -1]

# ==========================================================
# SAFE NA FILTERING
# ==========================================================

# Remove columns completely filled with NA
refl_samples <- refl_samples[, colSums(is.na(refl_samples)) < nrow(refl_samples), drop = FALSE]
trans_samples <- trans_samples[, colSums(is.na(trans_samples)) < nrow(trans_samples), drop = FALSE]

# Keep only common sample columns
common_samples <- intersect(colnames(refl_samples), colnames(trans_samples))
refl_samples <- refl_samples[, common_samples, drop = FALSE]
trans_samples <- trans_samples[, common_samples, drop = FALSE]

# ==========================================================
# Apply Savitzky-Golay smoothing (Safe Column-wise Processing)
# ==========================================================

# Using lapply instead of apply preserves data frame integrity and names
refl_smooth <- as.data.frame(lapply(refl_samples, function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sgolayfilt(x, p = 2, n = 11)
}))

trans_smooth <- as.data.frame(lapply(trans_samples, function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sgolayfilt(x, p = 2, n = 11)
}))

# Calculate Absorptance matrix upfront (Vectorized)
abs_smooth <- 1 - (refl_smooth + trans_smooth)

# ==========================================================
# Extract spectral metrics (Robust Array Processing)
# ==========================================================

samples <- colnames(refl_smooth)
n_samples <- length(samples)

# Pre-allocate output arrays with NA instead of 0 to catch omissions cleanly
max_refl <- wl_max_refl <- T_at_max_refl <- rep(NA_real_, n_samples)
max_trans <- wl_max_trans <- R_at_max_trans <- rep(NA_real_, n_samples)
min_abs <- wl_min_abs <- R_min_abs <- T_min_abs <- rep(NA_real_, n_samples)

for (i in seq_along(samples)) {
  r_val <- refl_smooth[[i]]
  t_val <- trans_smooth[[i]]
  a_val <- abs_smooth[[i]]
  
  # Identify indices where both R and T possess valid, non-NA figures
  valid_idx <- which(!is.na(r_val) & !is.na(t_val))
  if (length(valid_idx) == 0) next
  
  # Reflectance Metrics
  sub_idx_r <- which.max(r_val[valid_idx])
  actual_idx_r <- valid_idx[sub_idx_r]
  
  max_refl[i]       <- r_val[actual_idx_r]
  wl_max_refl[i]    <- wavelengths[actual_idx_r]
  T_at_max_refl[i]  <- t_val[actual_idx_r]
  
  # Transmittance Metrics
  sub_idx_t <- which.max(t_val[valid_idx])
  actual_idx_t <- valid_idx[sub_idx_t]
  
  max_trans[i]      <- t_val[actual_idx_t]
  wl_max_trans[i]   <- wavelengths[actual_idx_t]
  R_at_max_trans[i] <- r_val[actual_idx_t]
  
  # Absorptance Metrics
  sub_idx_a <- which.min(a_val[valid_idx])
  actual_idx_a <- valid_idx[sub_idx_a]
  
  min_abs[i]        <- a_val[actual_idx_a]
  wl_min_abs[i]     <- wavelengths[actual_idx_a]
  R_min_abs[i]      <- r_val[actual_idx_a]
  T_min_abs[i]      <- t_val[actual_idx_a]
}

# ==========================================================
# Build and Combine Output Table
# ==========================================================

max_LMA <- data.frame(
  Sample = samples,
  
  R_Max_Reflectance      = max_refl,
  T_At_Max_Reflectance   = T_at_max_refl,
  WL_Max_Reflectance     = wl_max_refl,
  
  R_At_Max_Transmittance = R_at_max_trans,
  T_Max_Transmittance    = max_trans,
  WL_Max_Transmittance   = wl_max_trans,
  
  R_Min_Absorptance      = R_min_abs,
  T_Min_Absorptance      = T_min_abs,
  Min_Absorptance        = min_abs,
  WL_Min_Absorptance     = wl_min_abs
)

# Join and structure downstream variables cleanly
max_LMA <- max_LMA %>%
  left_join(pigment, by = c("Sample" = "sampleid")) %>%
  select(
    Sample, Species,
    R_Max_Reflectance, T_At_Max_Reflectance, WL_Max_Reflectance,
    R_At_Max_Transmittance, T_Max_Transmittance, WL_Max_Transmittance,
    R_Min_Absorptance, T_Min_Absorptance, Min_Absorptance, WL_Min_Absorptance,
    LMA = LMA.mg.cm2., 
    EWT = Cw.mg.cm2., 
    Cab = CHLab_tot_ug_cm2
  )

# ==========================================================
# Export & Preview
# ==========================================================

write.table(max_LMA, "Data/Baringo_data/Merged_Data/max_R_T_LMA.txt", row.names = FALSE)
head(max_LMA)