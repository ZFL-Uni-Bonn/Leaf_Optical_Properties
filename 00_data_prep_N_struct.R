# ==========================================================
# Prepare data for N structural parameter optimization
# Includes:
# 1. Savitzky-Golay smoothing
# 2. Matching pigments using sample IDs
# ==========================================================

library(signal)
library(dplyr)

# ==========================================================
# Load data
# ==========================================================

reflect <- read.table(
  "Data/Baringo_data/Reflectance_merged.txt",
  header = TRUE
)

transmit <- read.table(
  "Data/Baringo_data/Transmittance_T_merged.txt",
  header = TRUE
)

pigment <- read.table(
  "Data/Baringo_data/Merged_output.csv",
  sep = ",",
  header = TRUE,
  fill = TRUE
)

# ==========================================================
# Restrict wavelength range
# ==========================================================

reflect_750_900 <- reflect[
  reflect$wavelength >= 750 &
    reflect$wavelength <= 900,
]

transmit_750_900 <- transmit[
  transmit$wavelength >= 750 &
    transmit$wavelength <= 900,
]

# ==========================================================
# Extract wavelengths
# ==========================================================

wavelengths <- reflect_750_900$wavelength

# Remove wavelength column
refl_samples <- reflect_750_900[, -1]
trans_samples <- transmit_750_900[, -1]

# ==========================================================
# Apply Savitzky-Golay smoothing
# ==========================================================

# Polynomial order = 3
# Window size = 11 (must be odd)

refl_smooth <- as.data.frame(
  apply(refl_samples, 2, function(x)
    sgolayfilt(x, p = 2, n = 11))
)

trans_smooth <- as.data.frame(
  apply(trans_samples, 2, function(x)
    sgolayfilt(x, p = 2, n = 11))
)

# Preserve column names
colnames(refl_smooth) <- colnames(refl_samples)
colnames(trans_smooth) <- colnames(trans_samples)

# ==========================================================
# Initialize output vectors
# ==========================================================

samples <- colnames(refl_smooth)

max_refl <- numeric(length(samples))
wl_max_refl <- numeric(length(samples))

max_trans <- numeric(length(samples))
wl_max_trans <- numeric(length(samples))

# ==========================================================
# Extract maxima
# ==========================================================

for(i in seq_along(samples)){
  
  # Reflectance
  reflect_values <- refl_smooth[[i]]
  
  max_refl[i] <- max(reflect_values, na.rm = TRUE)
  
  wl_max_refl[i] <- wavelengths[
    which.max(reflect_values)
  ]
  
  # Transmittance
  trans_values <- trans_smooth[[i]]
  
  max_trans[i] <- max(trans_values, na.rm = TRUE)
  
  wl_max_trans[i] <- wavelengths[
    which.max(trans_values)
  ]
}

# ==========================================================
# Create spectral metrics table
# ==========================================================

max_LMA <- data.frame(
  Sample = samples,
  Max_Reflectance = max_refl,
  WL_Max_Reflectance = wl_max_refl,
  Max_Transmittance = max_trans,
  WL_Max_Transmittance = wl_max_trans
)

# ==========================================================
# Match pigments by sample ID
# ==========================================================

# Replace "SampleID" with actual column name
# from pigment table if different

max_LMA <- max_LMA %>%
  left_join(
    pigment,
    by = c("Sample" = "sampleid")
  )

# ==========================================================
# Select desired pigment variables
# ==========================================================

max_LMA <- max_LMA %>%
  select(
    Sample,
    Max_Reflectance,
    WL_Max_Reflectance,
    Max_Transmittance,
    WL_Max_Transmittance,
    LMA = LMA.mg.cm2.,
    EWT = Cw.mg.cm2.,
    Cab = CHLab_tot_ug_cm2
  )

# ==========================================================
# Export
# ==========================================================

write.csv(
  max_LMA,
  "Data/Baringo_data/max_R_T_LMA.csv",
  row.names = FALSE
)
