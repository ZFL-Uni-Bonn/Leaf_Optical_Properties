# -----------------------------------------------------------
# Title: Spectra Smoothing and Plotting
# Description: Reads a spectra file, applies Savitzky-Golay smoothing, and plots original vs. smoothed.
# Author: Victor Korir
# -----------------------------------------------------------

# ---- Load Libraries ----
library(signal)
library(tidyverse)
library(tidyr)

# ---- Load External Functions ----
#source("Spectra_read_plot.R")  # Assumes this file is in the same folder

# ---- Define Helper Function ----
read_spectra <- function(file) {
  data <- read.table(file, skip = 41, header = FALSE, sep = ',')
  data$Spectrum <- tools::file_path_sans_ext(basename(file))  # Use file name as ID
  return(data)
}

# ---- Load and Visualize Spectra ----
spectra_file <- '/home/victor/Emap_R/exported_sorted_spectra/S26/prosopis_shaded/SPC.001.txt'
spc_004 <- read_spectra(spectra_file)

ggplot(spc_004, aes(x = V1, y = V2)) +
  geom_line(color = 'blue') +
  labs(x = 'Wavelength', y = 'Reflectance') +
  theme_minimal(base_size = 14)

# ---- Subset and Smooth the Data (1000???1300 nm) ----
subset_range <- spc_004$V1 >= 1000 & spc_004$V1 <= 1300
spc_sub <- spc_004[subset_range, ]

smoothed <- sgolayfilt(spc_sub$V2, p = 4, n = 31)

# ---- Update Smoothed Values in Full Dataset ----
spec_smth <- spc_004
spec_smth$V2[subset_range] <- smoothed

# ---- Prepare Data for Plotting ----
plot_df <- data.frame(
  wvl = spc_004$V1,
  original = spc_004$V2,
  smooth = spec_smth$V2
)

plot_df_long <- pivot_longer(plot_df, cols = c("original", "smooth"),
                             names_to = "Source", values_to = "Value")

# ---- Plot Original vs Smoothed ----
ggplot(plot_df_long, aes(x = wvl, y = Value, color = Source)) +
  geom_line(size = 0.7) +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance",
    color = "Source"
  ) +
  scale_x_continuous(limits = c(900, 1400)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black")
  )
