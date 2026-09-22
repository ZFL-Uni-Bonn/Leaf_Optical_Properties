#* =========================================================
#* Optimization of PROSPECT N structural parameter
#*
#* Merit function (over 3 selected wavelengths λR, λT, λA):
#*   J = Σ[(Rmes - Rmod)^2 + (Tmes - Tmod)^2]
#*
#* Each sample is keyed by Sample_ID.
#* Author: Victor Korir
#* =========================================================

library(prospect)
library(tidyverse)

# ---------------------------------------------------------
# Load and clean data
# ---------------------------------------------------------
prep_data <- read.table(
  "Data/Baringo_data/Merged_Data/max_R_T_LMA.txt",
  header = TRUE
) %>%
  drop_na(
    R_Max_Reflectance, T_At_Max_Reflectance,
    R_At_Max_Transmittance, T_Max_Transmittance,
    R_Min_Absorptance,  T_Min_Absorptance,
    LMA
  )

# Fail loudly if the ID column is missing rather than processing positionally
id_col <- "Sample"            # <-- set this to your actual identifier column
stopifnot(id_col %in% names(prep_data))

# ---------------------------------------------------------
# Spectral setup (computed ONCE, reused for every sample)
# ---------------------------------------------------------
wl_range   <- seq(750, 900)
adjust_vnir <- fit_spectral_data(lambda = wl_range)
spec        <- adjust_vnir$spec_prospect

# ---------------------------------------------------------
# Per-sample optimization, returning a tidy row
# ---------------------------------------------------------
optimize_one <- function(row) {
  
  lambda_r <- row$WL_Max_Reflectance
  lambda_t <- row$WL_Max_Transmittance
  lambda_a <- row$WL_Min_Absorptance
  
  R_measured <- c(row$R_Max_Reflectance,
                  row$R_At_Max_Transmittance,
                  row$R_Min_Absorptance)
  
  T_measured <- c(row$T_At_Max_Reflectance,
                  row$T_Max_Transmittance,
                  row$T_Min_Absorptance)
  
  LMA <- as.numeric(row$LMA)
  EWT <- as.numeric(row$EWT)
  CHL <- as.numeric(row$Cab)
  
  merit_function <- function(N_value) {
    
    sim <- tryCatch(
      prospect(
        spec_prospect = spec,
        n_struct      = N_value,
        lma           = LMA/1000,
      ),
      error = function(e) NULL
    )
    if (is.null(sim)) return(Inf)
    
    # nearest simulated wavelength to each target
    idx <- c(
      which.min(abs(sim$wvl - lambda_r)),
      which.min(abs(sim$wvl - lambda_t)),
      which.min(abs(sim$wvl - lambda_a))
    )
    
    R_sim <- sim$reflectance[idx]
    T_sim <- sim$transmittance[idx]
    
    sum((R_measured - R_sim)^2 + (T_measured - T_sim)^2)
  }
  
  opt <- optimize(merit_function, interval = c(1, 3), maximum = FALSE)
  
  tibble(
    Estimated_N        = opt$minimum,
    Objective_Function = opt$objective
  )
}

# ---------------------------------------------------------
# Run for every sample, keyed by ID, with progress
# ---------------------------------------------------------
results <- prep_data %>%
  mutate(.row = row_number()) %>%
  group_split(.row) %>%
  map_dfr(function(row) {
    out <- optimize_one(row)
    message(sprintf(
      "Sample %s | N = %.4f | J = %.6f",
      row[[id_col]], out$Estimated_N, out$Objective_Function
    ))
    bind_cols(row, out)
  }) %>%
  select(-.row)

# ---------------------------------------------------------
# Export
# ---------------------------------------------------------
write.csv(results, "Data/Baringo_data/estimated_N_results.csv",
          row.names = FALSE)

# ---------------------------------------------------------
# Distribution of estimated N
# ---------------------------------------------------------
ggplot(results, aes(x = Estimated_N)) +
  geom_histogram(bins = 30, fill = "skyblue",
                 color = "black", linewidth = 0.3) +
  geom_vline(aes(xintercept = mean(Estimated_N, na.rm = TRUE)),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of N values",
       x = "Estimated N", y = "Frequency") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title      = element_text(face = "bold")
  )
