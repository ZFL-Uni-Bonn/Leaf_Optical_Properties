#* =========================================================
#* Optimization of PROSPECT N and LMA parameters
#*
#* Author: Victor Korir
#* =========================================================

library(prospect)
library(tidyverse)
library(DEoptim)

# ==========================================================
# Load prepared data
# ==========================================================

prep_data <- read.table(
  "Data/Baringo_data/Merged_Data/max_R_T_LMA.txt", header = TRUE
)

# ==========================================================
# Remove incomplete rows
# ==========================================================

prep_data <- prep_data %>%
  dplyr::filter(
    !is.na(R_Max_Reflectance),
    !is.na(T_At_Max_Reflectance),
    !is.na(R_At_Max_Transmittance),
    !is.na(T_Max_Transmittance),
    !is.na(R_Min_Absorptance),
    !is.na(T_Min_Absorptance)
  )

# ==========================================================
# Spectral range
# ==========================================================

wl_range <- seq(750, 900)

adjust_vnir <- fit_spectral_data(
  lambda = wl_range
)

# ==========================================================
# Storage vectors
# ==========================================================

estimated_N_values   <- numeric(nrow(prep_data))
estimated_LMA_values <- numeric(nrow(prep_data))
objective_values     <- numeric(nrow(prep_data))

# ==========================================================
# Optimization loop
# ==========================================================

for(i in 1:nrow(prep_data)){
  
  # --------------------------------------------------------
  # Selected wavelengths
  # --------------------------------------------------------
  
  lambda_r <- prep_data$WL_Max_Reflectance[i]
  lambda_t <- prep_data$WL_Max_Transmittance[i]
  lambda_a <- prep_data$WL_Min_Absorptance[i]
  
  # --------------------------------------------------------
  # Measured reflectance
  # --------------------------------------------------------
  
  R_measured <- c(
    prep_data$R_Max_Reflectance[i],
    prep_data$R_At_Max_Transmittance[i],
    prep_data$R_Min_Absorptance[i]
  )
  
  # --------------------------------------------------------
  # Measured transmittance
  # --------------------------------------------------------
  
  T_measured <- c(
    prep_data$T_At_Max_Reflectance[i],
    prep_data$T_Max_Transmittance[i],
    prep_data$T_Min_Absorptance[i]
  )
  
  # ========================================================
  # Merit function (two parameters)
  # ========================================================
  
  merit_function <- function(params){
    
    N_value <- params[1]
    lma     <- params[2]
    
    sim <- tryCatch(
      
      prospect(
        spec_prospect = adjust_vnir$spec_prospect,
        n_struct = N_value,
        lma = lma
      ),
      
      error = function(e) return(Inf)
    )
    
    if(is.null(sim)) return(Inf)
    
    # Find closest wavelengths
    idx_r <- which.min(abs(sim$wvl - lambda_r))
    idx_t <- which.min(abs(sim$wvl - lambda_t))
    idx_a <- which.min(abs(sim$wvl - lambda_a))
    
    # Simulated reflectance
    R_sim <- c(
      sim$reflectance[idx_r],
      sim$reflectance[idx_t],
      sim$reflectance[idx_a]
    )
    
    # Simulated transmittance
    T_sim <- c(
      sim$transmittance[idx_r],
      sim$transmittance[idx_t],
      sim$transmittance[idx_a]
    )
    
    # Merit function
    J <- sum((R_measured - R_sim)^2 + (T_measured - T_sim)^2)
    
    return(J)
  }
  
  # ========================================================
  # Optimization with DEoptim
  # ========================================================
  
  optim_result <- DEoptim(
    fn      = merit_function,
    lower   = c(1.0, 0.001),
    upper   = c(3.0, 0.060),
    control = DEoptim.control(
      itermax = 150,
      trace   = FALSE,
      reltol  = 1e-6
    )
  )
  
  # ========================================================
  # Store results
  # ========================================================
  
  estimated_N_values[i]   <- optim_result$optim$bestmem[1]
  estimated_LMA_values[i] <- optim_result$optim$bestmem[2]
  objective_values[i]     <- optim_result$optim$bestval
  
  # ========================================================
  # Progress
  # ========================================================
  
  print(
    paste(
      "Sample", i,
      "| N =", round(optim_result$optim$bestmem[1], 4),
      "| LMA =", round(optim_result$optim$bestmem[2], 6),
      "| J =", round(optim_result$optim$bestval, 6)
    )
  )
}

# ==========================================================
# Combine results
# ==========================================================

results <- cbind(
  prep_data,
  Estimated_N   = estimated_N_values,
  Estimated_LMA = estimated_LMA_values,
  Objective     = objective_values
)

# ==========================================================
# Export
# ==========================================================

write.csv(
  results,
  "Data/Baringo_data/estimated_N_LMA_results.csv",
  row.names = FALSE
)

# ==========================================================
# Visualization
# ==========================================================

# N distribution
ggplot(results, aes(x = Estimated_N)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", linewidth = 0.3) +
  geom_vline(aes(xintercept = mean(Estimated_N)), linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Estimated N", x = "N", y = "Frequency") +
  theme_minimal(base_size = 14)

# LMA distribution
ggplot(results, aes(x = Estimated_LMA)) +
  geom_histogram(bins = 30, fill = "coral", color = "black", linewidth = 0.3) +
  geom_vline(aes(xintercept = mean(Estimated_LMA)), linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Estimated LMA", x = "LMA (g/cm²)", y = "Frequency") +
  theme_minimal(base_size = 14)

# N vs LMA scatter
ggplot(results, aes(x = Estimated_N, y = Estimated_LMA)) +
  geom_point(alpha = 0.6) +
  labs(title = "Estimated N vs LMA", x = "N", y = "LMA (g/cm²)") +
  theme_minimal(base_size = 14)
