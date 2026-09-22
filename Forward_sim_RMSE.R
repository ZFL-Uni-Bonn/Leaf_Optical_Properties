# ==========================================================
# Forward PROSPECT simulation using retrieved N
#
# Computes:
# 1. Simulated reflectance
# 2. Simulated transmittance
# 3. RMSE between measured and simulated spectra
#
# Author: Victor Korir
# ==========================================================

# ==========================================================
# Load required packages
# ==========================================================

library(prospect)
library(dplyr)

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
# Load retrieved parameters
# ==========================================================

results <- read.csv(
  "Data/Baringo_data/estimated_N_LMA_results.csv"
)

# ==========================================================
# Restrict wavelength range
# ==========================================================

reflect <- reflect[
  reflect$wavelength >= 800 &
    reflect$wavelength <= 900,
]

transmit <- transmit[
  transmit$wavelength >= 800 &
    transmit$wavelength <= 900,
]

# ==========================================================
# Extract wavelength vector
# ==========================================================

wavelengths <- reflect$wavelength

# ==========================================================
# Remove wavelength column
# ==========================================================

refl_samples <- reflect[, -1]

trans_samples <- transmit[, -1]

# ==========================================================
# Fit spectral domain
# ==========================================================

adjust_vnir <- fit_spectral_data(
  lambda = wavelengths
)

# ==========================================================
# Storage vectors
# ==========================================================

rmse_reflectance <- numeric(nrow(results))

rmse_transmittance <- numeric(nrow(results))

rmse_Absorptance <- numeric(nrow(results))

mean_rmse <- numeric(nrow(results))

# ==========================================================
# Loop through samples
# ==========================================================

for(i in 1:nrow(results)){
  
  sample_name <- results$Sample[i]
  
  cat("Processing:", sample_name, "\n")
  
  # --------------------------------------------------------
  # Skip missing samples
  # --------------------------------------------------------
  
  if(!(sample_name %in% colnames(refl_samples))){
    
    rmse_reflectance[i] <- NA
    rmse_transmittance[i] <- NA
    mean_rmse[i] <- NA
    
    next
  }
  
  # --------------------------------------------------------
  # Measured spectra
  # --------------------------------------------------------
  
  R_measured <- refl_samples[[sample_name]]
  
  T_measured <- trans_samples[[sample_name]]
  
  # --------------------------------------------------------
  # Retrieved parameters
  # --------------------------------------------------------
  
  N_value <- results$Estimated_N[i]
  
  LMA <- results$LMA[i]
  
  EWT <- results$EWT[i]
  
  CHL <- results$Cab[i]
  
  # --------------------------------------------------------
  # Forward PROSPECT simulation
  # --------------------------------------------------------
  
  sim <- tryCatch(
    
    prospect(
      spec_prospect = adjust_vnir$spec_prospect,
      
      n_struct = N_value,
      
      lma = LMA/1000,
      
      ewt = EWT/1000
    ),
    
    error = function(e) return(NULL)
  )
  
  # --------------------------------------------------------
  # Handle failed simulations
  # --------------------------------------------------------
  
  if(is.null(sim)){
    
    rmse_reflectance[i] <- NA
    rmse_transmittance[i] <- NA
    mean_rmse[i] <- NA
    
    next
  }
  
  # --------------------------------------------------------
  # Simulated spectra
  # --------------------------------------------------------
  
  R_sim <- sim$reflectance
  
  T_sim <- sim$transmittance
  
  A_sim <- 1-(R_sim+T_sim)
  
  A_measured  <- 1-(R_measured+T_measured)
  
  # --------------------------------------------------------
  # RMSE reflectance
  # --------------------------------------------------------
  
  rmse_R <- sqrt(
    
    mean(
      (R_measured - R_sim)^2,
      na.rm = TRUE
    )
  )
  
  # --------------------------------------------------------
  # RMSE transmittance
  # --------------------------------------------------------
  
  rmse_T <- sqrt(
    
    mean(
      (T_measured - T_sim)^2,
      na.rm = TRUE
    )
  )
  
  # --------------------------------------------------------
  # RMSE Absorptance
  # --------------------------------------------------------
  
  rmse_A <- sqrt(
    
    mean(
      (A_measured - A_sim)^2,
      na.rm = TRUE
    )
  )
  
  
  # --------------------------------------------------------
  # Mean RMSE
  # --------------------------------------------------------
  
  mean_rmse_sample <- mean(
    c(rmse_R, rmse_T),
    na.rm = TRUE
  )
  
  # --------------------------------------------------------
  # Store
  # --------------------------------------------------------
  
  rmse_reflectance[i] <- rmse_R
  
  rmse_transmittance[i] <- rmse_T
  
  rmse_Absorptance[i] <- rmse_A
  
  mean_rmse[i] <- mean_rmse_sample
}

# ==========================================================
# Add RMSE metrics to results
# ==========================================================

results$RMSE_Reflectance <- rmse_reflectance

results$RMSE_Transmittance <- rmse_transmittance

results$RMSE_Absorptance <- rmse_Absorptance

results$Mean_RMSE <- mean_rmse

# ==========================================================
# Export
# ==========================================================

write.csv(
  
  results,
  
  "Data/Baringo_data/Merged_Data/forward_simulation_rmse.csv",
  
  row.names = FALSE
)

# ==========================================================
# Summary
# ==========================================================

summary(results$Mean_RMSE)

# ==========================================================
# Histogram
# ==========================================================

#Filtering results based on 1)N-Number 2) Absorptance RMSE 3) Reflectance RMSE
# 4) Transmittance RMSE

#Filterning based on N number
hist(results$Estimated_N) 

results_n <- results[
  !is.na(results$Estimated_N) &
    results$Estimated_N <2.8,
]
hist(results_n$Estimated_N)
#Filtering based on Absorptance RMSE
hist(results_n$RMSE_Absorptance)  

results_n <- results_n[
  !is.na(results_n$RMSE_Absorptance) &
    results_n$RMSE_Absorptance <0.066,
]
hist(results_n$RMSE_Absorptance)
#Filtering based on Reflectance RMSE
hist(results_n$RMSE_Reflectance)  

results_n <- results_n[
  !is.na(results_n$RMSE_Reflectance) &
    results_n$RMSE_Reflectance <0.06,
]
hist(results_n$RMSE_Reflectance)
#Filtering based on Transmittance RMSE
hist(results_n$RMSE_Transmittance)  

results_n <- results_n[
  !is.na(results_n$RMSE_Transmittance) &
    results_n$RMSE_Transmittance <0.06,
]
hist(results_n$RMSE_Transmittance)



write.csv(
  
  results_n,
  
  "Data/Baringo_data/Merged_Data/results_filtered_N_LMA_TW.csv",
  
  row.names = FALSE
)
# Plotting

library(tidyverse)

bins = ceiling(log2(nrow(results)) + 1)
x <- results$Estimated_N

breaks <- seq(min(x, na.rm = TRUE),
              max(x, na.rm = TRUE),
              length.out = bins + 1)
plot<-ggplot(
  results_n,
  aes(x = Estimated_N)
) +
  
  geom_histogram(
    # bins = round(sqrt(nrow(results_n))),
    bins = ceiling(log2(nrow(results)) + 1),
    breaks = breaks,
    fill = "skyblue",
    color = "black",
    linewidth = 0.01
  ) +
  
  scale_x_continuous(
    breaks = breaks, labels = round(breaks, 3))+
  
  geom_vline(
  #   aes(xintercept = mean(RMSE_Absorptance, na.rm = TRUE)),
    aes(xintercept = 2.8),
    linetype = "dashed",
    linewidth = 1
  ) +
  
  labs(
    title = "Estimated N",
    x = "N Value",
    y = "Frequency"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(
      face = "bold"
    )
  )

plot
filtered <- results[!results$Sample %in% results_n$Sample,]
#------------------------------------------------------------------------------------
library(tidyverse)

#----------------------------------------------------------
# Clean species names
#----------------------------------------------------------
results_plot <- results_n%>%
  mutate(
    Species = trimws(Species),
    Species = gsub("\\s+", " ", Species)
  )

#----------------------------------------------------------
# Keep species with enough observations
#----------------------------------------------------------
species_counts <- table(results_plot$Species)

keep_species <- names(species_counts[species_counts >= 5])

results_plot <- results_plot %>%
  filter(Species %in% keep_species)

#----------------------------------------------------------
# Order species by overall Mean_RMSE
#----------------------------------------------------------
species_order <- results_plot %>%
  group_by(Species) %>%
  summarise(
    Mean_RMSE = mean(Mean_RMSE, na.rm = TRUE)
  ) %>%
  arrange(desc(Mean_RMSE)) %>%
  pull(Species)

results_plot$Species <- factor(
  results_plot$Species,
  levels = species_order
)

#----------------------------------------------------------
# Function to generate violin + boxplot
#----------------------------------------------------------
plot_rmse <- function(data, variable, title){
  
  ggplot(
    data,
    aes(
      x = Species,
      y = .data[[variable]],
      fill = Species
    )
  ) +
    
    geom_violin(
      trim = FALSE,
      alpha = 0.6,
      colour = NA
    ) +
    
    geom_boxplot(
      width = 0.12,
      alpha = 0.9,
      outlier.shape = 16,
      outlier.size = 1
    ) +
    
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3,
      fill = "white"
    ) +
    
    labs(
      title = title,
      x = "Species",
      y = "RMSE"
    ) +
    
    theme_minimal(base_size = 14) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

#----------------------------------------------------------
# Generate plots
#----------------------------------------------------------
p_reflectance <- plot_rmse(
  results_plot,
  "RMSE_Reflectance",
  "Reflectance RMSE by Species"
)

p_transmittance <- plot_rmse(
  results_plot,
  "RMSE_Transmittance",
  "Transmittance RMSE by Species"
)

p_absorptance <- plot_rmse(
  results_plot,
  "RMSE_Absorptance",
  "Absorptance RMSE by Species"
)

# Display individually
p_reflectance
p_transmittance
p_absorptance
