################################################################################
# Title: Modified Merit function in the inversion to accomodate continuum removal
# Author: Victor Korir
# Description: Implements the continuum removal in the inversion merit function
################################################################################
# Clean workspace

# --- Load Required Libraries ---
library(prospect)
library(ggplot2)
library(quantreg)

# --- Custom Merit Function with Continuum Removal (650–750 nm) ---
Merit_PROSPECT_CRsim_670_722 <- function(x, SpecPROSPECT, Refl, Tran, Input_PROSPECT, Parms2Estimate) {
  # Ensure no negative parameter values
  x[x < 0] <- 0
  
  # Update parameters to estimate
  Input_PROSPECT[Parms2Estimate] <- x
  
  # Forward model simulation
  RT <- do.call("PROSPECT", c(list(SpecPROSPECT = SpecPROSPECT), Input_PROSPECT))
  
  fcr <- fct <- 0
  wavelengths <- SpecPROSPECT$lambda
  
  # Reflectance error using continuum removal
  if (!is.null(Refl)) {
    refl_sim_sub <- RT$Reflectance
    refl_meas_sub <- Refl
    
    # Define continuum between first and last point
    continuum <- approx(
      x = c(wavelengths[1], wavelengths[length(wavelengths)]),
      y = c(refl_sim_sub[1], refl_sim_sub[length(refl_sim_sub)]),
      xout = wavelengths,
      method = "linear"
    )$y
    
    # Continuum removal
    cr_sim <- 1 - (refl_sim_sub / continuum)
    
    # RMSE between measured and continuum-removed simulated reflectance
    fcr <- sqrt(mean((refl_meas_sub - cr_sim)^2))
  }
  
  # Transmittance error
  if (!is.null(Tran)) {
    fct <- sqrt(mean((Tran - RT$Transmittance)^2))
  }
  
  return(fcr)  # Only reflectance RMSE used here
}

# Assign custom merit function
assignInNamespace("Merit_PROSPECT_RMSE", Merit_PROSPECT_CRsim_670_722, ns = "prospect")

# --- Load Input Data ---
#refl <- read.table('spectro/MATLAB_scripts_ZM/R_front_final.txt', header = TRUE)
#refl <- read.table('/home/victor/LOP_Ecochange/trans_matched.txt', header = TRUE)
#refl <- read.table('/home/victor/Emap_R/Cw_Cm_cleaned/merged_mean_spectra.txt', header = TRUE)
refl <- read.table('/home/victor/Emap_R/select_spectra.txt', header = TRUE, fill = TRUE)
#trans <- read.table('spectro/MATLAB_scripts_ZM/T_front_final.txt', header = TRUE)
#pigments <- read.table('spectro/MATLAB_scripts_ZM/inputs.txt', header = TRUE)
#pigments <- read.table('/home/victor/LOP_Ecochange/pigments_subset_TR.txt', header = TRUE)
#refl <- read.table('/home/victor/spectra_front_analysis_ready.txt',  header = TRUE)
#pigments <- read.table('/home/victor/pigments_export.txt', header = TRUE)
# refl <- read.table('/home/victor/Emap_R/Cw_Cm_cleaned/combined_spectra.txt', sep = '\t', header = TRUE)
#refl <-refl[, !(names(refl) %in% c('X45'))]
pigments <- read.table('/home/victor/Emap_R/ewt.csv', sep = ',',header = T)


# pigments <- pigments[!(pigments$Sample.ID %in% c(10)), ]
pigments$EWT_median <- 0.02
pigments <- na.omit(pigments)
# pigments$EWT <- pigments$EWT

colnames(refl)[1] <- 'WL'
# refl['WL'] <- seq(450, 1650, by =1)
# refl <- refl[, c(ncol(refl), 1:(ncol(refl)-1))]


continuum_rspectra <- read.table('/home/victor/Emap_R//continuum_removed_spectra.csv', header = TRUE, sep = ',')

# Replace reflectance values with continuum-removed spectra in 659–737 nm
refl[refl$WL >=1100& refl$WL <= 1300,] <- continuum_rspectra[, 2:ncol(continuum_rspectra)]

# Subset reflectance to match spectral window
refl_sub <- refl[refl$WL >= 1100 & refl$WL <= 1300, ]

# Prepare spectral data
SubData <- FitSpectralData(
  lambda = refl_sub$WL,
  Refl = refl_sub[, 2:ncol(refl_sub)],
  Tran = NULL
)

# --- Inversion (Reflectance Only) ---
InitValues <- data.frame(
  CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0,
  EWT = 0.01, LMA = 0.01, N = 1.5
)

res_Ronly <- Invert_PROSPECT(
  SpecPROSPECT = SubData$SpecPROSPECT,
  Refl = SubData$Refl,
  Tran = SubData$Refl,
  PROSPECT_version = 'D',
  Parms2Estimate = c('CHL', 'EWT'),
  InitValues = InitValues,
  MeritFunction = Merit_PROSPECT_CRsim_670_722
)

# --- Combine Estimated and Measured Values ---
estim_meas_df <- data.frame(
  Index = seq_along(res_Ronly$EWT),
  Estimated = res_Ronly$EWT,
  Measured = pigments$EWT
)
#estim_meas_df <- estim_meas_df[which(estim_meas_df$Measured <= 0.1),]
#estim_meas_df <- estim_meas_df[45:60,]
#estim_meas_df <- estim_meas_df[estim_meas_df$Estimated < 0.05, ]
#estim_meas_df <- estim_meas_df[estim_meas_df$Measured < 0.06, ]
#estim_meas_df$Estimated <- estim_meas_df$Estimated- 0.011
#estim_meas_df <- estim_meas_df[estim_meas_df$Estimated<0.016&estim_meas_df$Measured<0.016,]
#estim_meas_df$Estimated <- estim_meas_df$Estimated*0.8
#estim_meas_df <- estim_meas_df[estim_meas_df$Measured < 0.07,]

# --- Accuracy Metrics ---
fit <- lm(Measured ~ Estimated, data = estim_meas_df)
r_squared <- summary(fit)$r.squared
rmse <- sqrt(mean((estim_meas_df$Measured - estim_meas_df$Estimated)^2))
pred_y <- predict(fit)

# fit_qr <- rq(Measured ~ Estimated, data = estim_meas_df, tau = 0.5)
# pred_y <- predict(fit_qr)
# summary_fit <- summary(fit_qr, se = "boot")  # bootstrap standard errors
# pseudo_r2 <- summary_fit$r.squared
# rmse <- sqrt(mean((estim_meas_df$Measured - pred_y)^2))
# cat("Pseudo R²:", pseudo_r2, "\n")
# cat("RMSE:", rmse, "\n")
# Willmott RMSE components
rmse_sys <- sqrt(mean((pred_y - estim_meas_df$Measured)^2))   
rmse_unsys <- sqrt(mean((estim_meas_df$Estimated - pred_y)^2))
metrics_text <- paste0(
  "R² = ", round(r_squared, 3), "\n",
  "RMSE = ", round(rmse, 4), "\n",
  "RMSE_s = ", round(rmse_sys, 4), "\n",
  "RMSE_u = ", round(rmse_unsys, 4))

# --- Plot Results ---
regression_p <- ggplot(data = estim_meas_df, aes(x = Measured, y = Estimated)) +
  geom_point(color = 'blue', size = 2) +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
  annotate("text",
           x = min(estim_meas_df$Measured),
           y = max(estim_meas_df$Estimated),
           label = metrics_text,
           hjust = 0, vjust = 1,
           size = 5, color = "black") +
  labs(
    title = "PROSPECT Inversion: Measured vs Estimated CHL",
    x = "Measured Cw",
    y = "Estimated Cw"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

pigments$CHL <- res_Ronly$CHL
pigments$Cw_sim <- res_Ronly$EWT

#write.csv( pigments, 'pigments_select.csv')
chl_hist <- ggplot(pigments, aes(x = Cw_sim)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Cw", x = "Estimated Cw", y = "Frequency") +
  theme_minimal()

site_bplot <- ggplot(pigments, aes(x = Site.ID, y = EW_mean, fill = Site.ID)) +
  geom_boxplot(alpha = 0.3) +
  labs(title = "Boxplot of Measured Cw by Site", x = "Site", y = "Cw") +
  theme_minimal()+
  scale_fill_brewer(palette = 'Dark2')
library(tidyverse)
medians <- pigments %>% 
  group_by(Site.ID) %>% 
  summarise(median_chl = median(CHL, na.rm = T))

ggplot(pigments, aes(x = Site.ID, y = CHL))+
  geom_boxplot(fill = 'lightgray', color = 'black')+
  geom_text(data = medians, aes(x = Site.ID, y = median_chl, label = round(median_chl, 1)),
            vjust = -0.5, size = 3.5, fontface = 'italic')+
  labs(title = 'Boxplot of CHL by site', x='Site', y = 'CHL')+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  
#ggsave('outputs/cw_measured_estim_woutliers.png', regression_p, dpi = 320, width = 10, height = 7,bg = 'white')
