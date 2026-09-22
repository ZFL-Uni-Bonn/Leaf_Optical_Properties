################################################################################
# Title: Inversion of Leaf Optical Properties (LOP) using PROSPECT
# Author: Victor Korir
# Description: Inverts reflectance and transmittance spectra using PROSPECT model
#              to estimate biochemical traits such as CHL, CAR, EWT, and LMA.
################################################################################

# ---- Load Required Packages ----
library(prospect)
library(ggplot2)
library(prospect)

# ---- Load Input Data ----
refl <- read.table('spectro/MATLAB_scripts_ZM/R_front_final.txt', header = TRUE)
continuum_rspectra <- read.table('continuum_removed_spectra.csv', header = T,sep = ',')
trans <- read.table('spectro/MATLAB_scripts_ZM/T_front_final.txt', header = TRUE)
pigments <- read.table('spectro/MATLAB_scripts_ZM/inputs.txt', header = TRUE)

# Prior estimation of N using R only
sub_data <- FitSpectralData(lambda = continuum_rspectra$V1, Refl = continuum_rspectra[,3:53])
LeafDB <- download_LeafDB(dbName = 'ANGERS')

# ---- Prepare Spectral Data for PROSPECT Inversion ----
zm_data <- FitSpectralData(
  SpecPROSPECT = SpecPROSPECT_FullRange, 
  Refl = refl[,-1], 
  Tran = trans[,-1], 
  lambda = refl$WL, 
  UserDomain = refl$WL
)

# ---- Parameters to Estimate ----
Parms2Estimate <- c('CHL')
InitValues <- data.frame(CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0, 
                         EWT = 0.01, LMA = 0.01, N = 1.5)

# ---- Invert PROSPECT Using Optimal Wavelengths and continuum removed reflectance ----
ParmEst_R_OPT <- Invert_PROSPECT_OPT(SpecPROSPECT = sub_data$SpecPROSPECT, 
                                     lambda = sub_data$lambda, 
                                     Refl = sub_data$Refl, Tran = , 
                                     PROSPECT_version = 'D',
                                     Parms2Estimate = Parms2Estimate, 
                                     InitValues = InitValues)

# ---- Invert PROSPECT Using Optimal Wavelengths ----
res_WL <- Invert_PROSPECT(
  SpecPROSPECT  = zm_data$SpecPROSPECT, 
  Refl = zm_data$Refl,
  Tran = zm_data$Tran,
  PROSPECT_version = 'D',
  Parms2Estimate = params2estimate
)

# ---- Combine Estimated and Measured Values ----
estim_meas_df <- data.frame(
  Index = seq_along(res_opt_WL$CHL),
  Estimated = res_WL$CHL,
  Measured = pigments$Cab
)

# ---- Compute Accuracy Metrics ----
fit <- lm(Measured ~ Estimated, data = estim_meas_df)
r_squared <- summary(fit)$r.squared
rmse <- sqrt(mean((estim_meas_df$Measured - estim_meas_df$Estimated)^2))
metrics_text <- paste0("R² = ", round(r_squared, 3), "\nRMSE = ", round(rmse, 3))

# ---- Plot Estimated vs Measured CHL ----
ggplot(data = estim_meas_df, aes(x = Measured, y = Estimated)) +
  geom_point(color = 'steelblue', size = 2) +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
  annotate("text", 
           x = min(estim_meas_df$Measured), 
           y = max(estim_meas_df$Estimated), 
           label = metrics_text, 
           hjust = 0, vjust = 1,
           size = 5, color = "black") +
  labs(
    title = "PROSPECT Inversion: Measured vs Estimated CHL",
    x = 'Measured CHL',
    y = 'Estimated CHL'
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
