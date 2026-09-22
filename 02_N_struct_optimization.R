#*TODO: Optimization of N struct parameter
#*
#*Requires: Leaf optical properties(transmittance and reflectance), LMA
#*
#*References:
#* https://www.is.uni-freiburg.de/resources/computational-economics/5_OptimizationR.pdf&ved=2ahUKEwiYn4nXxsWMAxXL0gIHHbe1OP4QFnoECAkQAQ&usg=AOvVaw2ASyuSQH_8RD1Ahbx5N-rA
#* https://github.com/jbferet/prospect
#*
# Author: Victor Korir

################################################################################
#Loading required packages
library('prospect')


#* Load data: Assuming you have a table of the leaf optical 
#* properties(Transmittance and Reflectance) at the wavelengths λ_R, λ_T, λ_A  and LMA
#* You can use the script data_prep_N_struct to format the reflectance and trasmittance data
prep_data <- read.csv('Data/Baringo_data/max_R_T_LMA.csv')

#Create empty vector to store estimated N
estimated_N_values <- numeric(nrow(prep_data))
estimated_N_values_optim <- numeric(nrow(prep_data))
# Optimization loop
for (i in 1:nrow(prep_data)) {
  # Extract wavelengths and measured values
  lambda_r <- prep_data$WL_Max_Reflectance[i]
  lambda_t <- prep_data$WL_Max_Transmittance[i]

  
  R_measured <- prep_data$Max_Reflectance[i]
  T_measured <- prep_data$Max_Transmittance[i]
  LMA <- prep_data$LMA[i]
  CHL <- prep_data$Cab[i]
  EWT <- prep_data$EWT[i]
  
  
  
  # RMSE function
  rmse_N <- function(params) {
    N_value <- params[1]
    
    sim <- PROSPECT(SpecPROSPECT = prospect::SpecPROSPECT_FullRange, N=N_value, CHL = CHL, EWT = EWT, LMA = LMA)
    
    #Subsetting reflectance and transmittance at the specific wavelenghts where maximums occur between 750-950nm
    R_sim <- sim$Reflectance[which(sim$wvl==round(lambda_r))]
    T_sim <- sim$Transmittance[which(sim$wvl==round(lambda_t))]
    
    #RMSE computation from measured and simulated data
    rmse <- sqrt(sum((R_measured - R_sim)^2 + (T_measured - T_sim)^2))
 
    return(rmse)
  }
  # Optimization using base optim
  #option1 - using optimize
  optim_result <- optimize(
    f = rmse_N,
    maximum = FALSE,
    interval = c(1, 3)  # Typical N range
  )
  
  # #option2 - using optim function
   optim_result1 <- optim(
     par = 1.5,
     fn = rmse_N,
    method = "L-BFGS-B",
     lower = 1,
     upper = 3
   )
  
  estimated_N_values[i] <- optim_result$minimum
  print(paste("Sample", i, "Estimated N =", optim_result$minimum, "RMSE =", round(optim_result$objective, 5)))
  #print(optim_result$convergence == 0 )#check if the optimization converged  #uncomment for option2
  estimated_N_values_optim[i] <- optim_result1$par  #uncomment for option2
}

#Combine raw dara and view results
results <- cbind(prep_data, Estimated_N = estimated_N_values)
write.table(estimated_N_values_optim, file = 'spectro/R_scripts/estimated_N_optim.txt', sep = '\t', row.names = F, quote = F)
