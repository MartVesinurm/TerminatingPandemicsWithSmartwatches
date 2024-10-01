library(MASS)
library(fitdistrplus)
library(dplyr)
library(ggthemes)
library(tidyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)
library(patchwork)
library(metR)

#THIS FINDS THE ECDF VL AT TIMEPOINT
extract_value <- function(time_index, vl_vector) {
  if(time_index > length(vl_vector)) {
    return(NA) # Returns NA if the time index is greater than the vector length
  } else {
    return(vl_vector[time_index])
  }
}

set.seed(123)
######ANCESTRAL######
ancestral_days <- 0:20 #Start from Infection not symptom onset
ancestral_VL <- c(1,1,1, 51, 82, 87, 87, 84, 73, 57, 43, 34, 26, 20, 14, 10, 5, 1, 1, 1, 1)
ancestral_df <- data.frame(day = ancestral_days, VL = ancestral_VL)
plot(ancestral_df, type = "l")
ancestral_vec <- rep(ancestral_df$day,ancestral_df$VL) #Gamma magic (multiply day by viral load)
ancestral_vec <- ifelse(ancestral_vec > 0, ancestral_vec, 0.001)
ancestral_gamma <- fitdist(ancestral_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(ancestral_gamma)
plot(ancestral_gamma)
ancestral_shape_param <- ancestral_gamma$estimate["shape"]
ancestral_rate_param <- ancestral_gamma$estimate["rate"]

ancestral_cohort <- data.frame(
  ID = 1:1000,
  Withdrawal = runif(1000, min = 0.50, max =0.75),
  Symptomatic = rbinom(1000, 1, prob = 0.844),
  Day_of_symptom_onset = 60,
  Detected_by_smartwatch = rbinom(1000, 1, prob = 0.88),
  Day_of_smartwatch_detection = 30, #Adding a note mostly for myself, for clarity, this variable means "How many days before symptom onset decetion occurs divided by 10"
  Time_from_onset_to_testing = rnorm(1000, mean = 21.7, sd = 0.9),
  R = rgamma(1000, shape = 13.608908, rate = 5.339734)
)

######DELTA######
delta_days <- 0:20 #Start from Infection not symptom onset
delta_VL <- c(1,1, 37, 97, 115, 117, 114, 105, 87, 68, 52, 42, 32, 24, 16, 10, 4, 1, 1, 1, 1)
delta_df <- data.frame(day = delta_days, VL = delta_VL)
plot(delta_df, type = "l")
delta_vec <- rep(delta_df$day,delta_df$VL) #Gamma magic (multiply day by viral load)
delta_vec <- ifelse(delta_vec > 0, delta_vec, 0.001)
delta_gamma <- fitdist(delta_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(delta_gamma)
plot(delta_gamma)
delta_shape_param <- delta_gamma$estimate["shape"]
delta_rate_param <- delta_gamma$estimate["rate"]

delta_cohort <- data.frame(
  ID = 1:1000,
  Withdrawal = runif(1000, min = 0.50, max =0.75),
  Symptomatic = rbinom(1000, 1, prob = 0.916),
  Day_of_symptom_onset = 50,
  Detected_by_smartwatch = rbinom(1000, 1, prob = 0.88),
  Day_of_smartwatch_detection = 30,
  Time_from_onset_to_testing = rnorm(1000, mean = 21.7, sd = 0.9),
  R = rnorm(1000, mean = 1.55, sd= 0.194))


######OMICRON######
omicron_days <- 0:20 #Start from Infection not symptom onset
omicron_VL <- c(1,11, 77, 95, 100, 100, 95, 84, 72, 53, 40, 33, 27, 21, 16, 9, 5, 1, 1, 1, 1)
omicron_df <- data.frame(day = omicron_days, VL = omicron_VL)
plot(omicron_df, type = "l")
omicron_vec <- rep(omicron_df$day,omicron_df$VL) #Gamma magic (multiply day by viral load)
omicron_vec <- ifelse(omicron_vec > 0, omicron_vec, 0.001)
omicron_gamma <- fitdist(omicron_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(omicron_gamma)
plot(omicron_gamma)
omicron_shape_param <- omicron_gamma$estimate["shape"]
omicron_rate_param <- omicron_gamma$estimate["rate"]

omicron_cohort <- data.frame(
  ID = 1:1000,
  Withdrawal = runif(1000, min = 0.50, max =0.75),
  Symptomatic = rbinom(1000, 1, prob = 0.745),
  Day_of_symptom_onset = 50,
  Detected_by_smartwatch = rbinom(1000, 1, prob = 0.88),
  Day_of_smartwatch_detection = 30,
  Time_from_onset_to_testing = rnorm(1000, mean = 21.7, sd = 0.9),
  R = rnorm(1000, mean = 4.2, sd = 1.097) #Du et al. 2022
)

######PANDEMIC INFLUENZA######

#SYMPTOMATICS#
pandemic_asym_days <- 0:20 #Start from Infection not symptom onset
pandemic_asym_VL <- c(1,10, 69, 100, 76, 83, 52, 51, 29, 9, 1, 1, 1, 1,1,1,1,1,1,1,1) #(Dennis et al., 2016)
pandemic_asym_df <- data.frame(day = pandemic_asym_days, VL = pandemic_asym_VL)
plot(pandemic_asym_df, type = "l")
pandemic_asym_vec <- rep(pandemic_asym_df$day,pandemic_asym_df$VL) #Gamma magic (multiply day by viral load)
pandemic_asym_vec <- ifelse(pandemic_asym_vec > 0, pandemic_asym_vec, 0.001)
pandemic_asym_gamma <- fitdist(pandemic_asym_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(pandemic_asym_gamma)
plot(pandemic_asym_gamma)
pandemic_shape_param <- pandemic_asym_gamma$estimate["shape"]
pandemic_rate_param <- pandemic_asym_gamma$estimate["rate"]


#ASYMPTOMATICS#
pandemic_asym_days <- 0:20 #Start from Infection not symptom onset
pandemic_asym_VL <- c(1,12,100, 61, 40, 45, 16, 34, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1) #(Dennis et al., 2016)
pandemic_asym_df <- data.frame(day = pandemic_asym_days, VL = pandemic_asym_VL)
plot(pandemic_asym_df, type = "l")
pandemic_asym_vec <- rep(pandemic_asym_df$day,pandemic_asym_df$VL) #Gamma magic (multiply day by viral load)
pandemic_asym_vec <- ifelse(pandemic_asym_vec > 0, pandemic_asym_vec, 0.001)
pandemic_asym_gamma <- fitdist(pandemic_asym_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(pandemic_asym_gamma)
plot(pandemic_asym_gamma)
pandemic_asym_shape_param <- pandemic_asym_gamma$estimate["shape"]
pandemic_asym_rate_param <- pandemic_asym_gamma$estimate["rate"]


pandemic_cohort <- data.frame(
  ID = 1:1000,
  Withdrawal = runif(1000, min = 0.50, max =0.75),
  Symptomatic = rbinom(1000, 1, prob = 0.773),
  Day_of_symptom_onset = 30,
  Detected_by_smartwatch = rbinom(1000, 1, prob = 0.9),
  Day_of_smartwatch_detection = 15,
  Time_from_onset_to_testing = rnorm(1000, mean = 21.7, sd = 0.9),
  R = rgamma(1000, shape = 25.86939 , rate = 16.61437)
)


######seasonal INFLUENZA######

#SYMPTOMATICS#
seasonal_asym_days <- 0:20 #Start from Infection not symptom onset
seasonal_asym_VL <- c(1,10, 69, 100, 76, 83, 52, 51, 29, 9, 1, 1, 1, 1,1,1,1,1,1,1,1) #(Dennis et al., 2016)
seasonal_asym_df <- data.frame(day = seasonal_asym_days, VL = seasonal_asym_VL)
plot(seasonal_asym_df, type = "l")
seasonal_asym_vec <- rep(seasonal_asym_df$day,seasonal_asym_df$VL) #Gamma magic (multiply day by viral load)
seasonal_asym_vec <- ifelse(seasonal_asym_vec > 0, seasonal_asym_vec, 0.001)
seasonal_asym_gamma <- fitdist(seasonal_asym_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(seasonal_asym_gamma)
plot(seasonal_asym_gamma)
seasonal_shape_param <- seasonal_asym_gamma$estimate["shape"]
seasonal_rate_param <- seasonal_asym_gamma$estimate["rate"]


#ASYMPTOMATICS#
seasonal_asym_days <- 0:20 #Start from Infection not symptom onset
seasonal_asym_VL <- c(1,12,100, 61, 40, 45, 16, 34, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1) #(Dennis et al., 2016)
seasonal_asym_df <- data.frame(day = seasonal_asym_days, VL = seasonal_asym_VL)
plot(seasonal_asym_df, type = "l")
seasonal_asym_vec <- rep(seasonal_asym_df$day,seasonal_asym_df$VL) #Gamma magic (multiply day by viral load)
seasonal_asym_vec <- ifelse(seasonal_asym_vec > 0, seasonal_asym_vec, 0.001)
seasonal_asym_gamma <- fitdist(seasonal_asym_vec, distr = "gamma", method = "mle") #Fit gamma with Maxmimum Likelihood
print(seasonal_asym_gamma)
plot(seasonal_asym_gamma)
seasonal_asym_shape_param <- seasonal_asym_gamma$estimate["shape"]
seasonal_asym_rate_param <- seasonal_asym_gamma$estimate["rate"]


seasonal_cohort <- data.frame(
  ID = 1:1000,
  Withdrawal = runif(1000, min = 0.50, max =0.75),
  Symptomatic = rbinom(1000, 1, prob = 0.79),
  Day_of_symptom_onset = 30,
  Detected_by_smartwatch = rbinom(1000, 1, prob = 0.94),
  Day_of_smartwatch_detection = 9.6,
  Time_from_onset_to_testing = rnorm(1000, mean = 16.3, sd = 0.96),
  R = rgamma(1000, shape = 91.99575, rate = 71.77707)
)



######SETTINGS FOR UNCERTAINTY PLOT1######
##WITHDRAWAL##
set.seed(123)
ancestral_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
delta_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
omicron_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
pandemic_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
seasonal_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS



#ANCESTRAL#
ancestral_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = ancestral_rate_param, shape = ancestral_shape_param)), simplify = TRUE)
ancestral_cohort$VL_ecdf <- lapply(ancestral_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
ancestral_cohort$day_of_smartwatch_withdrawal <- ancestral_cohort$Day_of_symptom_onset - ancestral_cohort$Day_of_smartwatch_detection
ancestral_cohort$time_at_natural_withdrawal <- round(ancestral_cohort$Day_of_symptom_onset + ancestral_cohort$Time_from_onset_to_testing,0)

#First calculate ECDF up to point of withdrawal
ancestral_cohort$VL_at_natural_withdrawal <- mapply(extract_value, ancestral_cohort$time_at_natural_withdrawal, ancestral_cohort$VL_ecdf)
ancestral_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, ancestral_cohort$day_of_smartwatch_withdrawal, ancestral_cohort$VL_ecdf)

#Then include the rest with the withdrawal coefficient, and without for those 
ancestral_cohort$VL_at_natural_withdrawal <- ifelse(ancestral_cohort$Symptomatic == 1, ancestral_cohort$VL_at_natural_withdrawal + (1-ancestral_cohort$VL_at_natural_withdrawal)*(1-ancestral_cohort$Withdrawal), 1)
ancestral_cohort$VL_at_smartwatch_withdrawal <- ifelse(ancestral_cohort$Detected_by_smartwatch == 1, ancestral_cohort$VL_at_smartwatch_withdrawal + (1-ancestral_cohort$VL_at_smartwatch_withdrawal)*(1-ancestral_cohort$Withdrawal), ancestral_cohort$VL_at_natural_withdrawal) 

#DELTA#
delta_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = delta_rate_param, shape = delta_shape_param)), simplify = TRUE)
delta_cohort$VL_ecdf <- lapply(delta_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
delta_cohort$day_of_smartwatch_withdrawal <- delta_cohort$Day_of_symptom_onset - delta_cohort$Day_of_smartwatch_detection
delta_cohort$time_at_natural_withdrawal <- round(delta_cohort$Day_of_symptom_onset + delta_cohort$Time_from_onset_to_testing,0)

#First calculate ECDF up to point of withdrawal
delta_cohort$VL_at_natural_withdrawal <- mapply(extract_value, delta_cohort$time_at_natural_withdrawal, delta_cohort$VL_ecdf)
delta_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, delta_cohort$day_of_smartwatch_withdrawal, delta_cohort$VL_ecdf)

#Then include the rest with the withdrawal coefficient, and without for those 
delta_cohort$VL_at_natural_withdrawal <- ifelse(delta_cohort$Symptomatic == 1, delta_cohort$VL_at_natural_withdrawal + (1-delta_cohort$VL_at_natural_withdrawal)*(1-delta_cohort$Withdrawal), 1)
delta_cohort$VL_at_smartwatch_withdrawal <- ifelse(delta_cohort$Detected_by_smartwatch == 1, delta_cohort$VL_at_smartwatch_withdrawal + (1-delta_cohort$VL_at_smartwatch_withdrawal)*(1-delta_cohort$Withdrawal), delta_cohort$VL_at_natural_withdrawal) 


#OMICRON#

omicron_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = omicron_rate_param, shape = omicron_shape_param)), simplify = TRUE)
omicron_cohort$VL_ecdf <- lapply(omicron_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
omicron_cohort$day_of_smartwatch_withdrawal <- omicron_cohort$Day_of_symptom_onset - omicron_cohort$Day_of_smartwatch_detection
omicron_cohort$time_at_natural_withdrawal <- round(omicron_cohort$Day_of_symptom_onset + omicron_cohort$Time_from_onset_to_testing,0)

#First calculate ECDF up to point of withdrawal
omicron_cohort$VL_at_natural_withdrawal <- mapply(extract_value, omicron_cohort$time_at_natural_withdrawal, omicron_cohort$VL_ecdf)
omicron_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, omicron_cohort$day_of_smartwatch_withdrawal, omicron_cohort$VL_ecdf)

#Then include the rest with the withdrawal coefficient, and without for those 
omicron_cohort$VL_at_natural_withdrawal <- ifelse(omicron_cohort$Symptomatic == 1, omicron_cohort$VL_at_natural_withdrawal + (1-omicron_cohort$VL_at_natural_withdrawal)*(1-omicron_cohort$Withdrawal), 1)
omicron_cohort$VL_at_smartwatch_withdrawal <- ifelse(omicron_cohort$Detected_by_smartwatch == 1, omicron_cohort$VL_at_smartwatch_withdrawal + (1-omicron_cohort$VL_at_smartwatch_withdrawal)*(1-omicron_cohort$Withdrawal), omicron_cohort$VL_at_natural_withdrawal) 

#PANDEMIC INFLUENZA#

pandemic_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = pandemic_rate_param, shape = pandemic_shape_param)), simplify = TRUE)
pandemic_cohort$VL[pandemic_cohort$Symptomatic == 0] <- replicate(n = length(pandemic_cohort$VL[pandemic_cohort$Symptomatic == 0]), list(rgamma(200, rate = pandemic_asym_rate_param, shape = pandemic_asym_shape_param)), simplify = TRUE)

pandemic_cohort$VL_ecdf <- lapply(pandemic_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
pandemic_cohort$VL_ecdf[pandemic_cohort$Symptomatic == 0] <- lapply(pandemic_cohort$VL_ecdf[pandemic_cohort$Symptomatic == 0], function(vl) vl * (2.7 / 4.7)) #Weights from Mean viral load from Table 3 of Dennis et al., 2016

pandemic_cohort$day_of_smartwatch_withdrawal <- pandemic_cohort$Day_of_symptom_onset - pandemic_cohort$Day_of_smartwatch_detection
pandemic_cohort$time_at_natural_withdrawal <- round(pandemic_cohort$Day_of_symptom_onset + pandemic_cohort$Time_from_onset_to_testing,0)

#First calculate ECDF up to point of withdrawal
pandemic_cohort$VL_at_natural_withdrawal <- mapply(extract_value, pandemic_cohort$time_at_natural_withdrawal, pandemic_cohort$VL_ecdf)
pandemic_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, pandemic_cohort$day_of_smartwatch_withdrawal, pandemic_cohort$VL_ecdf)


#Then include the rest with the withdrawal coefficient
pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1] <- pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1] + (1-pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1])*(1-pandemic_cohort$Withdrawal[pandemic_cohort$Symptomatic == 1])
pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 0] <- sapply(pandemic_cohort$VL_ecdf[pandemic_cohort$Symptomatic == 0], function(x) x[length(x)])

pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 0] <- pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 0]
pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 0] <- pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 0]
pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1] <- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1] + (1- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1]) * (1-pandemic_cohort$Withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1])
pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1] <- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1] + ((2.7 / 4.7)- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1]) * (1-pandemic_cohort$Withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1])

#SEASONAL INFLUENZA#
seasonal_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = seasonal_rate_param, shape = seasonal_shape_param)), simplify = TRUE)
seasonal_cohort$VL[seasonal_cohort$Symptomatic == 0] <- replicate(n = length(seasonal_cohort$VL[seasonal_cohort$Symptomatic == 0]), list(rgamma(200, rate = seasonal_asym_rate_param, shape = seasonal_asym_shape_param)), simplify = TRUE)

seasonal_cohort$VL_ecdf <- lapply(seasonal_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
seasonal_cohort$VL_ecdf[seasonal_cohort$Symptomatic == 0] <- lapply(seasonal_cohort$VL_ecdf[seasonal_cohort$Symptomatic == 0], function(vl) vl * (4.0 / 5.6)) #Weights from Mean viral load from Table 3 of Dennis et al., 2016

seasonal_cohort$day_of_smartwatch_withdrawal <- seasonal_cohort$Day_of_symptom_onset - seasonal_cohort$Day_of_smartwatch_detection
seasonal_cohort$time_at_natural_withdrawal <- round(seasonal_cohort$Day_of_symptom_onset + seasonal_cohort$Time_from_onset_to_testing,0)

#First calculate ECDF up to point of withdrawal
seasonal_cohort$VL_at_natural_withdrawal <- mapply(extract_value, seasonal_cohort$time_at_natural_withdrawal, seasonal_cohort$VL_ecdf)
seasonal_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, seasonal_cohort$day_of_smartwatch_withdrawal, seasonal_cohort$VL_ecdf)


#Then include the rest with the withdrawal coefficient
seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1] <- seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1] + (1-seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1])*(1-seasonal_cohort$Withdrawal[seasonal_cohort$Symptomatic == 1])
seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 0] <- sapply(seasonal_cohort$VL_ecdf[seasonal_cohort$Symptomatic == 0], function(x) x[length(x)])

seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 0] <- seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 0]
seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 0] <- seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 0]
seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1] <- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1] + (1- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1]) * (1-seasonal_cohort$Withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1])
seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1] <- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1] + ((4.0 / 5.6)- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1]) * (1-seasonal_cohort$Withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1])



####R's with SMARTWATCH

ancestral_cohort$R_sw <- (ancestral_cohort$VL_at_smartwatch_withdrawal * ancestral_cohort$R)/ ancestral_cohort$VL_at_natural_withdrawal
delta_cohort$R_sw <- (delta_cohort$VL_at_smartwatch_withdrawal * delta_cohort$R)/ delta_cohort$VL_at_natural_withdrawal
omicron_cohort$R_sw <- (omicron_cohort$VL_at_smartwatch_withdrawal * omicron_cohort$R)/ omicron_cohort$VL_at_natural_withdrawal
pandemic_cohort$R_sw <- (pandemic_cohort$VL_at_smartwatch_withdrawal * pandemic_cohort$R)/ pandemic_cohort$VL_at_natural_withdrawal
seasonal_cohort$R_sw <- (seasonal_cohort$VL_at_smartwatch_withdrawal * seasonal_cohort$R)/ seasonal_cohort$VL_at_natural_withdrawal


plot(density(ancestral_cohort$R_sw))
plot(density(delta_cohort$R_sw))
plot(density(omicron_cohort$R_sw))
plot(density(pandemic_cohort$R_sw))
plot(density(seasonal_cohort$R_sw))

#ANCESTRAL PROBS
ancestral_cohort$probs <- lapply(1:nrow(ancestral_cohort), function(x) rep(NA, 251))
for(i in 1:1000){
  if(ancestral_cohort$Detected_by_smartwatch[i] == 1) {
    ancestral_cohort$probs[[i]] <- ancestral_cohort$VL_ecdf[[i]] + (1- ancestral_cohort$VL_ecdf[[i]]) * (1-ancestral_cohort$Withdrawal[i]) #IF DE
  } else if (ancestral_cohort$Symptomatic[i] == 1) {
    ancestral_cohort$probs[[i]] <-   rep(ancestral_cohort$VL_at_natural_withdrawal[i], 251)   
  } else {
    ancestral_cohort$probs[[i]] <- rep(1, 251)
  }                                     
  ancestral_cohort$probs[[i]] <- ancestral_cohort$probs[[i]] * ancestral_cohort$R[i] / ancestral_cohort$VL_at_natural_withdrawal[i]
}
ancestral_cohort$probs_by_day <- lapply(ancestral_cohort$probs, function(x) x[seq(1, length(x), by = 10)])

#DELTA PROBS
delta_cohort$probs <- lapply(1:nrow(delta_cohort), function(x) rep(NA, 251))
for(i in 1:1000){
  if(delta_cohort$Detected_by_smartwatch[i] == 1) {
    delta_cohort$probs[[i]] <- delta_cohort$VL_ecdf[[i]] + (1- delta_cohort$VL_ecdf[[i]]) * (1-delta_cohort$Withdrawal[i]) #IF DE
  } else if (delta_cohort$Symptomatic[i] == 1) {
    delta_cohort$probs[[i]] <-   rep(delta_cohort$VL_at_natural_withdrawal[i], 251)   
  } else {
    delta_cohort$probs[[i]] <- rep(1, 251)
  }                                     
  delta_cohort$probs[[i]] <- delta_cohort$probs[[i]] * delta_cohort$R[i] / delta_cohort$VL_at_natural_withdrawal[i]
}
delta_cohort$probs_by_day <- lapply(delta_cohort$probs, function(x) x[seq(1, length(x), by = 10)])

#OMICRON PROBS
omicron_cohort$probs <- lapply(1:nrow(omicron_cohort), function(x) rep(NA, 251))
for(i in 1:1000){
  if(omicron_cohort$Detected_by_smartwatch[i] == 1) {
    omicron_cohort$probs[[i]] <- omicron_cohort$VL_ecdf[[i]] + (1- omicron_cohort$VL_ecdf[[i]]) * (1-omicron_cohort$Withdrawal[i]) #IF DE
  } else if (omicron_cohort$Symptomatic[i] == 1) {
    omicron_cohort$probs[[i]] <-   rep(omicron_cohort$VL_at_natural_withdrawal[i], 251)   
  } else {
    omicron_cohort$probs[[i]] <- rep(1, 251)
  }                                     
  omicron_cohort$probs[[i]] <- omicron_cohort$probs[[i]] * omicron_cohort$R[i] / omicron_cohort$VL_at_natural_withdrawal[i]
}
omicron_cohort$probs_by_day <- lapply(omicron_cohort$probs, function(x) x[seq(1, length(x), by = 10)])

#PANDEMIC PROBS
pandemic_cohort$probs <- lapply(1:nrow(pandemic_cohort), function(x) rep(NA, 251))
for(i in 1:1000){
  if(pandemic_cohort$Detected_by_smartwatch[i] == 1) {
    pandemic_cohort$probs[[i]] <- pandemic_cohort$VL_ecdf[[i]] + ( (2.7/4.7) - pandemic_cohort$VL_ecdf[[i]]) * (1-pandemic_cohort$Withdrawal[i]) #IF DE
  } else if (pandemic_cohort$Symptomatic[i] == 1) {
    pandemic_cohort$probs[[i]] <-   rep(pandemic_cohort$VL_at_natural_withdrawal[i], 251)   
  } else {
    pandemic_cohort$probs[[i]] <- rep(1, 251)
  }                                     
  pandemic_cohort$probs[[i]] <- pandemic_cohort$probs[[i]] * pandemic_cohort$R[i] / pandemic_cohort$VL_at_natural_withdrawal[i]
}
pandemic_cohort$probs_by_day <- lapply(pandemic_cohort$probs, function(x) x[seq(1, length(x), by = 10)])


#SEASONAL PROBS
seasonal_cohort$probs <- lapply(1:nrow(seasonal_cohort), function(x) rep(NA, 251))
for(i in 1:1000){
  if(seasonal_cohort$Detected_by_smartwatch[i] == 1) {
    seasonal_cohort$probs[[i]] <- seasonal_cohort$VL_ecdf[[i]] + ( (4.0 / 5.6) - seasonal_cohort$VL_ecdf[[i]]) * (1-seasonal_cohort$Withdrawal[i]) #IF DE
  } else if (seasonal_cohort$Symptomatic[i] == 1) {
    seasonal_cohort$probs[[i]] <-   rep(seasonal_cohort$VL_at_natural_withdrawal[i], 251)   
  } else {
    seasonal_cohort$probs[[i]] <- rep(1, 251)
  }                                     
  seasonal_cohort$probs[[i]] <- seasonal_cohort$probs[[i]] * seasonal_cohort$R[i] / seasonal_cohort$VL_at_natural_withdrawal[i]
}
seasonal_cohort$probs_by_day <- lapply(seasonal_cohort$probs, function(x) x[seq(1, length(x), by = 10)])






ancestral_uncertainty_analysis <- data.frame(t(as.data.frame(ancestral_cohort$probs_by_day)))
rownames(ancestral_uncertainty_analysis) <- 1:1000
colnames(ancestral_uncertainty_analysis) <- 0:25
ancestral_uncertainty_analysis <- ifelse(ancestral_uncertainty_analysis < 1, 1, 0)
ancestral_uncertainty_analysis <- colSums(ancestral_uncertainty_analysis) / nrow(ancestral_uncertainty_analysis)
ancestral_uncertainty_analysis <- data.frame(Ancestral = ancestral_uncertainty_analysis)


delta_uncertainty_analysis <- data.frame(t(as.data.frame(delta_cohort$probs_by_day)))
rownames(delta_uncertainty_analysis) <- 1:1000
colnames(delta_uncertainty_analysis) <- 0:25
delta_uncertainty_analysis <- ifelse(delta_uncertainty_analysis < 1, 1, 0)
delta_uncertainty_analysis <- colSums(delta_uncertainty_analysis) / nrow(delta_uncertainty_analysis)
delta_uncertainty_analysis <- data.frame(Delta = delta_uncertainty_analysis)


omicron_uncertainty_analysis <- data.frame(t(as.data.frame(omicron_cohort$probs_by_day)))
rownames(omicron_uncertainty_analysis) <- 1:1000
colnames(omicron_uncertainty_analysis) <- 0:25
omicron_uncertainty_analysis <- ifelse(omicron_uncertainty_analysis < 1, 1, 0)
omicron_uncertainty_analysis <- colSums(omicron_uncertainty_analysis) / nrow(omicron_uncertainty_analysis)
omicron_uncertainty_analysis <- data.frame(Omicron = omicron_uncertainty_analysis)


pandemic_uncertainty_analysis <- data.frame(t(as.data.frame(pandemic_cohort$probs_by_day)))
rownames(pandemic_uncertainty_analysis) <- 1:1000
colnames(pandemic_uncertainty_analysis) <- 0:25
pandemic_uncertainty_analysis <- ifelse(pandemic_uncertainty_analysis < 1, 1, 0)
pandemic_uncertainty_analysis <- colSums(pandemic_uncertainty_analysis) / nrow(pandemic_uncertainty_analysis)
pandemic_uncertainty_analysis <- data.frame(Pandemic = pandemic_uncertainty_analysis)

seasonal_uncertainty_analysis <- data.frame(t(as.data.frame(seasonal_cohort$probs_by_day)))
rownames(seasonal_uncertainty_analysis) <- 1:1000
colnames(seasonal_uncertainty_analysis) <- 0:25
seasonal_uncertainty_analysis <- ifelse(seasonal_uncertainty_analysis < 1, 1, 0)
seasonal_uncertainty_analysis <- colSums(seasonal_uncertainty_analysis) / nrow(seasonal_uncertainty_analysis)
seasonal_uncertainty_analysis <- data.frame(Seasonal = seasonal_uncertainty_analysis)




uncertainty_analysis <- data.frame(ancestral = ancestral_uncertainty_analysis,
                                   delta = delta_uncertainty_analysis,
                                   omicron = omicron_uncertainty_analysis,
                                   pandemic = pandemic_uncertainty_analysis,
                                   seasonal = seasonal_uncertainty_analysis)

uncertainty_covid <- data.frame(ancestral = ancestral_uncertainty_analysis,
                                delta = delta_uncertainty_analysis,
                                omicron = omicron_uncertainty_analysis)

uncertainty_influenza<- data.frame(pandemic = pandemic_uncertainty_analysis,
                                   seasonal = seasonal_uncertainty_analysis)


covid_onset <- data.frame(x = c(5,5,6),
                          y = c(0.01, 0.275, 0.065))


uncertainty_covid_ggplot <- ggplot() +
  geom_smooth(data = uncertainty_covid, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Ancestral, linetype="Ancestral"), color = "black")+
  geom_smooth(data = uncertainty_covid, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Delta, linetype="Delta"),color = "black") +
  geom_smooth(data = uncertainty_covid, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Omicron, linetype="Omicron"),color = "black") +
  geom_point(aes(x=2,y=0.048,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=2,y=0.875,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=3,y=0.253,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=5,y=0.01,colour="Symptom onset"), size=4)+
  geom_point(aes(x=5,y=0.268,colour="Symptom onset"), size=4)+
  geom_point(aes(x=6,y=0.06,colour="Symptom onset"), size=4)+
  scale_color_manual(values = c("Smartwatch detection" = "#88CCEE", "Symptom onset" = "#CC3311", "black" = "black"),
                     labels = c("Smartwatch detection" = "Smartwatch detection", "Symptom onset" = "Symptom onset")) +
  scale_x_continuous(expand= c(0,0), limits = c(0, 10), breaks= c(0,2,4,6,8,10)) +
  scale_y_continuous(expand= c(0,0), limits = c(0, 1), labels = scales::percent)  +  
  labs(title = "",
       x = "Time from Exposure (Days)",
       y = "Probability of the effective reproduction
       number being under 1 (%)") +
  guides(
    color = guide_legend(order = 2, title = "Point Types"),
    linetype = guide_legend(order = 1, title = "Line Types")
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        legend.title=element_blank(),
        legend.position = c(0.75, 0.8),
        text = element_text(size = 25),  # Adjust text size for readability
        plot.title = element_text(size = 25, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 25, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 25, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 25)  # Caption styling
  )

uncertainty_influenza_ggplot <- ggplot() +
  geom_smooth(data = uncertainty_influenza, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Pandemic, linetype="Pandemic"), color = "black")+
  geom_smooth(data = uncertainty_influenza, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Seasonal, linetype="Seasonal"),color = "black") +
  geom_point(aes(x=0.5,y=0.95,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=1,y=0.98,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=2,y=0.93,colour="Symptom onset"), size=4)+
  geom_point(aes(x=2,y=0.82,colour="Symptom onset"), size=4)+
  
  scale_color_manual(values = c("Smartwatch detection" = "#88CCEE", "Symptom onset" = "#CC3311", "black" = "black"),
                     labels = c("Smartwatch detection" = "Smartwatch detection", "Symptom onset" = "Symptom onset")) +
  
  
  scale_x_continuous(expand= c(0,0), limits = c(0, 10), breaks= c(0,2,4,6,8,10)) +
  scale_y_continuous(expand= c(0,0), limits = c(0, 1), labels = scales::percent)  +  
  labs(title = "",
       x = "Time from Exposure (Days)",
       y = "") +
  guides(
    color = guide_legend(order = 2, title = "Point Types"),
    linetype = guide_legend(order = 1, title = "Line Types")
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        legend.title=element_blank(),
        legend.position = c(0.75, 0.8),
        text = element_text(size = 25),  # Adjust text size for readability
        plot.title = element_text(size = 25, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 25, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 25, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 25)  # Caption styling
  )

Figure3 <- uncertainty_covid_ggplot + uncertainty_influenza_ggplot
Figure3 # Run this for figure 3
#ggsave("240929_Figure3.jpeg")
#Check for numerical results in "uncertainty_influenza" and "uncertainty_covid"