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

######BEGIN######

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
  Day_of_smartwatch_detection = 30,
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


##########CALIBRATE MODEL IN THIS SECTION##########

######SETTINGS FOR BOXPLOT 1######
##WITHDRAWAL##
set.seed(123)
ancestral_cohort$Withdrawal <- 0.66 # RUN THIS FOR THE FIRST BOXPLOT FIGURE
delta_cohort$Withdrawal <- 0.66 # RUN THIS FOR THE FIRST BOXPLOT FIGURE
omicron_cohort$Withdrawal <- 0.66 # RUN THIS FOR THE FIRST BOXPLOT FIGURE
pandemic_cohort$Withdrawal <- 0.66 # RUN THIS FOR THE FIRST BOXPLOT FIGURE
seasonal_cohort$Withdrawal <- 0.66 # RUN THIS FOR THE FIRST BOXPLOT FIGURE

######SETTINGS FOR BOXPLOT 2######
##WITHDRAWAL##
set.seed(123)
ancestral_cohort$Withdrawal <- 0.75 # RUN THIS FOR THE SECOND BOXPLOT FIGURE
delta_cohort$Withdrawal <- 0.75 # RUN THIS FOR THE SECOND BOXPLOT FIGURE
omicron_cohort$Withdrawal <- 0.75 # RUN THIS FOR THE SECOND BOXPLOT FIGURE
pandemic_cohort$Withdrawal <- 0.75 # RUN THIS FOR THE SECOND BOXPLOT FIGURE
seasonal_cohort$Withdrawal <- 0.75 # RUN THIS FOR THE SECOND BOXPLOT FIGURE


######SETTINGS FOR UNCERTAINTY PLOT1######
##WITHDRAWAL##
set.seed(123)
ancestral_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
delta_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
omicron_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
pandemic_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS
seasonal_cohort$Withdrawal <- runif(1000, min = 0.66, max =0.75) # RUN THIS FOR UNCERTAINTY ANALYSIS



##################################################


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















#FIGURE 2

#Ancestral Boxplot
ancestral_boxplot_n <- data.frame(R = ancestral_cohort$R,
                                  Withdrawal = "Natural")
ancestral_boxplot_sw <- data.frame(R = ancestral_cohort$R_sw,
                                  Withdrawal = "Smartwatch")
ancestral_boxplot <- rbind(ancestral_boxplot_n,ancestral_boxplot_sw)
ancestral_boxplot$disease <- "Ancestral"

#Delta Boxplot
delta_boxplot_n <- data.frame(R = delta_cohort$R,
                                  Withdrawal = "Natural")
delta_boxplot_sw <- data.frame(R = delta_cohort$R_sw,
                                   Withdrawal = "Smartwatch")
delta_boxplot <- rbind(delta_boxplot_n,delta_boxplot_sw)
delta_boxplot$disease <- "Delta"

#Omicron Boxplot
omicron_boxplot_n <- data.frame(R = omicron_cohort$R,
                              Withdrawal = "Natural")
omicron_boxplot_sw <- data.frame(R = omicron_cohort$R_sw,
                               Withdrawal = "Smartwatch")
omicron_boxplot <- rbind(omicron_boxplot_n,omicron_boxplot_sw)
omicron_boxplot$disease <- "Omicron"

#Pandemic Boxplot
pandemic_boxplot_n <- data.frame(R = pandemic_cohort$R,
                                Withdrawal = "Natural")
pandemic_boxplot_sw <- data.frame(R = pandemic_cohort$R_sw,
                                 Withdrawal = "Smartwatch")
pandemic_boxplot <- rbind(pandemic_boxplot_n,pandemic_boxplot_sw)
pandemic_boxplot$disease <- "Pandemic"

#Seasonal Boxplot
seasonal_boxplot_n <- data.frame(R = seasonal_cohort$R,
                                 Withdrawal = "Natural")
seasonal_boxplot_sw <- data.frame(R = seasonal_cohort$R_sw,
                                  Withdrawal = "Smartwatch")
seasonal_boxplot <- rbind(seasonal_boxplot_n,seasonal_boxplot_sw)
seasonal_boxplot$disease <- "Seasonal"


#COMBINE BOXPLOTS#
boxplot_dataframe <- do.call(rbind, list(ancestral_boxplot,
                                         delta_boxplot,
                                         omicron_boxplot,
                                         pandemic_boxplot,
                                         seasonal_boxplot))

#COMBINE COVID BOXPLOTS#
boxplot_covid <- do.call(rbind, list(ancestral_boxplot,
                                         delta_boxplot,
                                         omicron_boxplot))

#COMBINE INFLUENZAS#
boxplot_influenza <- do.call(rbind, list(pandemic_boxplot,
                                         seasonal_boxplot))



###PLOT BOXPLOT###
covid_ggplot <- ggplot(boxplot_covid, aes(x = disease, y = R, fill = Withdrawal)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept=1, color="darkred", linetype="dashed", linewidth =1) +
  labs(title = "",
       x = "Covid Variants",
       y = "Reproduction Number") +
  scale_fill_manual(values = c("Natural" = "#CC3311", "Smartwatch" = "#88CCEE")) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        legend.title=element_blank(),
        legend.position = "none",
        text = element_text(size = 15),  # Adjust text size for readability
        plot.title = element_text(size = 15, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 15, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 15, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 15)  # Caption styling
  )  + scale_y_continuous(expand = c(0,0), breaks=seq(0,7, by=1), limits = c(0,7.2))

influenza_ggplot <- ggplot(boxplot_influenza, aes(x = disease, y = R, fill = Withdrawal)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept=1, color="darkred", linetype="dashed", linewidth =1) +
  labs(title = "",
       x = "Influenza Variants",
       y = "") +
  scale_fill_manual(values = c("Natural" = "#CC3311", "Smartwatch" = "#88CCEE")) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),        
        legend.title=element_blank(),
        legend.position = c(0.75, 0.9),
        text = element_text(size = 15),  # Adjust text size for readability
        plot.title = element_text(size = 15, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 15, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 15, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 15)  # Caption styling
  ) + scale_y_continuous(expand = c(0,0), breaks=seq(0,7, by=1), limits = c(0,7.2))


covid_ggplot + influenza_ggplot





###RESULTS REPORTING#### NOTE the values here are the same as in _cohort DFs
mean(ancestral_boxplot_n$R)
IQR(ancestral_boxplot_n$R)
quantile(ancestral_boxplot_n$R)
mean(ancestral_boxplot_sw$R)
IQR(ancestral_boxplot_sw$R)
quantile(ancestral_boxplot_sw$R)
(mean(ancestral_boxplot_n$R)-mean(ancestral_boxplot_sw$R))/mean(ancestral_boxplot_n$R) # Percentage reduction
t.test(ancestral_cohort$R, ancestral_cohort$R_sw, paired = TRUE)

mean(delta_boxplot_n$R)
IQR(delta_boxplot_n$R)
quantile(delta_boxplot_n$R)
mean(delta_boxplot_sw$R)
IQR(delta_boxplot_sw$R)
quantile(delta_boxplot_sw$R)
(mean(delta_boxplot_n$R)-mean(delta_boxplot_sw$R))/mean(delta_boxplot_n$R)
t.test(delta_cohort$R, delta_cohort$R_sw, paired = TRUE)

mean(omicron_boxplot_n$R)
IQR(omicron_boxplot_n$R)
quantile(omicron_boxplot_n$R)
mean(omicron_boxplot_sw$R)
IQR(omicron_boxplot_sw$R)
quantile(omicron_boxplot_sw$R)
(mean(omicron_boxplot_n$R)-mean(omicron_boxplot_sw$R))/mean(omicron_boxplot_n$R)
t.test(omicron_cohort$R, omicron_cohort$R_sw, paired = TRUE)

mean(pandemic_boxplot_n$R)
IQR(pandemic_boxplot_n$R)
quantile(pandemic_boxplot_n$R)
mean(pandemic_boxplot_sw$R)
IQR(pandemic_boxplot_sw$R)
quantile(pandemic_boxplot_sw$R)
(mean(pandemic_boxplot_n$R)-mean(pandemic_boxplot_sw$R))/mean(pandemic_boxplot_n$R)
t.test(pandemic_cohort$R, pandemic_cohort$R_sw, paired = TRUE)

mean(seasonal_boxplot_n$R)
IQR(seasonal_boxplot_n$R)
quantile(seasonal_boxplot_n$R)
mean(seasonal_boxplot_sw$R)
IQR(seasonal_boxplot_sw$R)
quantile(seasonal_boxplot_sw$R)
(mean(seasonal_boxplot_n$R)-mean(seasonal_boxplot_sw$R))/mean(seasonal_boxplot_n$R)
t.test(seasonal_cohort$R, seasonal_cohort$R_sw, paired = TRUE)


##FIGURE3##


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
  geom_point(aes(x=2,y=0.05,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=2,y=0.87,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=3,y=0.26,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=5,y=0.01,colour="Symptom onset"), size=4)+
  geom_point(aes(x=5,y=0.27,colour="Symptom onset"), size=4)+
  geom_point(aes(x=6,y=0.065,colour="Symptom onset"), size=4)+
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
        text = element_text(size = 15),  # Adjust text size for readability
        plot.title = element_text(size = 15, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 15, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 15, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 15)  # Caption styling
  )

uncertainty_influenza_ggplot <- ggplot() +
  geom_smooth(data = uncertainty_influenza, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Pandemic, linetype="Pandemic"), color = "black")+
  geom_smooth(data = uncertainty_influenza, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:25, y = Seasonal, linetype="Seasonal"),color = "black") +
  geom_point(aes(x=0.5,y=0.95,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=1,y=0.98,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=2,y=0.925,colour="Symptom onset"), size=4)+
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
        text = element_text(size = 15),  # Adjust text size for readability
        plot.title = element_text(size = 15, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 15, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 15, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 15)  # Caption styling
  )

uncertainty_covid_ggplot + uncertainty_influenza_ggplot











##FIGURE 4
set.seed(1)
calculate_needed_withdrawal <- function(df, withdrawal) {
  test <- df
  test$probs <- lapply(1:nrow(test), function(x) rep(NA, 251))
  for(i in 1:1000){
    if(test$Detected_by_smartwatch[i] == 1) {
      test$probs[[i]] <- test$VL_ecdf[[i]] + (1- test$VL_ecdf[[i]]) * (1-withdrawal) #IF DE
    } else if (test$Symptomatic[i] == 1) {
      test$probs[[i]] <-   rep(test$VL_at_natural_withdrawal[i], 251)   
    } else {
      test$probs[[i]] <- rep(1, 251)
    }                                     
    test$probs[[i]] <- test$probs[[i]] * test$R[i] / test$VL_at_natural_withdrawal[i]
  }
  test$probs_by_day <- lapply(test$probs, function(x) x[seq(1, length(x), by = 10)])
  
  test_uncertainty <- data.frame(t(as.data.frame(test$probs_by_day)))
  rownames(test_uncertainty) <- 1:1000
  colnames(test_uncertainty) <- 0:25
  test_uncertainty <- ifelse(test_uncertainty < 1, 1, 0)
  test_uncertainty <- colSums(test_uncertainty) / nrow(test_uncertainty)
  test_uncertainty <- data.frame(result = test_uncertainty)
  return(test_uncertainty)
}

####75% DF####
ancestral_withdrawal_requirement_75 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)
threshold <- 0.75

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(ancestral_cohort, j)
    if(temp$result[i] >= threshold) {
      ancestral_withdrawal_requirement_75[i] <- j
      break
    }
  }
}

delta_withdrawal_requirement_75 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(delta_cohort, j)
    if(temp$result[i] >= threshold) {
      delta_withdrawal_requirement_75[i] <- j
      break
    }
  }
}

omicron_withdrawal_requirement_75 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(omicron_cohort, j)
    if(temp$result[i] >= threshold) {
      omicron_withdrawal_requirement_75[i] <- j
      break
    }
  }
}

pandemic_withdrawal_requirement_75 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(pandemic_cohort, j)
    if(temp$result[i] >= threshold) {
      pandemic_withdrawal_requirement_75[i] <- j
      break
    }
  }
}


seasonal_withdrawal_requirement_75 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(seasonal_cohort, j)
    if(temp$result[i] >= threshold) {
      seasonal_withdrawal_requirement_75[i] <- j
      break
    }
  }
}


#
ancestral_withdrawal_requirement_75
delta_withdrawal_requirement_75
omicron_withdrawal_requirement_75
pandemic_withdrawal_requirement_75
seasonal_withdrawal_requirement_75



####50% DF####
ancestral_withdrawal_requirement_50 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)
threshold <- 0.50

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(ancestral_cohort, j)
    if(temp$result[i] >= threshold) {
      ancestral_withdrawal_requirement_50[i] <- j
      break
    }
  }
}

delta_withdrawal_requirement_50 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(delta_cohort, j)
    if(temp$result[i] >= threshold) {
      delta_withdrawal_requirement_50[i] <- j
      break
    }
  }
}

omicron_withdrawal_requirement_50 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(omicron_cohort, j)
    if(temp$result[i] >= threshold) {
      omicron_withdrawal_requirement_50[i] <- j
      break
    }
  }
}

pandemic_withdrawal_requirement_50 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(pandemic_cohort, j)
    if(temp$result[i] >= threshold) {
      pandemic_withdrawal_requirement_50[i] <- j
      break
    }
  }
}


seasonal_withdrawal_requirement_50 <- rep(0, 10)
withdrawals <- seq(0.4, 1, by=0.01)

for(i in 1:11) {
  for(j in withdrawals) {
    temp <- calculate_needed_withdrawal(seasonal_cohort, j)
    if(temp$result[i] >= threshold) {
      seasonal_withdrawal_requirement_50[i] <- j
      break
    }
  }
}


#
ancestral_withdrawal_requirement_50
delta_withdrawal_requirement_50
omicron_withdrawal_requirement_50
pandemic_withdrawal_requirement_50
seasonal_withdrawal_requirement_50


withdrawal_requirement_covids_50 <- data.frame(Ancestral = ancestral_withdrawal_requirement_50,
                                      Delta = delta_withdrawal_requirement_50,
                                      Omicron = omicron_withdrawal_requirement_50)

withdrawal_requirement_covids_75 <- data.frame(Ancestral = ancestral_withdrawal_requirement_75,
                                        Delta = delta_withdrawal_requirement_75,
                                        Omicron = omicron_withdrawal_requirement_75)


withdrawal_requirement_covids_50[withdrawal_requirement_covids_50 == 0] <- 2
withdrawal_requirement_covids_75[withdrawal_requirement_covids_75 == 0] <- 2





withdrawal_requirement_influenzas_50 <- data.frame(Pandemic = pandemic_withdrawal_requirement_50,
                                            Seasonal = seasonal_withdrawal_requirement_50)
withdrawal_requirement_influenzas_75 <- data.frame(Pandemic = pandemic_withdrawal_requirement_75,
                                                   Seasonal = seasonal_withdrawal_requirement_75)


withdrawal_requirement_influenzas_50[withdrawal_requirement_influenzas_50 == 0] <- 2
withdrawal_requirement_influenzas_75[withdrawal_requirement_influenzas_75 == 0] <- 2





######## FIGURE 4#########

withdrawal_covid_ggplot <- ggplot() +
  geom_smooth(data = withdrawal_requirement_covids_75, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Ancestral, linetype="Ancestral", color = "75% probability of elimination"))+
  geom_smooth(data = withdrawal_requirement_covids_75, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Delta, linetype="Delta",color = "75% probability of elimination")) +
  geom_smooth(data = withdrawal_requirement_covids_75, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Omicron, linetype="Omicron",color = "75% probability of elimination")) +
  geom_smooth(data = withdrawal_requirement_covids_50, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Ancestral, linetype="Ancestral", color = "50% probability of elimination"))+
  geom_smooth(data = withdrawal_requirement_covids_50, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Delta, linetype="Delta",color = "50% probability of elimination")) +
  geom_smooth(data = withdrawal_requirement_covids_50, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Omicron, linetype="Omicron",color = "50% probability of elimination")) +
  
  geom_point(aes(x=2,y=0.602,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=2,y=0.661,colour="Smartwatch detection"), size=4)+  
  
  geom_point(aes(x=2,y=0.899,colour="Smartwatch detection"), size=4)+  
  geom_point(aes(x=2,y=0.91,colour="Smartwatch detection"), size=4)+  
  
  geom_point(aes(x=3,y=0.837,colour="Smartwatch detection"), size=4)+  
  geom_point(aes(x=3,y=0.808,colour="Smartwatch detection"), size=4)+  
  

  geom_point(aes(x=5,y=0.743,colour="Symptom onset"), size=4)+
  geom_point(aes(x=5,y=0.811,colour="Symptom onset"), size=4)+
  
  scale_color_manual(values = c("Smartwatch detection" = "#88CCEE", "Symptom onset" = "#CC3311",
                                "50% probability of elimination" = "black", "75% probability of elimination" = "darkgray"),
                     labels = c("Smartwatch detection" = "Smartwatch detection", "Symptom onset" = "Symptom onset",
                                "50% probability of elimination" = "50% probability of elimination", "75% probability of elimination" = "75% probability of elimination")) +
  scale_x_continuous(expand= c(0,0), limits = c(0,8), breaks= c(0,2,4,6,8)) +
  scale_y_continuous(expand= c(0,0), limits = c(0.4, 1), labels = scales::percent)  +
  labs(title = "",
       x = "Time from exposure (days)",
       y = "Required reduction in contacts (%)") +
  guides(
    color = guide_legend(order = 2, title = "Point Types"),
    linetype = guide_legend(order = 1, title = "Line Types")
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.6),
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  )

withdrawal_covid_ggplot







withdrawal_influenza_ggplot <- ggplot() +
  geom_smooth(data = withdrawal_requirement_influenzas_75, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Pandemic, linetype="Pandemic", color = "75% probability of elimination"))+
  geom_smooth(data = withdrawal_requirement_influenzas_75, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Seasonal, linetype="Seasonal",color = "75% probability of elimination")) +
  geom_smooth(data = withdrawal_requirement_influenzas_50, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Pandemic, linetype="Pandemic", color = "50% probability of elimination"))+
  geom_smooth(data = withdrawal_requirement_influenzas_50, method= "glm", method.args = list(family = "binomial"), se = FALSE, aes(x = 0:9, y = Seasonal, linetype="Seasonal",color = "50% probability of elimination")) +
  
  geom_point(aes(x=0.5,y=0.551,colour="Smartwatch detection"), size=4)+
  geom_point(aes(x=0.5,y=0.65,colour="Smartwatch detection"), size=4)+  
  geom_point(aes(x=1,y=0.502,colour="Smartwatch detection"), size=4)+  
  geom_point(aes(x=1,y=0.557,colour="Smartwatch detection"), size=4)+  
  
  
  geom_point(aes(x=2,y=0.598,colour="Symptom onset"), size=4)+
  geom_point(aes(x=2,y=0.668,colour="Symptom onset"), size=4)+
  geom_point(aes(x=2,y=0.723,colour="Symptom onset"), size=4)+
  geom_point(aes(x=2,y=0.759,colour="Symptom onset"), size=4)+
  
  scale_color_manual(values = c("Smartwatch detection" = "#88CCEE", "Symptom onset" = "#CC3311",
                                "50% probability of elimination" = "black", "75% probability of elimination" = "darkgray"),
                     labels = c("Smartwatch detection" = "Smartwatch detection", "Symptom onset" = "Symptom onset",
                                "50% probability of elimination" = "50% probability of elimination", "75% probability of elimination" = "75% probability of elimination")) +
  scale_x_continuous(expand= c(0,0), limits = c(0,8), breaks= c(0,2,4,6,8)) +
  scale_y_continuous(expand= c(0,0), limits = c(0.4, 1), labels = scales::percent)  +  
  labs(title = "",
       x = "Time from exposure (days)",
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
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.6),
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  )

withdrawal_influenza_ggplot
withdrawal_covid_ggplot + withdrawal_influenza_ggplot





ggsave("240916_Figure4.jpeg")######### SENSITIVITY ANALYSES ############

##ANCESTRAL SENSITIVITY: DETECTION X COMPLIANCCE##
set.seed(123)
ancestral_detection <- data.frame(detection = seq(from = 10 , to = 70, by= 10))
ancestral_compliance <- data.frame(compliance = seq(from = 0.6 , to = 1, by= 0.1))
ancestral_combinations <- expand.grid(detection=ancestral_detection$detection , compliance=ancestral_compliance$compliance)

for(i in ancestral_detection$detection) { 
  for(j in ancestral_compliance$compliance) {
  ancestral_cohort$compliance <- rbinom(1000, 1, prob = j)
  temp <-  mapply(extract_value, i, ancestral_cohort$VL_ecdf)
  ancestral_cohort[[paste0("ECDF_at_day_", i)]] <- ifelse(ancestral_cohort$Detected_by_smartwatch == 1 & ancestral_cohort$compliance == 1, temp + (1-temp)*(1-ancestral_cohort$Withdrawal), ancestral_cohort$VL_at_natural_withdrawal)
  ancestral_combinations$ECDF[ancestral_combinations$detection == i & ancestral_combinations$compliance == j] <- mean(ancestral_cohort[[paste0("ECDF_at_day_", i)]])
  }
}

ancestral_combinations$R_original <- mean(ancestral_cohort$R)
ancestral_combinations$R_eff <- (ancestral_combinations$ECDF * ancestral_combinations$R_original) / mean(ancestral_cohort$VL_at_natural_withdrawal)
ancestral_combinations$R_reduced <- abs(ancestral_combinations$R_eff - ancestral_combinations$R_original)
ancestral_combinations$R_redc_perc <- ancestral_combinations$R_reduced / ancestral_combinations$R_original

ancestral_combinations$detection <- ancestral_combinations$detection/10
ancestral_combinations$R_redc_perc <- ancestral_combinations$R_redc_perc* 100

# Heatmap 
ancestral_heatmap <- ggplot(ancestral_combinations, aes(detection, compliance, z = R_redc_perc)) +
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_text_contour(aes(z = R_redc_perc), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 20, limits=c(0,60))+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Compliance to smartwatch (%)",
    fill = "% reduction in effective reproduction number"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  ) + ggtitle('Covid (Ancestral)')

ancestral_heatmap


##DELTA SENSITIVITY: DETECTION X COMPLIANCCE##
set.seed(123)
delta_detection <- data.frame(detection = seq(from = 10 , to = 70, by= 10))
delta_compliance <- data.frame(compliance = seq(from = 0.6 , to = 1, by= 0.1))
delta_combinations <- expand.grid(detection=delta_detection$detection , compliance=delta_compliance$compliance)

for(i in delta_detection$detection) { 
  for(j in delta_compliance$compliance) {
    delta_cohort$compliance <- rbinom(1000, 1, prob = j)
    temp <-  mapply(extract_value, i, delta_cohort$VL_ecdf)
    delta_cohort[[paste0("ECDF_at_day_", i)]] <- ifelse(delta_cohort$Detected_by_smartwatch == 1 & delta_cohort$compliance == 1, temp + (1-temp)*(1-delta_cohort$Withdrawal), delta_cohort$VL_at_natural_withdrawal)
    delta_combinations$ECDF[delta_combinations$detection == i & delta_combinations$compliance == j] <- mean(delta_cohort[[paste0("ECDF_at_day_", i)]])
  }
}

delta_combinations$R_original <- mean(delta_cohort$R)
delta_combinations$R_eff <- (delta_combinations$ECDF * delta_combinations$R_original) / mean(delta_cohort$VL_at_natural_withdrawal)
delta_combinations$R_reduced <- abs(delta_combinations$R_eff - delta_combinations$R_original)
delta_combinations$R_redc_perc <- delta_combinations$R_reduced / delta_combinations$R_original

delta_combinations$detection <- delta_combinations$detection/10
delta_combinations$R_redc_perc <- delta_combinations$R_redc_perc* 100

# Heatmap 
delta_heatmap <- ggplot(delta_combinations, aes(detection, compliance, z = R_redc_perc)) + 
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_text_contour(aes(z = R_redc_perc), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 20, limits=c(0,60))+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Compliance to smartwatch (%)",
    fill = "% reduction in effective reproduction number"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  )  + ggtitle('Covid (Delta)')

delta_heatmap


##OMICRON SENSITIVITY: DETECTION X COMPLIANCCE##
set.seed(123)
omicron_detection <- data.frame(detection = seq(from = 10 , to = 70, by= 10))
omicron_compliance <- data.frame(compliance = seq(from = 0.6 , to = 1, by= 0.05))
omicron_combinations <- expand.grid(detection=omicron_detection$detection , compliance=omicron_compliance$compliance)

for(i in omicron_detection$detection) { 
  for(j in omicron_compliance$compliance) {
    omicron_cohort$compliance <- rbinom(1000, 1, prob = j)
    temp <-  mapply(extract_value, i, omicron_cohort$VL_ecdf)
    omicron_cohort[[paste0("ECDF_at_day_", i)]] <- ifelse(omicron_cohort$Detected_by_smartwatch == 1 & omicron_cohort$compliance == 1, temp + (1-temp)*(1-omicron_cohort$Withdrawal), omicron_cohort$VL_at_natural_withdrawal)
    omicron_combinations$ECDF[omicron_combinations$detection == i & omicron_combinations$compliance == j] <- mean(omicron_cohort[[paste0("ECDF_at_day_", i)]])
  }
}

omicron_combinations$R_original <- mean(omicron_cohort$R)
omicron_combinations$R_eff <- (omicron_combinations$ECDF * omicron_combinations$R_original) / mean(omicron_cohort$VL_at_natural_withdrawal)
omicron_combinations$R_reduced <- abs(omicron_combinations$R_eff - omicron_combinations$R_original)
omicron_combinations$R_redc_perc <- omicron_combinations$R_reduced / omicron_combinations$R_original

omicron_combinations$detection <- omicron_combinations$detection/10
omicron_combinations$R_redc_perc <- omicron_combinations$R_redc_perc* 100

# Heatmap 
omicron_heatmap <- ggplot(omicron_combinations, aes(detection, compliance, z = R_redc_perc)) + 
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_text_contour(aes(z = R_redc_perc), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 20, limits=c(0,60))+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Compliance to smartwatch (%)",
    fill = "% reduction in effective reproduction number"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  )  + ggtitle('Covid (Omicron)')
omicron_heatmap 




##pandemic SENSITIVITY: DETECTION X COMPLIANCCE##
pandemic_detection <- data.frame(detection = seq(from = 10 , to = 50, by= 10))
pandemic_compliance <- data.frame(compliance = seq(from = 0.6 , to = 1, by= 0.05))
pandemic_combinations <- expand.grid(detection=pandemic_detection$detection , compliance=pandemic_compliance$compliance)

for(i in pandemic_detection$detection) { 
  for(j in pandemic_compliance$compliance) {
    pandemic_cohort$compliance <- rbinom(1000, 1, prob = j)
    temp <-  mapply(extract_value, i, pandemic_cohort$VL_ecdf)
    pandemic_cohort[[paste0("ECDF_at_day_", i)]] <- ifelse(pandemic_cohort$Detected_by_smartwatch == 1 & pandemic_cohort$compliance == 1, temp + (1-temp)*(1-pandemic_cohort$Withdrawal), pandemic_cohort$VL_at_natural_withdrawal)
    pandemic_combinations$ECDF[pandemic_combinations$detection == i & pandemic_combinations$compliance == j] <- mean(pandemic_cohort[[paste0("ECDF_at_day_", i)]])
  }
}

pandemic_combinations$R_original <- mean(pandemic_cohort$R)
pandemic_combinations$R_eff <- (pandemic_combinations$ECDF * pandemic_combinations$R_original) / mean(pandemic_cohort$VL_at_natural_withdrawal)
pandemic_combinations$R_reduced <- abs(pandemic_combinations$R_eff - pandemic_combinations$R_original)
pandemic_combinations$R_redc_perc <- pandemic_combinations$R_reduced / pandemic_combinations$R_original

pandemic_combinations$detection <- pandemic_combinations$detection/10
pandemic_combinations$R_redc_perc <- pandemic_combinations$R_redc_perc* 100

# Heatmap 
pandemic_heatmap <- ggplot(pandemic_combinations, aes(detection, compliance, z = R_redc_perc)) + 
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_text_contour(aes(z = R_redc_perc), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 20, limits=c(0,60))+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Compliance to smartwatch (%)",
    fill = "% reduction in effective reproduction number"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  )  + ggtitle('Pandemic Influenza')
pandemic_heatmap 



##seasonal SENSITIVITY: DETECTION X COMPLIANCCE##
seasonal_detection <- data.frame(detection = seq(from = 10 , to = 50, by= 10))
seasonal_compliance <- data.frame(compliance = seq(from = 0.6 , to = 1, by= 0.05))
seasonal_combinations <- expand.grid(detection=seasonal_detection$detection , compliance=seasonal_compliance$compliance)

for(i in seasonal_detection$detection) { 
  for(j in seasonal_compliance$compliance) {
    seasonal_cohort$compliance <- rbinom(1000, 1, prob = j)
    temp <-  mapply(extract_value, i, seasonal_cohort$VL_ecdf)
    seasonal_cohort[[paste0("ECDF_at_day_", i)]] <- ifelse(seasonal_cohort$Detected_by_smartwatch == 1 & seasonal_cohort$compliance == 1, temp + (1-temp)*(1-seasonal_cohort$Withdrawal), seasonal_cohort$VL_at_natural_withdrawal)
    seasonal_combinations$ECDF[seasonal_combinations$detection == i & seasonal_combinations$compliance == j] <- mean(seasonal_cohort[[paste0("ECDF_at_day_", i)]])
  }
}

seasonal_combinations$R_original <- mean(seasonal_cohort$R)
seasonal_combinations$R_eff <- (seasonal_combinations$ECDF * seasonal_combinations$R_original) / mean(seasonal_cohort$VL_at_natural_withdrawal)
seasonal_combinations$R_reduced <- abs(seasonal_combinations$R_eff - seasonal_combinations$R_original)
seasonal_combinations$R_redc_perc <- seasonal_combinations$R_reduced / seasonal_combinations$R_original

seasonal_combinations$detection <- seasonal_combinations$detection/10
seasonal_combinations$R_redc_perc <- seasonal_combinations$R_redc_perc* 100

# Heatmap 
seasonal_heatmap <- ggplot(seasonal_combinations, aes(detection, compliance, z = R_redc_perc)) + 
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_text_contour(aes(z = R_redc_perc), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 20, limits=c(0,60))+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Compliance to smartwatch (%)",
    fill = "% reduction in effective reproduction number"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),       
        text = element_text(size = 20),  # Adjust text size for readability
        plot.title = element_text(size = 20, face = "bold"),  # Title styling
        plot.subtitle = element_text(size = 20, face = "italic"),  # Subtitle styling
        axis.title = element_text(size = 20, face = "bold"),  # Axis title styling
        plot.caption = element_text(size = 20)  # Caption styling
  )  + ggtitle('Seasonal Influenza')
seasonal_heatmap 




(ancestral_heatmap|
delta_heatmap|
omicron_heatmap|
pandemic_heatmap|
seasonal_heatmap) + plot_layout(axis_titles = "collect")


















