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
library(gghighlight)
library(tidyverse)
library(shadowtext)


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


##WITHDRAWAL X DAY OF DETECTION##
#ANCESTRAL#
set.seed(123)
ancestral_detection <- data.frame(detection = seq(from = 10 , to = 70, by= 10))
ancestral_withdrawal <- data.frame(withdrawal = seq(from = 0.5 , to = 1, by= 0.01))
ancestral_combinations <- expand.grid(detection=ancestral_detection$detection , withdrawal=ancestral_withdrawal$withdrawal)

#ANCESTRAL#
ancestral_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = ancestral_rate_param, shape = ancestral_shape_param)), simplify = TRUE)
ancestral_cohort$VL_ecdf <- lapply(ancestral_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))

for(w in ancestral_withdrawal$withdrawal) {
  for(d in ancestral_detection$detection) {
    
ancestral_cohort$day_of_smartwatch_withdrawal <- d
ancestral_cohort$Withdrawal <- w
ancestral_cohort$time_at_natural_withdrawal <- round(ancestral_cohort$Day_of_symptom_onset + ancestral_cohort$Time_from_onset_to_testing,0)
#First calculate ECDF up to point of withdrawal
ancestral_cohort$VL_at_natural_withdrawal <- mapply(extract_value, ancestral_cohort$time_at_natural_withdrawal, ancestral_cohort$VL_ecdf)
ancestral_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, ancestral_cohort$day_of_smartwatch_withdrawal, ancestral_cohort$VL_ecdf)
#Then include the rest with the withdrawal coefficient, and without for those 
ancestral_cohort$VL_at_natural_withdrawal <- ifelse(ancestral_cohort$Symptomatic == 1, ancestral_cohort$VL_at_natural_withdrawal + (1-ancestral_cohort$VL_at_natural_withdrawal)*(1-ancestral_cohort$Withdrawal), 1)
ancestral_cohort$VL_at_smartwatch_withdrawal <- ifelse(ancestral_cohort$Detected_by_smartwatch == 1, ancestral_cohort$VL_at_smartwatch_withdrawal + (1-ancestral_cohort$VL_at_smartwatch_withdrawal)*(1-ancestral_cohort$Withdrawal), ancestral_cohort$VL_at_natural_withdrawal) 
####R's with SMARTWATCH
ancestral_cohort$R_sw <- (ancestral_cohort$VL_at_smartwatch_withdrawal * ancestral_cohort$R)/ ancestral_cohort$VL_at_natural_withdrawal
ancestral_combinations$R_new[ancestral_combinations$detection == d & ancestral_combinations$withdrawal == w] <- mean(ancestral_cohort$R_sw)
}}

ancestral_combinations$R_original <- mean(ancestral_cohort$R)
ancestral_combinations$R_reduced <- abs(ancestral_combinations$R_new - ancestral_combinations$R_original)
ancestral_combinations$R_redc_perc <- ancestral_combinations$R_reduced / ancestral_combinations$R_original

ancestral_combinations$detection <- ancestral_combinations$detection/10
ancestral_combinations$R_redc_perc <- ancestral_combinations$R_redc_perc* 100


# Heatmap 
ancestral_heatmap <- ggplot(ancestral_combinations, aes(detection, withdrawal*100, z = R_redc_perc)) +
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_contour(aes(z = R_redc_perc), breaks = 60, colour = "black", linewidth = 2) +
  geom_text_contour(aes(z = R_redc_perc), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 30)+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_continuous(expand = c(0, 0), n.breaks = 7)+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Effective reduction of social contacts (%)",
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




##WITHDRAWAL X DAY OF DETECTION##
#DELTA#
set.seed(123)
delta_detection <- data.frame(detection = seq(from = 10 , to = 70, by= 10))
delta_withdrawal <- data.frame(withdrawal = seq(from = 0.5 , to = 1, by= 0.01))
delta_combinations <- expand.grid(detection=delta_detection$detection , withdrawal=delta_withdrawal$withdrawal)

#delta#
delta_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = delta_rate_param, shape = delta_shape_param)), simplify = TRUE)
delta_cohort$VL_ecdf <- lapply(delta_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))

for(w in delta_withdrawal$withdrawal) {
  for(d in delta_detection$detection) {
    
    delta_cohort$day_of_smartwatch_withdrawal <- d
    delta_cohort$Withdrawal <- w
    delta_cohort$time_at_natural_withdrawal <- round(delta_cohort$Day_of_symptom_onset + delta_cohort$Time_from_onset_to_testing,0)
    #First calculate ECDF up to point of withdrawal
    delta_cohort$VL_at_natural_withdrawal <- mapply(extract_value, delta_cohort$time_at_natural_withdrawal, delta_cohort$VL_ecdf)
    delta_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, delta_cohort$day_of_smartwatch_withdrawal, delta_cohort$VL_ecdf)
    #Then include the rest with the withdrawal coefficient, and without for those 
    delta_cohort$VL_at_natural_withdrawal <- ifelse(delta_cohort$Symptomatic == 1, delta_cohort$VL_at_natural_withdrawal + (1-delta_cohort$VL_at_natural_withdrawal)*(1-delta_cohort$Withdrawal), 1)
    delta_cohort$VL_at_smartwatch_withdrawal <- ifelse(delta_cohort$Detected_by_smartwatch == 1, delta_cohort$VL_at_smartwatch_withdrawal + (1-delta_cohort$VL_at_smartwatch_withdrawal)*(1-delta_cohort$Withdrawal), delta_cohort$VL_at_natural_withdrawal) 
    ####R's with SMARTWATCH
    delta_cohort$R_sw <- (delta_cohort$VL_at_smartwatch_withdrawal * delta_cohort$R)/ delta_cohort$VL_at_natural_withdrawal
    delta_combinations$R_new[delta_combinations$detection == d & delta_combinations$withdrawal == w] <- mean(delta_cohort$R_sw)
  }}

delta_combinations$R_original <- mean(delta_cohort$R)
delta_combinations$R_reduced <- abs(delta_combinations$R_new - delta_combinations$R_original)
delta_combinations$R_redc_perc <- delta_combinations$R_reduced / delta_combinations$R_original

delta_combinations$detection <- delta_combinations$detection/10
delta_combinations$R_redc_perc <- delta_combinations$R_redc_perc* 100


# Heatmap 
delta_heatmap <- ggplot(delta_combinations, aes(detection, withdrawal*100, z = R_redc_perc)) +
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_contour(aes(z = R_redc_perc), breaks = 35, colour = "black", linewidth = 2) +
  geom_text_contour(aes(z = R_redc_perc),  breaks = c(10,20,30,35,40,50,60), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 30)+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_continuous(expand = c(0, 0), n.breaks = 7)+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Effective reduction of social contacts (%)",
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
  ) + ggtitle('Covid (Delta)')

delta_heatmap






##WITHDRAWAL X DAY OF DETECTION##
#OMICRON#
set.seed(123)
omicron_detection <- data.frame(detection = seq(from = 10 , to = 70, by= 10))
omicron_withdrawal <- data.frame(withdrawal = seq(from = 0.5 , to = 1, by= 0.01))
omicron_combinations <- expand.grid(detection=omicron_detection$detection , withdrawal=omicron_withdrawal$withdrawal)

#omicron#
omicron_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = omicron_rate_param, shape = omicron_shape_param)), simplify = TRUE)
omicron_cohort$VL_ecdf <- lapply(omicron_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))

for(w in omicron_withdrawal$withdrawal) {
  for(d in omicron_detection$detection) {
    
    omicron_cohort$day_of_smartwatch_withdrawal <- d
    omicron_cohort$Withdrawal <- w
    omicron_cohort$time_at_natural_withdrawal <- round(omicron_cohort$Day_of_symptom_onset + omicron_cohort$Time_from_onset_to_testing,0)
    #First calculate ECDF up to point of withdrawal
    omicron_cohort$VL_at_natural_withdrawal <- mapply(extract_value, omicron_cohort$time_at_natural_withdrawal, omicron_cohort$VL_ecdf)
    omicron_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, omicron_cohort$day_of_smartwatch_withdrawal, omicron_cohort$VL_ecdf)
    #Then include the rest with the withdrawal coefficient, and without for those 
    omicron_cohort$VL_at_natural_withdrawal <- ifelse(omicron_cohort$Symptomatic == 1, omicron_cohort$VL_at_natural_withdrawal + (1-omicron_cohort$VL_at_natural_withdrawal)*(1-omicron_cohort$Withdrawal), 1)
    omicron_cohort$VL_at_smartwatch_withdrawal <- ifelse(omicron_cohort$Detected_by_smartwatch == 1, omicron_cohort$VL_at_smartwatch_withdrawal + (1-omicron_cohort$VL_at_smartwatch_withdrawal)*(1-omicron_cohort$Withdrawal), omicron_cohort$VL_at_natural_withdrawal) 
    ####R's with SMARTWATCH
    omicron_cohort$R_sw <- (omicron_cohort$VL_at_smartwatch_withdrawal * omicron_cohort$R)/ omicron_cohort$VL_at_natural_withdrawal
    omicron_combinations$R_new[omicron_combinations$detection == d & omicron_combinations$withdrawal == w] <- mean(omicron_cohort$R_sw)
  }}

omicron_combinations$R_original <- mean(omicron_cohort$R)
omicron_combinations$R_reduced <- abs(omicron_combinations$R_new - omicron_combinations$R_original)
omicron_combinations$R_redc_perc <- omicron_combinations$R_reduced / omicron_combinations$R_original

omicron_combinations$detection <- omicron_combinations$detection/10
omicron_combinations$R_redc_perc <- omicron_combinations$R_redc_perc* 100


# Heatmap 
omicron_heatmap <- ggplot(omicron_combinations, aes(detection, withdrawal*100, z = R_redc_perc)) +
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_contour(aes(z = R_redc_perc), breaks = 75, colour = "black", linewidth = 2) +
  geom_text_contour(aes(z = R_redc_perc),breaks = c(10,20,30,40,50,60,75), stroke = 0.15, rotate=FALSE, skip=0, size=7) +
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 30)+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_continuous(expand = c(0, 0), n.breaks = 7)+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Effective reduction of social contacts (%)",
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
  ) + ggtitle('Covid (Omicron)')

omicron_heatmap






##WITHDRAWAL X DAY OF DETECTION##
#PANDEMIC#
set.seed(123)
pandemic_detection <- data.frame(detection = seq(from = 10 , to = 50, by= 10))
pandemic_withdrawal <- data.frame(withdrawal = seq(from = 0.5 , to = 1, by= 0.01))
pandemic_combinations <- expand.grid(detection=pandemic_detection$detection , withdrawal=pandemic_withdrawal$withdrawal)

#pandemic#
pandemic_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = pandemic_rate_param, shape = pandemic_shape_param)), simplify = TRUE)
pandemic_cohort$VL[pandemic_cohort$Symptomatic == 0] <- replicate(n = length(pandemic_cohort$VL[pandemic_cohort$Symptomatic == 0]), list(rgamma(200, rate = pandemic_asym_rate_param, shape = pandemic_asym_shape_param)), simplify = TRUE)

pandemic_cohort$VL_ecdf <- lapply(pandemic_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
pandemic_cohort$VL_ecdf[pandemic_cohort$Symptomatic == 0] <- lapply(pandemic_cohort$VL_ecdf[pandemic_cohort$Symptomatic == 0], function(vl) vl * (2.7 / 4.7)) #Weights from Mean viral load from Table 3 of Dennis et al., 2016


for(w in pandemic_withdrawal$withdrawal) {
  for(d in pandemic_detection$detection) {
    
    pandemic_cohort$day_of_smartwatch_withdrawal <- d
    pandemic_cohort$Withdrawal <- w
    pandemic_cohort$time_at_natural_withdrawal <- round(pandemic_cohort$Day_of_symptom_onset + pandemic_cohort$Time_from_onset_to_testing,0)
    #First calculate ECDF up to point of withdrawal
    pandemic_cohort$VL_at_natural_withdrawal <- mapply(extract_value, pandemic_cohort$time_at_natural_withdrawal, pandemic_cohort$VL_ecdf)
    pandemic_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, pandemic_cohort$day_of_smartwatch_withdrawal, pandemic_cohort$VL_ecdf)
    #Then include the rest with the withdrawal coefficient, and without for those 
    pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1] <- pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1] + (1-pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1])*(1-pandemic_cohort$Withdrawal[pandemic_cohort$Symptomatic == 1])
    pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 0] <- sapply(pandemic_cohort$VL_ecdf[pandemic_cohort$Symptomatic == 0], function(x) x[length(x)])
    
    pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 0] <- pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 0]
    pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 0] <- pandemic_cohort$VL_at_natural_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 0]
    pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1] <- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1] + (1- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1]) * (1-pandemic_cohort$Withdrawal[pandemic_cohort$Symptomatic == 1 & pandemic_cohort$Detected_by_smartwatch == 1])
    pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1] <- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1] + ((2.7 / 4.7)- pandemic_cohort$VL_at_smartwatch_withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1]) * (1-pandemic_cohort$Withdrawal[pandemic_cohort$Symptomatic == 0 & pandemic_cohort$Detected_by_smartwatch == 1])
    ####R's with SMARTWATCH
    pandemic_cohort$R_sw <- (pandemic_cohort$VL_at_smartwatch_withdrawal * pandemic_cohort$R)/ pandemic_cohort$VL_at_natural_withdrawal
    pandemic_combinations$R_new[pandemic_combinations$detection == d & pandemic_combinations$withdrawal == w] <- mean(pandemic_cohort$R_sw)
  }}

pandemic_combinations$R_original <- mean(pandemic_cohort$R)
pandemic_combinations$R_reduced <- abs(pandemic_combinations$R_new - pandemic_combinations$R_original)
pandemic_combinations$R_redc_perc <- pandemic_combinations$R_reduced / pandemic_combinations$R_original

pandemic_combinations$detection <- pandemic_combinations$detection/10
pandemic_combinations$R_redc_perc <- pandemic_combinations$R_redc_perc* 100


# Heatmap 
pandemic_heatmap <- ggplot(pandemic_combinations, aes(detection, withdrawal*100, z = R_redc_perc)) +
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_contour(aes(z = R_redc_perc), breaks = 35, colour = "black", linewidth = 2) +
  geom_text_contour(aes(z = R_redc_perc),  breaks = c(10,20,30,35,40,50), stroke = 0.15, rotate=FALSE, skip=0, size=7)+
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 30)+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_continuous(expand = c(0, 0), n.breaks = 5)+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Effective reduction of social contacts (%)",
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
  ) + ggtitle('Influenza (Pandemic)')

pandemic_heatmap





##WITHDRAWAL X DAY OF DETECTION##
#SEASONAL#
set.seed(123)
seasonal_detection <- data.frame(detection = seq(from = 10 , to = 50, by= 10))
seasonal_withdrawal <- data.frame(withdrawal = seq(from = 0.5 , to = 1, by= 0.01))
seasonal_combinations <- expand.grid(detection=seasonal_detection$detection , withdrawal=seasonal_withdrawal$withdrawal)

#seasonal#
seasonal_cohort$VL <- replicate(n = 1000, list(rgamma(200, rate = seasonal_rate_param, shape = seasonal_shape_param)), simplify = TRUE)
seasonal_cohort$VL[seasonal_cohort$Symptomatic == 0] <- replicate(n = length(seasonal_cohort$VL[seasonal_cohort$Symptomatic == 0]), list(rgamma(200, rate = seasonal_asym_rate_param, shape = seasonal_asym_shape_param)), simplify = TRUE)

seasonal_cohort$VL_ecdf <- lapply(seasonal_cohort$VL, function(vl) ecdf(vl)(seq(0,25, by=0.1)))
seasonal_cohort$VL_ecdf[seasonal_cohort$Symptomatic == 0] <- lapply(seasonal_cohort$VL_ecdf[seasonal_cohort$Symptomatic == 0], function(vl) vl * (4.0 / 5.6)) #Weights from Mean viral load from Table 3 of Dennis et al., 2016


for(w in seasonal_withdrawal$withdrawal) {
  for(d in seasonal_detection$detection) {
    
    seasonal_cohort$day_of_smartwatch_withdrawal <- d
    seasonal_cohort$Withdrawal <- w
    seasonal_cohort$time_at_natural_withdrawal <- round(seasonal_cohort$Day_of_symptom_onset + seasonal_cohort$Time_from_onset_to_testing,0)
    #First calculate ECDF up to point of withdrawal
    seasonal_cohort$VL_at_natural_withdrawal <- mapply(extract_value, seasonal_cohort$time_at_natural_withdrawal, seasonal_cohort$VL_ecdf)
    seasonal_cohort$VL_at_smartwatch_withdrawal <- mapply(extract_value, seasonal_cohort$day_of_smartwatch_withdrawal, seasonal_cohort$VL_ecdf)
    #Then include the rest with the withdrawal coefficient, and without for those 
    seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1] <- seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1] + (1-seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1])*(1-seasonal_cohort$Withdrawal[seasonal_cohort$Symptomatic == 1])
    seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 0] <- sapply(seasonal_cohort$VL_ecdf[seasonal_cohort$Symptomatic == 0], function(x) x[length(x)])
    
    seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 0] <- seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 0]
    seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 0] <- seasonal_cohort$VL_at_natural_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 0]
    seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1] <- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1] + (1- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1]) * (1-seasonal_cohort$Withdrawal[seasonal_cohort$Symptomatic == 1 & seasonal_cohort$Detected_by_smartwatch == 1])
    seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1] <- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1] + ((4.0 / 5.6)- seasonal_cohort$VL_at_smartwatch_withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1]) * (1-seasonal_cohort$Withdrawal[seasonal_cohort$Symptomatic == 0 & seasonal_cohort$Detected_by_smartwatch == 1])
    
    ####R's with SMARTWATCH
    seasonal_cohort$R_sw <- (seasonal_cohort$VL_at_smartwatch_withdrawal * seasonal_cohort$R)/ seasonal_cohort$VL_at_natural_withdrawal
    seasonal_combinations$R_new[seasonal_combinations$detection == d & seasonal_combinations$withdrawal == w] <- mean(seasonal_cohort$R_sw)
  }}

seasonal_combinations$R_original <- mean(seasonal_cohort$R)
seasonal_combinations$R_reduced <- abs(seasonal_combinations$R_new - seasonal_combinations$R_original)
seasonal_combinations$R_redc_perc <- seasonal_combinations$R_reduced / seasonal_combinations$R_original

seasonal_combinations$detection <- seasonal_combinations$detection/10
seasonal_combinations$R_redc_perc <- seasonal_combinations$R_redc_perc* 100


# Heatmap 
seasonal_heatmap <- ggplot(seasonal_combinations, aes(detection, withdrawal*100, z = R_redc_perc)) +
  geom_raster(aes(fill=R_redc_perc), interpolate = TRUE)  +
  geom_contour(colour="black") +
  geom_contour(aes(z = R_redc_perc), breaks = 25, colour = "black", linewidth = 2) +
  geom_text_contour(aes(z = R_redc_perc),  breaks = c(10,20,25,30,40,50), stroke = 0.15, rotate=FALSE, skip=0, size=7)+
  scale_fill_gradient2('pi0', low = "#CC3311", mid = "#AA8080", high = "#88CCEE", midpoint = 30)+ ##ADJUST COLOR
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_continuous(expand = c(0, 0), n.breaks = 5)+
  labs(
    x = "Time from exposure to smartwatch detection (days)",
    y = "Effective reduction of social contacts (%)",
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
  ) + ggtitle('Influenza (Seasonal)')

seasonal_heatmap

Figure5 <- (ancestral_heatmap|
    delta_heatmap|
    omicron_heatmap|
    pandemic_heatmap|
    seasonal_heatmap) + plot_layout(axis_titles = "collect")
Figure5
#ggsave("240930_Figure5.jpeg")