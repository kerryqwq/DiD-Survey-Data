library(tidyverse)
library(knitr)
library(nnet)

source("Code/bootfunction.R")

real_data = read.csv("Data/real_data.csv", sep = ",", header = TRUE)

fit <- multinom(comb_group_label ~ sex + age + bmi + race4, data = real_data, trace=FALSE) 

real_data$psvalue <- predict(fit, type = "probs") # find the propensity scores
real_data$psvalue <- ifelse(real_data$psvalue < 0.01, 0.01, ifelse(real_data$psvalue > 0.99, 0.99, real_data$psvalue))

real_data <- real_data %>%
  mutate(ps_weight = case_when( # calculate the propensity score weights by groups
    comb_group_label == 1 ~ psvalue[,"1"] / psvalue[,"1"], 
    comb_group_label == 2 ~ psvalue[,"1"] / psvalue[,"2"],
    comb_group_label == 3 ~ psvalue[,"1"] / psvalue[,"3"],
    comb_group_label == 4 ~ psvalue[,"1"] / psvalue[,"4"] 
  )) # calculate the final weights

# calculate the target estimate
real_data$final_weight = real_data$ps_weight * real_data$weight

# DiD IPW estimator with the final weights
did_1 <- with(real_data, sum(soda_usage[comb_group_label == 2]*final_weight[comb_group_label == 2])/ (sum(weight[comb_group_label == 1])) -
                sum(soda_usage[comb_group_label == 1]*final_weight[comb_group_label == 1])/(sum(weight[comb_group_label == 1]))   -
                sum(soda_usage[comb_group_label == 4]*final_weight[comb_group_label == 4])/(sum(weight[comb_group_label == 1]))  +
                sum(soda_usage[comb_group_label == 3]*final_weight[comb_group_label == 3])/(sum(weight[comb_group_label == 1])) )

# DiD IPW estimator with the propensity score weights
did_2 <- with(real_data, sum(soda_usage[comb_group_label == 2]*ps_weight[comb_group_label == 2])/(sum(comb_group_label == 1)) -
                sum(soda_usage[comb_group_label == 1]*ps_weight[comb_group_label == 1])/(sum(comb_group_label == 1))   -
                sum(soda_usage[comb_group_label == 4]*ps_weight[comb_group_label == 4])/(sum(comb_group_label == 1))  +
                sum(soda_usage[comb_group_label == 3]*ps_weight[comb_group_label == 3])/(sum(comb_group_label == 1)) )

# DiD IPW estimator with the sampling weights
did_3 <- with(real_data, sum(soda_usage[comb_group_label == 2]*weight[comb_group_label == 2])/ (sum(weight[comb_group_label == 2])) -
                sum(soda_usage[comb_group_label == 1]*weight[comb_group_label == 1])/(sum(weight[comb_group_label == 1]))   -
                sum(soda_usage[comb_group_label == 4]*weight[comb_group_label == 4])/(sum(weight[comb_group_label == 4]))  +
                sum(soda_usage[comb_group_label == 3]*weight[comb_group_label == 3])/(sum(weight[comb_group_label == 3])))


# Perform the bootstrap resampling with replacement with 100 replicates
bootstrap_results = matrix(NA, nboot, 3)
for(r in 1:nboot){
  boot_dat <- real_data[sample(1:nrow(real_data), nrow(real_data), replace = TRUE),]
  bootstrap_results[r,] <- did_functions_real(boot_dat)
}

did_1; c(did_1 - 1.96*sd(bootstrap_results[,1]), did_1 + 1.96*sd(bootstrap_results[,1]))
did_2; c(did_2 - 1.96*sd(bootstrap_results[,2]), did_2 + 1.96*sd(bootstrap_results[,2]))
did_3; c(did_3 - 1.96*sd(bootstrap_results[,3]), did_3 + 1.96*sd(bootstrap_results[,3]))                       
