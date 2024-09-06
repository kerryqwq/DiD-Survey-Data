library(tidyverse)
library(knitr)
library(nnet)

source("Code/bootfunction.R")

N <- 2000 # population number
max.time <- 1 # maximum time period: two period
trt.time <- 1 # the second period is the post-trt period
nsims <- 100 # simulation times

t0 <- 0 # for Y0_0
t1 <- 1 # for Y0_1
num_S <- 0

ipw_estimates_1 <- rep(NA, nsims)
ipw_estimates_2 <- rep(NA, nsims)
ipw_estimates_3 <- rep(NA, nsims)

bias1_1 <- rep(NA, nsims)
bias1_2 <- rep(NA, nsims)
bias1_3 <- rep(NA, nsims)

bias2_1 <- rep(NA, nsims)
bias2_2 <- rep(NA, nsims)
bias2_3 <- rep(NA, nsims)

bias3_1 <- rep(NA, nsims)
bias3_2 <- rep(NA, nsims)
bias3_3 <- rep(NA, nsims)

sd_estimates_1 <- 0
sd_estimates_2 <- 0
sd_estimates_3 <- 0

sd_did_1 <- rep(NA, nsims)
sd_did_2 <- rep(NA, nsims)
sd_did_3 <- rep(NA, nsims)

sd_boot_1 <- 0
sd_boot_2 <- 0
sd_boot_3 <- 0


patt_trt = patt_g1 = satt_g1 = rep(NA, nsims)


gamma0 = c(-1, 0.5); gamma1 = c(1, -0.5)
eta0 = c(0.5,-0.5); eta1 = c(0.5,-0.5)

for (i in 1:nsims) {
  set.seed(i)
 
  # create dataset for pre-treatment period
  orig_dat1 <- data.frame(id = 1:N, tp = 0) %>% arrange(id, tp) %>% group_by(id) %>%
    mutate(int = rnorm(1, 0, sd = 0.25), # random intercept
           x1 = rnorm(1, mean = 0.50, sd = 1), #time-invariant covariate
           trt = rbinom(1, 1, plogis(gamma0[1] + gamma0[2]*x1)), # treated units
           post = I(tp >= trt.time), # indicator of post-treatment period
           treated = I(post == 1 & trt == 1) # interaction
    )%>% 
    ungroup()
  
  orig_dat1 <- orig_dat1 %>%
    group_by(id) %>%
    mutate(S = rbinom(1, 1, plogis(eta0[1] + eta0[2]*x1)), # sampling indicator based on the sampling probability
           y0_0 = int + x1 + 0.5*x1^2 + trt + rnorm(1), # potential outcome under control in the pre period
           y0_1 = int + x1 + 0.5*x1^2 + trt + 0.5 + rnorm(1),# potential outcome under control in the post period
           y1_0 = y0_0 + 1 - 0.5*x1, # potential outcome under treatment in the pre period
           y1_1 = y0_1 + 1 - 0.5*x1, # potential outcome under treatment in the post period
           treated_label = case_when(
             trt ==1 & post == 0 ~ 1, # label as group 1 if the treated unit is in pre-treatment period
             trt ==1 & post == 1 ~ 2, # label as group 2 if the treated unit is in post-treatment period
             trt ==0 & post == 0 ~ 3, # label as group 3 if the control unit is in pre-treatment period
             trt ==0 & post == 1 ~ 4  # label as group 4 if the control unit is in post-treatment period
           ),
           y = case_when(
             treated_label == 1 ~ y0_0, # the observed outcome of group 1 equivalents to y0_0
             treated_label == 2 ~ y1_1, # the observed outcome of group 2 equivalents to y1_1
             treated_label == 3 ~ y0_0, # the observed outcome of group 3 equivalents to y0_0
             treated_label == 4 ~ y0_1, # the observed outcome of group 4 equivalents to y0_1
           ))
  
  dat1 <- orig_dat1 %>%
    filter(S==1) # sample the data from the population  
  
  # ------------------------------------------------------------------------------------------
  # create dataset for post-treatment period
  
  orig_dat2 <- data.frame(id = (N+1):(2*N), tp = 1) %>% arrange(id, tp) %>% group_by(id)  %>%
    mutate(int = rnorm(1, 0, sd = 0.25), # random intercept
           x1 = rnorm(1, mean = 0.50, sd = 1), #time-invariant covariate
           trt = rbinom(1, 1, plogis(gamma1[1] + gamma1[2]*x1)), # treated units
           post = I(tp >= trt.time), # indicator of post-treatment period
           treated = I(post == 1 & trt == 1) # interaction
    )%>%
    ungroup()
  
  
  orig_dat2 <- orig_dat2 %>%
    group_by(id) %>%
    mutate(S = rbinom(1, 1, plogis(eta1[1] + eta1[2]*x1)), # sampling indicator based on the sampling probability
           y0_0 = int + x1 + 0.5*x1^2 + trt + rnorm(1), # potential outcome under control in the pre period
           y0_1 = int + x1 + 0.5*x1^2 + trt + 0.5 + rnorm(1),# potential outcome under control in the post period
           y1_0 = y0_0 + 1 - 0.5*x1, # potential outcome under treatment in the pre period
           y1_1 = y0_1 + 1 - 0.5*x1, # potential outcome under treatment in the post period
           treated_label = case_when(
             trt ==1 & post == 0 ~ 1, # label as group 1 if the treated unit is in pre-treatment period
             trt ==1 & post == 1 ~ 2, # label as group 2 if the treated unit is in post-treatment period
             trt ==0 & post == 0 ~ 3, # label as group 3 if the control unit is in pre-treatment period
             trt ==0 & post == 1 ~ 4  # label as group 4 if the control unit is in post-treatment period
           ),
           y = case_when(
             treated_label == 1 ~ y0_0, # the observed outcome of group 1 equivalents to y0_0
             treated_label == 2 ~ y1_1, # the observed outcome of group 2 equivalents to y1_1
             treated_label == 3 ~ y0_0, # the observed outcome of group 3 equivalents to y0_0
             treated_label == 4 ~ y0_1, # the observed outcome of group 4 equivalents to y0_1
           ))
  
  dat2 <- orig_dat2 %>%
    filter(S==1) # sample the data from the population 
  
  
  orig_dat <- rbind(orig_dat1, orig_dat2) %>% # combine the original dataset from pre and post
    arrange(id, tp)
  
  dat <- rbind(dat1, dat2) %>% # combine the sample dataset from pre and post
    arrange(id, tp)
  
 
  patt_trt[i] <- mean(orig_dat$y1_1[orig_dat$treated_label == 2]) - mean(orig_dat$y0_0[orig_dat$treated_label == 1]) - 
    (mean(orig_dat$y0_1[orig_dat$treated_label == 4]) - mean(orig_dat$y0_0[orig_dat$treated_label == 3]))
  
  patt_g1[i] <- mean(orig_dat$y1_1[orig_dat$treated_label ==1]) -
    mean(orig_dat$y0_1[orig_dat$treated_label ==1])
  
  satt_g1[i] <- mean(orig_dat$y1_1[orig_dat$treated_label ==1 & orig_dat$S ==1]) -
    mean(orig_dat$y0_1[orig_dat$treated_label ==1 & orig_dat$S ==1])
  
  # estimate propensity scores by using multinomial logistic regression
  fit <- multinom(treated_label ~ x1, data = dat, trace=FALSE) 
  
  dat$psvalue <- predict(fit, type = "probs") # find the propensity scores
  dat$psvalue <- ifelse(dat$psvalue < 0.01, 0.01, ifelse(dat$psvalue > 0.99, 0.99, dat$psvalue))
  
  dat <- dat %>%
    mutate(ps_weight = case_when( # calculate the propensity score weights by groups
      treated_label == 1 ~ psvalue[,"1"] / psvalue[,"1"], 
      treated_label == 2 ~ psvalue[,"1"] / psvalue[,"2"],
      treated_label == 3 ~ psvalue[,"1"] / psvalue[,"3"],
      treated_label == 4 ~ psvalue[,"1"] / psvalue[,"4"] 
    ),
    samp_prob = case_when( 
      tp == 0 ~ plogis(eta0[1] + eta0[2]*x1),
      tp == 1 ~ plogis(eta1[1] + eta1[2]*x1)
    ),
    samp_weight = 1 / samp_prob # calculate the sampling weights
    ) 
  
  dat$final_weight = dat$ps_weight * dat$samp_weight   # calculate the final weights
  
  # 
  did_1 <- with(dat, sum(y[treated_label == 2]*final_weight[treated_label == 2])/ (sum(samp_weight[treated_label == 1])) -
                  sum(y[treated_label == 1]*final_weight[treated_label == 1])/(sum(samp_weight[treated_label == 1]))   -
                  sum(y[treated_label == 4]*final_weight[treated_label == 4])/(sum(samp_weight[treated_label == 1]))  +
                  sum(y[treated_label == 3]*final_weight[treated_label == 3])/(sum(samp_weight[treated_label == 1])) )
  
  # 
  did_2 <- with(dat, sum(y[treated_label == 2]*ps_weight[treated_label == 2])/(sum(treated_label == 1)) -
                  sum(y[treated_label == 1]*ps_weight[treated_label == 1])/(sum(treated_label == 1))   -
                  sum(y[treated_label == 4]*ps_weight[treated_label == 4])/(sum(treated_label == 1))  +
                  sum(y[treated_label == 3]*ps_weight[treated_label == 3])/(sum(treated_label == 1)) )
  
  # 
  did_3 <- with(dat, sum(y[treated_label == 2]*samp_weight[treated_label == 2])/ (sum(samp_weight[treated_label == 2])) -
                  sum(y[treated_label == 1]*samp_weight[treated_label == 1])/(sum(samp_weight[treated_label == 1]))   -
                  sum(y[treated_label == 4]*samp_weight[treated_label == 4])/(sum(samp_weight[treated_label == 4]))  +
                  sum(y[treated_label == 3]*samp_weight[treated_label == 3])/(sum(samp_weight[treated_label == 3])))
  
  
  bootstrap_results = matrix(NA, nboot, 3)
  for(r in 1:nboot){
    boot_dat <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE),]
    bootstrap_results[r,] <- did_functions(boot_dat)
  }
  
  # Calculate the bootstrap standard deviation across 100 replicates 
  sd_did_1[i] <- sd(bootstrap_results[,1])
  sd_did_2[i] <- sd(bootstrap_results[,2])
  sd_did_3[i] <- sd(bootstrap_results[,3])
  
  # calculate the IPW estimates
  ipw_estimates_1[i] <- did_1
  ipw_estimates_2[i] <- did_2
  ipw_estimates_3[i] <- did_3
  
  bias1_1[i] <- did_1 - patt_g1[i] 
  bias1_2[i] <- did_1 - satt_g1[i] 
  bias1_3[i] <- did_1 - patt_trt[i] 
  
  bias2_1[i] <- did_2 - patt_g1[i] 
  bias2_2[i] <- did_2 - satt_g1[i] 
  bias2_3[i] <- did_2 - patt_trt[i] 
  
  bias3_1[i] <- did_3 - patt_g1[i] 
  bias3_2[i] <- did_3 - satt_g1[i] 
  bias3_3[i] <- did_3 - patt_trt[i] 
  
  # calculate the number of units in the sample
  num_S[i] <- length(unique(dat$id))
  
}

# calculate the bias in 500 simulation times
mean_bias1_1 <- mean(bias1_1, na.rm=TRUE) 
mean_bias1_2 <- mean(bias1_2, na.rm=TRUE) 
mean_bias1_3 <- mean(bias1_3, na.rm=TRUE) 

mean_bias2_1 <- mean(bias2_1, na.rm=TRUE)
mean_bias2_2 <- mean(bias2_2, na.rm=TRUE)
mean_bias2_3 <- mean(bias2_3, na.rm=TRUE)

mean_bias3_1 <- mean(bias3_1, na.rm=TRUE)
mean_bias3_2 <- mean(bias3_2, na.rm=TRUE)
mean_bias3_3 <- mean(bias3_3, na.rm=TRUE)

# calculate the standard deviation for the ipw estimate across 500 simulation times
sd_estimates_1 <- sd(ipw_estimates_1)
sd_estimates_2 <- sd(ipw_estimates_2)
sd_estimates_3 <- sd(ipw_estimates_3)


## confidence interval coverage rates 
mean(ipw_estimates_1 - 1.96*sd_estimates_1 <= patt_g1 & patt_g1 <= ipw_estimates_1 + 1.96*sd_estimates_1)
mean(ipw_estimates_2 - 1.96*sd_estimates_2 <= patt_g1 & patt_g1 <= ipw_estimates_2 + 1.96*sd_estimates_2)
mean(ipw_estimates_3 - 1.96*sd_estimates_3 <= patt_g1 & patt_g1 <= ipw_estimates_3 + 1.96*sd_estimates_3)

mean(ipw_estimates_1 - 1.96*sd_estimates_1 <= satt_g1 & satt_g1 <= ipw_estimates_1 + 1.96*sd_estimates_1)
mean(ipw_estimates_2 - 1.96*sd_estimates_2 <= satt_g1 & satt_g1 <= ipw_estimates_2 + 1.96*sd_estimates_2)
mean(ipw_estimates_3 - 1.96*sd_estimates_3 <= satt_g1 & satt_g1 <= ipw_estimates_3 + 1.96*sd_estimates_3)

mean(ipw_estimates_1 - 1.96*sd_estimates_1 <= patt_trt & patt_trt <= ipw_estimates_1 + 1.96*sd_estimates_1)
mean(ipw_estimates_2 - 1.96*sd_estimates_2 <= patt_trt & patt_trt <= ipw_estimates_2 + 1.96*sd_estimates_2)
mean(ipw_estimates_3 - 1.96*sd_estimates_3 <= patt_trt & patt_trt <= ipw_estimates_3 + 1.96*sd_estimates_3)

# calculate the average bootstrap standard deviation across 500 simulation times
sd_boot_1 <- mean(sd_did_1)
sd_boot_2 <- mean(sd_did_2)
sd_boot_3 <- mean(sd_did_3)