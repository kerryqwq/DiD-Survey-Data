library(tidyverse)
library(knitr)
library(nnet)
library(MASS)

source("Code/bootfunctione.R")

N <- 500 # population number
max.time <- 1 # maximum time period: two period
trt.time <- 1 # the second period is the post-trt period
nsims <- 500 # simulation times
nboot <- 100
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

smd.x1.unweighted = smd.x2.unweighted = smd.x3.unweighted= smd.x4.unweighted = rep(NA, nsims)
smd.x1.ps_weighted = smd.x2.ps_weighted = smd.x3.ps_weighted = smd.x4.ps_weighted = rep(NA, nsims)
smd.x1.sample_weighted = smd.x2.sample_weighted = smd.x3.sample_weighted = smd.x4.sample_weighted = rep(NA, nsims)
smd.x1.final_weighted = smd.x2.final_weighted = smd.x3.final_weighted = smd.x4.final_weighted = rep(NA, nsims)

gamma2 = c(1, -0.5, -0.5, 0.0, -1.0)
gamma4 = c(-1, 0.5, 1.00, 0.0, -0.5)
gamma3 = c(0, 1, 0.2, 0.5, 0.5)
eta = c(0.5, 0.50, -1.00, 1.00, 0.0)


mu <- c(0, 0)  
Sigma <- matrix(c(2, 0.3, 0.3, 1), ncol = 2)  # Covariance matrix

for (i in 1:nsims) {
  
  print(i)
  set.seed(i)
  G = S = D = rep(NA, N)
  
  mvndat = mvrnorm(N, mu = mu, Sigma = Sigma)
  x1 = mvndat[,1]; x2 = mvndat[,2]
  x3 = rbinom(N, 1, 0.5)
  x4 = runif(N, 0, 1)
  S = rbinom(N, 1, plogis(eta[1] + eta[2]*x1 + eta[3]*x2 + eta[4]*x3 +eta[5]*x4))
 
  
  exp_logits <- cbind(1, exp(gamma2[1] + gamma2[2]*x1 + gamma2[3]*x2 + gamma2[4]*x3 + gamma2[5]*x4), 
                      exp(gamma3[1] + gamma3[2]*x1 + gamma3[3]*x2 + gamma3[4]*x3 + gamma3[5]*x4), 
                      exp(gamma4[1] + gamma4[2]*x1 + gamma4[3]*x2 + gamma4[4]*x3 + gamma4[5]*x4)) 
  probs = exp_logits[S==1,] / rowSums(exp_logits[S==1,])
  G[S==1] = apply(probs, 1, function(p) which(rmultinom(1,1,p)==1))
  probs2 = exp_logits[S==0,] / rowSums(exp_logits[S==0,])
  G[S==0] = apply(probs2, 1, function(p) which(rmultinom(1, 1, p) == 1))
  D = ifelse(G <= 2, 1, 2)
  
  # generate the potential outcomes 
  y0_t0 = 0.5*x1 + 0.05*x2 + 0.2*x3 + 0.15*x4 + D + rnorm(N) 
  y0_t1 = y0_t0 + 0.5 + rnorm(N)
  y1_t0 = y0_t0 + 1 + 0.5*x1 - 0.30*x2 + 0.5*x3 
  y1_t1 = y0_t1 + 1 + 0.5*x1 - 0.30*x2 + 0.5*x3
  
  observed.Y = ifelse(G == 1, y0_t0, ifelse(G == 2, y1_t1, ifelse(G ==3, y0_t0, y0_t1)))
   
  orig_dat = data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, S = S, G = G, D = D,
                        y0_t0 = y0_t0, y0_t1 = y0_t1, y1_t0 = y1_t0, y1_t1 = y1_t1, 
                        observed.Y = observed.Y)
  sample_dat = orig_dat[orig_dat$S == 1,]
  
  # target estimand
  patt_trt[i] <- mean(orig_dat$y1_t1[orig_dat$G == 2]) - mean(orig_dat$y0_t0[orig_dat$G == 1]) - 
    (mean(orig_dat$y0_t1[orig_dat$G == 4]) - mean(orig_dat$y0_t0[orig_dat$G == 3]))
  
  patt_g1[i] <- mean(orig_dat$y1_t1[orig_dat$G ==1]) - mean(orig_dat$y0_t1[orig_dat$G ==1])
  
  satt_g1[i] <- mean(orig_dat$y1_t1[orig_dat$G ==1 & orig_dat$S ==1]) - mean(orig_dat$y0_t1[orig_dat$G ==1 & orig_dat$S ==1])
  
  # estimate propensity scores by using multinomial logistic regression
  fit <- multinom(G ~ x1 + x2 + x3 + x4, data = sample_dat, trace=FALSE) 
  
  sample_dat$psvalue <- predict(fit, type = "probs") # find the propensity scores
  sample_dat$psvalue <- ifelse(sample_dat$psvalue < 0.01, 0.01, ifelse(sample_dat$psvalue > 0.99, 0.99, sample_dat$psvalue))
  
  sample_dat <- sample_dat %>%
    mutate(ps_weight = case_when( # calculate the propensity score weights by groups
      G == 1 ~ psvalue[,"1"] / psvalue[,"1"], 
      G == 2 ~ psvalue[,"1"] / psvalue[,"2"],
      G == 3 ~ psvalue[,"1"] / psvalue[,"3"],
      G == 4 ~ psvalue[,"1"] / psvalue[,"4"] 
    ),
    samp_prob = plogis(eta[1] + eta[2]*x1 + eta[3]*x2 + eta[4]*x3 +eta[5]*x4),
    samp_weight = 1 / samp_prob # calculate the sampling weights
    ) 
  
  sample_dat$final_weight = sample_dat$ps_weight * sample_dat$samp_weight   # calculate the final weights
  
  # our proposed estimator
  did_1 <- with(sample_dat, sum(observed.Y[G == 2]*final_weight[G == 2])/(sum(samp_weight[G == 1])) -
                  sum(observed.Y[G == 1]*final_weight[G == 1])/(sum(samp_weight[G == 1])) -
                  sum(observed.Y[G == 4]*final_weight[G == 4])/(sum(samp_weight[G == 1])) +
                  sum(observed.Y[G == 3]*final_weight[G == 3])/(sum(samp_weight[G == 1])) )
  
  # tau hat pw
  did_2 <- with(sample_dat, sum(observed.Y[G == 2]*ps_weight[G == 2])/(sum(ps_weight[G == 1])) -
                  sum(observed.Y[G == 1]*ps_weight[G == 1])/(sum(ps_weight[G == 1]))   -
                  sum(observed.Y[G == 4]*ps_weight[G == 4])/(sum(ps_weight[G == 1]))  +
                  sum(observed.Y[G == 3]*ps_weight[G == 3])/(sum(ps_weight[G == 1])) )
  
  # tau hat sw
  did_3 <- with(sample_dat, sum(observed.Y[G == 2]*samp_weight[G == 2])/ (sum(samp_weight[G == 2])) -
                  sum(observed.Y[G == 1]*samp_weight[G == 1])/(sum(samp_weight[G == 1]))   -
                  sum(observed.Y[G == 4]*samp_weight[G == 4])/(sum(samp_weight[G == 4]))  +
                  sum(observed.Y[G == 3]*samp_weight[G == 3])/(sum(samp_weight[G == 3])))
  

  # Weighted histogram
  smd.x1.unweighted[i] = (mean(orig_dat$x1[orig_dat$G==1])-mean(sample_dat$x1))/sqrt(var(orig_dat$x1[orig_dat$G == 1]))
  smd.x1.ps_weighted[i] = (mean(orig_dat$x1[orig_dat$G==1])-sum(sample_dat$x1*sample_dat$ps_weight)/sum(sample_dat$ps_weight))/sqrt(var(orig_dat$x1[orig_dat$G == 1]))
  smd.x1.sample_weighted[i] = (mean(orig_dat$x1[orig_dat$G==1])-sum(sample_dat$x1*sample_dat$samp_weight)/sum(sample_dat$samp_weight))/sqrt(var(orig_dat$x1[orig_dat$G == 1]))
  smd.x1.final_weighted[i] = (mean(orig_dat$x1[orig_dat$G==1])-sum(sample_dat$x1*sample_dat$final_weight)/sum(sample_dat$final_weight))/sqrt(var(orig_dat$x1[orig_dat$G == 1]))
  
  smd.x2.unweighted[i] = (mean(orig_dat$x2[orig_dat$G==1])-mean(sample_dat$x2))/sqrt(var(orig_dat$x2[orig_dat$G == 1]))
  smd.x2.ps_weighted[i] = (mean(orig_dat$x2[orig_dat$G==1])-sum(sample_dat$x2*sample_dat$ps_weight)/sum(sample_dat$ps_weight))/sqrt(var(orig_dat$x2[orig_dat$G == 1]))
  smd.x2.sample_weighted[i] = (mean(orig_dat$x2[orig_dat$G==1])-sum(sample_dat$x2*sample_dat$samp_weight)/sum(sample_dat$samp_weight))/sqrt(var(orig_dat$x2[orig_dat$G == 1]))
  smd.x2.final_weighted[i] = (mean(orig_dat$x2[orig_dat$G==1])-sum(sample_dat$x2*sample_dat$final_weight)/sum(sample_dat$final_weight))/sqrt(var(orig_dat$x2[orig_dat$G == 1]))
  
  smd.x3.unweighted[i] = (mean(orig_dat$x3[orig_dat$G==1])-mean(sample_dat$x3))/sqrt(var(orig_dat$x3[orig_dat$G == 1]))
  smd.x3.ps_weighted[i] = (mean(orig_dat$x3[orig_dat$G==1])-sum(sample_dat$x3*sample_dat$ps_weight)/sum(sample_dat$ps_weight))/sqrt(var(orig_dat$x3[orig_dat$G == 1]))
  smd.x3.sample_weighted[i] = (mean(orig_dat$x3[orig_dat$G==1])-sum(sample_dat$x3*sample_dat$samp_weight)/sum(sample_dat$samp_weight))/sqrt(var(orig_dat$x3[orig_dat$G == 1]))
  smd.x3.final_weighted[i] = (mean(orig_dat$x3[orig_dat$G==1])-sum(sample_dat$x3*sample_dat$final_weight)/sum(sample_dat$final_weight))/sqrt(var(orig_dat$x3[orig_dat$G == 1]))
  
  smd.x4.unweighted[i] = (mean(orig_dat$x4[orig_dat$G==1])-mean(sample_dat$x4))/sqrt(var(orig_dat$x4[orig_dat$G == 1]))
  smd.x4.ps_weighted[i] = (mean(orig_dat$x4[orig_dat$G==1])-sum(sample_dat$x4*sample_dat$ps_weight)/sum(sample_dat$ps_weight))/sqrt(var(orig_dat$x4[orig_dat$G == 1]))
  smd.x4.sample_weighted[i] = (mean(orig_dat$x4[orig_dat$G==1])-sum(sample_dat$x4*sample_dat$samp_weight)/sum(sample_dat$samp_weight))/sqrt(var(orig_dat$x4[orig_dat$G == 1]))
  smd.x4.final_weighted[i] = (mean(orig_dat$x4[orig_dat$G==1])-sum(sample_dat$x4*sample_dat$final_weight)/sum(sample_dat$final_weight))/sqrt(var(orig_dat$x4[orig_dat$G == 1]))
  

  bootstrap_results = matrix(NA, nboot, 3)
  for(r in 1:nboot){
    boot_dat <- sample_dat[sample(1:nrow(sample_dat), nrow(sample_dat), replace = TRUE),]
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


## Table 1 
tab = matrix(NA, nrow = 4, ncol = 4)
tab[,1] = c(paste0(round(mean(smd.x1.unweighted),2), " (", round(sd(smd.x1.unweighted), 2) ,")"),
            paste0(round(mean(smd.x1.ps_weighted),2), " (", round(sd(smd.x1.ps_weighted), 2) ,")"),
            paste0(round(mean(smd.x1.sample_weighted),2), " (", round(sd(smd.x1.sample_weighted), 2) ,")"),
            paste0(round(mean(smd.x1.final_weighted),2), " (", round(sd(smd.x1.final_weighted), 2) ,")"))

tab[,2] = c(paste0(round(mean(smd.x2.unweighted),2), " (", round(sd(smd.x2.unweighted), 2) ,")"),
            paste0(round(mean(smd.x2.ps_weighted),2), " (", round(sd(smd.x2.ps_weighted), 2) ,")"),
            paste0(round(mean(smd.x2.sample_weighted),2), " (", round(sd(smd.x2.sample_weighted), 2) ,")"),
            paste0(round(mean(smd.x2.final_weighted),2), " (", round(sd(smd.x2.final_weighted), 2) ,")"))

tab[,3] = c(paste0(round(mean(smd.x3.unweighted),2), " (", round(sd(smd.x3.unweighted), 2) ,")"),
            paste0(round(mean(smd.x3.ps_weighted),2), " (", round(sd(smd.x3.ps_weighted), 2) ,")"),
            paste0(round(mean(smd.x3.sample_weighted),2), " (", round(sd(smd.x3.sample_weighted), 2) ,")"),
            paste0(round(mean(smd.x3.final_weighted),2), " (", round(sd(smd.x3.final_weighted), 2) ,")"))

tab[,4] = c(paste0(round(mean(smd.x4.unweighted),2), " (", round(sd(smd.x4.unweighted), 2) ,")"),
            paste0(round(mean(smd.x4.ps_weighted),2), " (", round(sd(smd.x4.ps_weighted), 2) ,")"),
            paste0(round(mean(smd.x4.sample_weighted),2), " (", round(sd(smd.x4.sample_weighted), 2) ,")"),
            paste0(round(mean(smd.x4.final_weighted),2), " (", round(sd(smd.x4.final_weighted), 2) ,")"))

xtable(tab)



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


## Table 2 (N=500): target: patt_g1
tab2 = matrix(NA, nrow = 3, ncol = 5)
tab2[1,] = c(round(mean(bias1_1, na.rm=TRUE),2), round(sqrt(mean(bias1_1^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_1), 
            mean(ipw_estimates_1 - 1.96*sd(ipw_estimates_1)<= patt_g1 & patt_g1 <= ipw_estimates_1 + 1.96*sd(ipw_estimates_1)),
            mean(ipw_estimates_1 - 1.96*sd_did_1<= patt_g1 & patt_g1 <= ipw_estimates_1 + 1.96*sd_did_1))
tab2[2,] = c(round(mean(bias2_1, na.rm=TRUE),2), round(sqrt(mean(bias2_1^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_2), 
             mean(ipw_estimates_2 - 1.96*sd(ipw_estimates_2)<= patt_g1 & patt_g1 <= ipw_estimates_2 + 1.96*sd(ipw_estimates_2)),
             mean(ipw_estimates_2 - 1.96*sd_did_2<= patt_g1 & patt_g1 <= ipw_estimates_2 + 1.96*sd_did_2))
tab2[3,] = c(round(mean(bias3_1, na.rm=TRUE),2), round(sqrt(mean(bias3_1^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_3), 
             mean(ipw_estimates_3 - 1.96*sd(ipw_estimates_3)<= patt_g1 & patt_g1 <= ipw_estimates_3 + 1.96*sd(ipw_estimates_3)),
             mean(ipw_estimates_3 - 1.96*sd_did_3<= patt_g1 & patt_g1 <= ipw_estimates_3 + 1.96*sd_did_3))
xtable(tab2)


## Table A1: 
tab3 = matrix(NA, nrow = 3, ncol = 5)
tab3[1,] = c(round(mean(bias1_2, na.rm=TRUE),2), round(sqrt(mean(bias1_2^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_1), 
             mean(ipw_estimates_1 - 1.96*sd(ipw_estimates_1)<= satt_g1 & satt_g1 <= ipw_estimates_1 + 1.96*sd(ipw_estimates_1)),
             mean(ipw_estimates_1 - 1.96*sd_did_1<= satt_g1 & satt_g1 <= ipw_estimates_1 + 1.96*sd_did_1))
tab3[2,] = c(round(mean(bias2_2, na.rm=TRUE),2), round(sqrt(mean(bias2_2^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_2), 
             mean(ipw_estimates_2 - 1.96*sd(ipw_estimates_2)<= satt_g1 & satt_g1 <= ipw_estimates_2 + 1.96*sd(ipw_estimates_2)),
             mean(ipw_estimates_2 - 1.96*sd_did_2<= satt_g1 & satt_g1 <= ipw_estimates_2 + 1.96*sd_did_2))
tab3[3,] = c(round(mean(bias3_2, na.rm=TRUE),2), round(sqrt(mean(bias3_2^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_3), 
             mean(ipw_estimates_3 - 1.96*sd(ipw_estimates_3)<= satt_g1 & satt_g1 <= ipw_estimates_3 + 1.96*sd(ipw_estimates_3)),
             mean(ipw_estimates_3 - 1.96*sd_did_3<= satt_g1 & satt_g1 <= ipw_estimates_3 + 1.96*sd_did_3))
xtable(tab3)

## Table A2: 
tab4 = matrix(NA, nrow = 3, ncol = 5)
tab4[1,] = c(round(mean(bias1_3, na.rm=TRUE),2), round(sqrt(mean(bias1_3^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_1), 
             mean(ipw_estimates_1 - 1.96*sd(ipw_estimates_1)<= patt_trt & patt_trt <= ipw_estimates_1 + 1.96*sd(ipw_estimates_1)),
             mean(ipw_estimates_1 - 1.96*sd_did_1<= patt_trt & patt_trt <= ipw_estimates_1 + 1.96*sd_did_1))
tab4[2,] = c(round(mean(bias2_3, na.rm=TRUE),2), round(sqrt(mean(bias2_3^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_2), 
             mean(ipw_estimates_2 - 1.96*sd(ipw_estimates_2)<= patt_trt & patt_trt <= ipw_estimates_2 + 1.96*sd(ipw_estimates_2)),
             mean(ipw_estimates_2 - 1.96*sd_did_2<= patt_trt & patt_trt <= ipw_estimates_2 + 1.96*sd_did_2))
tab4[3,] = c(round(mean(bias3_3, na.rm=TRUE),2), round(sqrt(mean(bias3_3^2, na.rm = TRUE)),2), mean(1.96*2*sd_did_3), 
             mean(ipw_estimates_3 - 1.96*sd(ipw_estimates_3)<= patt_trt & patt_trt <= ipw_estimates_3 + 1.96*sd(ipw_estimates_3)),
             mean(ipw_estimates_3 - 1.96*sd_did_3<= patt_trt & patt_trt <= ipw_estimates_3 + 1.96*sd_did_3))
xtable(tab4)
