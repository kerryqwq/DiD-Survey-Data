did_functions = function(data){
  # estimate propensity scores by using multinomial logistic regression
  fit <- multinom(treated_label ~ x1, data = data, trace=FALSE) 
  
  dat$psvalue <- predict(fit, type = "probs") # find the propensity scores
  dat$psvalue <- ifelse(dat$psvalue < 0.01, 0.01, ifelse(dat$psvalue > 0.99, 0.99, dat$psvalue))
  
  data <- data %>%
    mutate(ps_weight = case_when( # calculate the propensity score weights by groups
      treated_label == 1 ~ psvalue[,"1"] / psvalue[,"1"], 
      treated_label == 2 ~ psvalue[,"1"] / psvalue[,"2"],
      treated_label == 3 ~ psvalue[,"1"] / psvalue[,"3"],
      treated_label == 4 ~ psvalue[,"1"] / psvalue[,"4"] 
    ),
    psweight_true = case_when(
      treated_label == 1 ~ 1,
      treated_label == 2 ~ plogis(gamma0[1] + gamma0[2]*x1)/plogis(gamma1[1] + gamma1[2]*x1),
      treated_label == 3 ~ plogis(gamma0[1] + gamma0[2]*x1)/(1-plogis(gamma0[1] + gamma0[2]*x1)),
      treated_label == 4 ~ plogis(gamma0[1] + gamma0[2]*x1)/(1-plogis(gamma1[1] + gamma1[2]*x1))
    ),
    samp_prob = case_when( 
      tp == 0 ~ plogis(eta0[1] + eta0[2]*x1),
      tp == 1 ~ plogis(eta1[1] + eta1[2]*x1)
    ),
    samp_weight = 1 / samp_prob # calculate the sampling weights
    ) # calculate the final weights
  pihat = mean(dat$treated_label == 1)
  # calculate the target estimate
  
  #dat$samp_weight = 2*N*dat$samp_weight/sum(dat$samp_weight)
  data$final_weight = data$ps_weight * data$samp_weight
  
  # DiD IPW estimator with the final weights
  did_1 <- with(data, sum(y[treated_label == 2]*final_weight[treated_label == 2])/ (sum(samp_weight[treated_label == 1])) -
                  sum(y[treated_label == 1]*final_weight[treated_label == 1])/(sum(samp_weight[treated_label == 1]))   -
                  sum(y[treated_label == 4]*final_weight[treated_label == 4])/(sum(samp_weight[treated_label == 1]))  +
                  sum(y[treated_label == 3]*final_weight[treated_label == 3])/(sum(samp_weight[treated_label == 1])) )
  
  # DiD IPW estimator with the propensity score weights
  did_2 <- with(data, sum(y[treated_label == 2]*ps_weight[treated_label == 2])/(sum(treated_label == 1)) -
                  sum(y[treated_label == 1]*ps_weight[treated_label == 1])/(sum(treated_label == 1))   -
                  sum(y[treated_label == 4]*ps_weight[treated_label == 4])/(sum(treated_label == 1))  +
                  sum(y[treated_label == 3]*ps_weight[treated_label == 3])/(sum(treated_label == 1)) )
  
  # DiD IPW estimator with the sampling weights
  did_3 <- with(data, sum(y[treated_label == 2]*samp_weight[treated_label == 2])/ (sum(samp_weight[treated_label == 2])) -
                  sum(y[treated_label == 1]*samp_weight[treated_label == 1])/(sum(samp_weight[treated_label == 1]))   -
                  sum(y[treated_label == 4]*samp_weight[treated_label == 4])/(sum(samp_weight[treated_label == 4]))  +
                  sum(y[treated_label == 3]*samp_weight[treated_label == 3])/(sum(samp_weight[treated_label == 3])))
  
  return(c(did_1, did_2, did_3))
  
}



did_functions_real = function(data){
  # estimate propensity scores by using multinomial logistic regression
  fit <- multinom(comb_group_label ~ sex + age + bmi + race4, data = data, trace=FALSE) 
  
  data$psvalue <- predict(fit, type = "probs") # find the propensity scores
  data$psvalue <- ifelse(data$psvalue < 0.01, 0.01, ifelse(data$psvalue > 0.99, 0.99, data$psvalue))
  
  data <- data %>%
    mutate(ps_weight = case_when( # calculate the propensity score weights by groups
      comb_group_label == 1 ~ psvalue[,"1"] / psvalue[,"1"], 
      comb_group_label == 2 ~ psvalue[,"1"] / psvalue[,"2"],
      comb_group_label == 3 ~ psvalue[,"1"] / psvalue[,"3"],
      comb_group_label == 4 ~ psvalue[,"1"] / psvalue[,"4"] 
    )) # calculate the final weights
  
  # calculate the target estimate
  data$final_weight = data$ps_weight * data$weight
  
  # DiD IPW estimator with the final weights
  did_1 <- with(data, sum(soda_usage[comb_group_label == 2]*final_weight[comb_group_label == 2])/ (sum(weight[comb_group_label == 1])) -
                  sum(soda_usage[comb_group_label == 1]*final_weight[comb_group_label == 1])/(sum(weight[comb_group_label == 1]))   -
                  sum(soda_usage[comb_group_label == 4]*final_weight[comb_group_label == 4])/(sum(weight[comb_group_label == 1]))  +
                  sum(soda_usage[comb_group_label == 3]*final_weight[comb_group_label == 3])/(sum(weight[comb_group_label == 1])) )
  
  # DiD IPW estimator with the propensity score weights
  did_2 <- with(data, sum(soda_usage[comb_group_label == 2]*ps_weight[comb_group_label == 2])/(sum(comb_group_label == 1)) -
                  sum(soda_usage[comb_group_label == 1]*ps_weight[comb_group_label == 1])/(sum(comb_group_label == 1))   -
                  sum(soda_usage[comb_group_label == 4]*ps_weight[comb_group_label == 4])/(sum(comb_group_label == 1))  +
                  sum(soda_usage[comb_group_label == 3]*ps_weight[comb_group_label == 3])/(sum(comb_group_label == 1)) )
  
  # DiD IPW estimator with the sampling weights
  did_3 <- with(data, sum(soda_usage[comb_group_label == 2]*weight[comb_group_label == 2])/ (sum(weight[comb_group_label == 2])) -
                  sum(soda_usage[comb_group_label == 1]*weight[comb_group_label == 1])/(sum(weight[comb_group_label == 1]))   -
                  sum(soda_usage[comb_group_label == 4]*weight[comb_group_label == 4])/(sum(weight[comb_group_label == 4]))  +
                  sum(soda_usage[comb_group_label == 3]*weight[comb_group_label == 3])/(sum(weight[comb_group_label == 3])))
  
  return(c(did_1, did_2, did_3))
  
}


