did_functions = function(data){
  
  fit <- multinom(G ~ x1 + x2 + x3 + x4, data = data, trace=FALSE) 
  
  data$psvalue <- predict(fit, type = "probs") # find the propensity scores
  data$psvalue <- ifelse(data$psvalue < 0.01, 0.01, ifelse(data$psvalue > 0.99, 0.99, data$psvalue))
  
  data <- data %>%
    mutate(ps_weight = case_when( # calculate the propensity score weights by groups
      G == 1 ~ psvalue[,"1"] / psvalue[,"1"], 
      G == 2 ~ psvalue[,"1"] / psvalue[,"2"],
      G == 3 ~ psvalue[,"1"] / psvalue[,"3"],
      G == 4 ~ psvalue[,"1"] / psvalue[,"4"] 
    ),
    samp_prob = plogis(eta[1] + eta[2]*x1 + eta[3]*x2 + eta[4]*x3 +eta[5]*x4),
    samp_weight = 1 / samp_prob # calculate the sampling weights
    )
   
  data$final_weight = data$ps_weight * data$samp_weight 
  
  # DiD IPW estimator with the final weights
  did_1 <- with(data, sum(observed.Y[G == 2]*final_weight[G == 2])/(sum(samp_weight[G == 1])) -
                  sum(observed.Y[G == 1]*final_weight[G == 1])/(sum(samp_weight[G == 1])) -
                  sum(observed.Y[G == 4]*final_weight[G == 4])/(sum(samp_weight[G == 1])) +
                  sum(observed.Y[G == 3]*final_weight[G == 3])/(sum(samp_weight[G == 1])) )
  
  # tau hat pw
  did_2 <- with(data, sum(observed.Y[G == 2]*ps_weight[G == 2])/(sum(ps_weight[G == 1])) -
                  sum(observed.Y[G == 1]*ps_weight[G == 1])/(sum(ps_weight[G == 1]))   -
                  sum(observed.Y[G == 4]*ps_weight[G == 4])/(sum(ps_weight[G == 1]))  +
                  sum(observed.Y[G == 3]*ps_weight[G == 3])/(sum(ps_weight[G == 1])) )
  
  # tau hat sw
  did_3 <- with(data, sum(observed.Y[G == 2]*samp_weight[G == 2])/ (sum(samp_weight[G == 2])) -
                  sum(observed.Y[G == 1]*samp_weight[G == 1])/(sum(samp_weight[G == 1]))   -
                  sum(observed.Y[G == 4]*samp_weight[G == 4])/(sum(samp_weight[G == 4]))  +
                  sum(observed.Y[G == 3]*samp_weight[G == 3])/(sum(samp_weight[G == 3])))
  
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



