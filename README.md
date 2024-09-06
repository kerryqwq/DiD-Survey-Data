# Difference-in-Differences Analysis with Survey Data

## Overview
Difference-in-differences (DiD) approach with the underlying counterfactual parallel trends assumption is one of the most widely used approaches for policy effect evaluation. However, traditional DiD methods may not be directly applicable in the absence of population-level panel data, particularly when the target population could substantially differ from the treated units observed during the post-intervention period. In this work, we focus on addressing the following two
challenges in using DiD methods with repeated cross-sectional (RCS) survey data: (1) the heterogeneous compositions of study samples across different time points, and (2) the availability of only a subset of the population. We formally introduce the policy-relevant target estimand and establish its identification conditions. We propose a new weighting approach that incorporates both estimated propensity scores and survey weights. We establish the theoretical properties of the proposed method and examine its finite sample performance through simulations. Finally, we apply our proposed method to a real-world data application, estimating the effect of a beverage tax on adolescent soda consumption in Philadelphia.

## Data

Our data application of the beverage tax on adolescent soda consumption is based on this publicly available data----Youth Risk Behavior Surveillance (YRBS) data. The YRBS data provide a school-district level biennial survey from the YRBSS, managed by the Centers for Disease Control and
Prevention (CDC). 
The data used in this study come from the following sources:

- **Original Data File:** `Data/SADCQ.csv`
- **Processed Data for Analysis:** `Data/realdata.csv`

The pre-processing procedures are demonstrated in the script `Code/did_preprocessing.R`, and the main analysis is conducted in `Code/real_data.R`. A data dictionary is also available for reference.

**Variables**

The key variables used in the analysis are as follows:

- **year:** Year when the survey was conducted. The year values are `2013`, `2015`, `2017`, and `2019`.
    
- **sex:** Binary indicator of the adolescent's sex.  
  - `1` = Female  
  - `2` = Male  

- **age:** Age of the adolescent.  
  - `1` = ≤ 12 years  
  - `2` = 13 years  
  - `3` = 14 years  
  - `4` = 15 years  
  - `5` = 16 years  
  - `6` = 17 years  
  - `7` = ≥ 18 years  

- **bmi:** Adolescent's Body Mass Index (BMI).

- **race4:** Race of the adolescent.  
  - `1` = White  
  - `2` = Black or African American  
  - `3` = Hispanic/Latino  
  - `4` = All Other Races  

- **soda_usage:** Average weekly soda consumption by the adolescent, categorized as follows:  
  - `0` = 0 sodas per week  
  - `2` = 1-3 sodas per week  
  - `5` = 4-6 sodas per week  
  - `7` = 7 sodas per week  
  - `14` = 14 sodas per week  
  - `21` = 21 sodas per week  
  - `28` = 28 or more sodas per week  

- **comb_group_label:** Group G indicator based on the treatment groups and time periods.  
  - `1` = Treated unit in the pre-treatment period  
  - `2` = Treated unit in the post-treatment period  
  - `3` = Control unit in the pre-treatment period  
  - `4` = Control unit in the post-treatment period  

- **weight:** Survey weights used in the analysis.
  
- **interaction:** Binary indicator representing the interaction term between the time period and the treatment group.  
  - `1` if the unit is from the treated group and in the post-treatment period  
  - `0` otherwise

## Code 
-`simulation.R`: In this simulation code, we demonstrate the finite sample performance of the proposed estimator, using the simulated RCS survey data. We consider total population of size N (500, 1000, 2000), divided equally between two time periods (t=0,1).  We compare our proposed estimator with other two IPW estimators. 
 
-`bootfunction.R`:In this bootstrap function code, we define the bootstrap function for both the simulation and the real data application.
 
-`real_data.R`: In this real data application code, we apply our proposed estimator, to the YRBS data to examine the effect of the beverage tax on soda consumption among high school students in Philadelphia.
