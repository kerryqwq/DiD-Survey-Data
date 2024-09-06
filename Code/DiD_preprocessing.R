# Load the packages
library(tidyverse)


# Load the data
did_penn <- read.csv("/Users/kerryqwq/Downloads/SADCQ.csv")
head(did_penn)

# Subset the data for specific district names, years, and variables
did_penn_sub_w_na <-did_penn %>%
  group_by(sitename) %>%
  filter(sitename %in% c('Philadelphia, PA (PH)', 'New York City, NY (NYC)', 'Orange County, FL (OL)', 'Palm Beach County, FL (PB)',
                         'Broward County, FL (FT)', 'San Diego, CA (SA)', 'Los Angeles, CA (LO)')) %>%
  filter(year %in% c(2013,2015,2017,2019)) %>%
  select(sitename, sitecode, year, sex, age, bmi, race4, q75, weight, stratum)

#rescale the outcomes, create post-treatment period labels and treatment group labels.
did_penn_sub_w_na <- did_penn_sub_w_na %>%
  rename(outcome=q75) %>%
  mutate(soda_usage=case_when(outcome==1 ~ 0,outcome==2 ~ 2, outcome==3 ~ 5,outcome==4 ~ 7, outcome==5 ~ 14, 
                              outcome==6 ~ 21, outcome==7 ~ 28)) %>%
  filter(!is.na(soda_usage)) %>%
  mutate(year_label = case_when(year %in% c(2017,2019) ~ 1, year %in% c(2013,2015) ~ 0),
         sitecode_label = case_when(sitecode %in% "PH" ~ 1, sitecode %in% c("FT","LO","NYC","OL","PB","SA") ~ 0),
         interaction = year_label * sitecode_label)

# Create a combined group label for all four groups based on treatment and period
did_penn_sub_w_na <- did_penn_sub_w_na %>% 
  mutate(comb_group_label = case_when(
    sitecode_label ==1 & year_label == 0 ~ 1, # label as group 1 if the treated unit is in pre-treatment period
    sitecode_label ==1 & year_label == 1 ~ 2, # label as group 2 if the treated unit is in post-treatment period
    sitecode_label ==0 & year_label == 0 ~ 3, # label as group 3 if the control unit is in pre-treatment period
    sitecode_label ==0 & year_label == 1 ~ 4  # label as group 4 if the control unit is in post-treatment period
  ))


# Factor the categorical columns and remove rows with missing values
did_penn_sub <- did_penn_sub_w_na %>%
  mutate_at(c("sex", "race4", "year_label", "sitecode_label", "interaction", "comb_group_label"), factor) %>%
  mutate(comb_group_label = factor(comb_group_label, levels = c("1", "2", "3", "4"))) %>%
  drop_na()

