### Calculate the peak age for each call type for 2nd Gen FLUT study subjects (by ID!) and include all predictors ####
### XXX, 2024 ########################################################################################################

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(tidyr)
library(dplyr)

set.seed(23)

# what is the age range to be considered?
MAX_AGE <- 180
MIN_AGE <- 1

PROP_data <- readRDS('data/PROP_data.rds')
B_prop <- readRDS('models/B_prop.rds')

################################################################################
################################################################################

#### FUNCTIONS ####
find_peak_age_CI <- function(data, n_bootstrap = 1000, confidence_level = 0.95) {
  calculate_peak_age <- function(data) {
    tryCatch(
      {
        # Create interpolation functions for median estimates using spline interpolation
        interp_data <- with(data, splinefun(REC_AGE_D, Normalized_prop))
        
        # Create a sequence of REC_AGE_D values above the minimum age
        rec_age_seq <- seq(min(data$REC_AGE_D),  max(data$REC_AGE_D))
        
        # Calculate Call_prop values
        call_prop <- interp_data(rec_age_seq)
        
        peak_age_index <- which.max(call_prop)
        
        
        # Determine the point estimate
        peak_age <- rec_age_seq[peak_age_index]
        
        return(peak_age)
      },
      error = function(e) {
        cat("Error occurred:", conditionMessage(e), "\n")
        return(NA)
      }
    )
  }
  
  # Bootstrap resampling
  peak_ages <- replicate(n_bootstrap, {
    indices <- sample(nrow(data), replace = TRUE)
    calculate_peak_age(data[indices, ])
  })
  
  # Calculate point estimate
  point_estimate <- median(peak_ages, na.rm = TRUE)
  
  # Calculate confidence intervals
  ci_lower <- quantile(peak_ages, (1 - confidence_level) / 2, na.rm = TRUE)
  ci_upper <- quantile(peak_ages, 1 - (1 - confidence_level) / 2, na.rm = TRUE)
  
  return(list(peak_age = point_estimate, ci_lower = ci_lower, ci_upper = ci_upper))
}

PEAK_AGE <- function(row){
  print(row)
  rec_age_c <- seq(MIN_AGE, MAX_AGE, by=1)
  age_z_vals <- (rec_age_c - mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D)
  
  # predict call proportion for every day 
  PROP_pred <- B_prop %>% 
    epred_draws(newdata = expand_grid(TREATMENT = factor(row$TREATMENT),
                                      SEX = factor(row$SEX),
                                      AGE_z = age_z_vals,
                                      WEIGHT_z = row$Mean_WEIGHT_z,
                                      COMP_NORM_z = row$Mean_COMP_z,
                                      GS_z = row$Mean_GS_z,
                                      RAIN_z = row$Mean_RAIN_z,
                                      MaternalRank = row$MaternalRank,
                                      LITTER_CODE = row$LITTER_CODE,
                                      ID = row$ID,
                                      Total_calls=1),
                re_formula = NULL, allow_new_levels = F, robust = T)
  
  PROP_pred$Call_prop <- PROP_pred$.epred
  PROP_pred$Call_type <- as.factor(PROP_pred$.category)
  PROP_pred$TREATMENT <- as.factor(PROP_pred$TREATMENT)
  PROP_pred$SEX <- as.factor(PROP_pred$SEX)
  PROP_pred$REC_AGE_D <- PROP_pred$AGE_z*sd(PROP_data$REC_AGE_D) + mean(PROP_data$REC_AGE_D)
  PROP_pred$WEIGHT <- row$Mean_WEIGHT_DIFF
  PROP_pred$COMP <- row$Mean_COMP_NORM
  PROP_pred$GS <- row$Mean_GROUPSIZE
  PROP_pred$RAIN <- row$Mean_RAIN
  PROP_pred$MAT_STAT <- row$MaternalRank
  
  #normalise proportions:
  PROP_pred_sum <- PROP_pred %>%
    group_by(TREATMENT, REC_AGE_D, SEX, WEIGHT, COMP, GS, RAIN, MAT_STAT, ID, .draw) %>%
    mutate(Sum_prop = sum(Call_prop))
  
  rm(PROP_pred)
  
  PROP_pred_normalized <- PROP_pred_sum %>%
    mutate(Normalized_prop = Call_prop / Sum_prop)
  
  rm(PROP_pred_sum)
  
  # Calculate median estimates and credible intervals
  median_estimates <- PROP_pred_normalized %>%
    group_by(TREATMENT, REC_AGE_D, SEX, WEIGHT, COMP, GS, MAT_STAT, RAIN, ID, Call_type) %>%
    median_qi(Normalized_prop, .width =  .95)
  
  rm(PROP_pred_normalized)
  
  sumDIG_data <- median_estimates %>%
    filter(Call_type == "SumDIG")
  sumREP_data <- median_estimates %>%
    filter(Call_type == "SumBEG")
  sumCC_data <- median_estimates %>%
    filter(Call_type == "SumCC")
  
  treatment <- unique(median_estimates$TREATMENT)
  sex <- unique(median_estimates$SEX)
  weight <- unique(median_estimates$WEIGHT)
  comp <- unique(median_estimates$COMP)
  gs <- unique(median_estimates$GS)
  rain <- unique(median_estimates$RAIN)
  id <- unique(median_estimates$ID)
  mat_stat <- unique(median_estimates$MAT_STAT)
  
  rm(median_estimates)
  
  peak_age <- data.frame(TREATMENT = character(),
                         SEX = character(),
                         WEIGHT = numeric(), 
                         COMP = numeric(),
                         GS = numeric(),
                         RAIN = numeric(),
                         MAT_STAT = character(),
                         REP = numeric(),
                         DIG = numeric(),
                         CC = numeric(),
                         ID = numeric(),
                         stringsAsFactors = FALSE)
  
  sumDIG_data_subset <- sumDIG_data %>% filter(TREATMENT == treatment, SEX == sex, WEIGHT == weight, COMP == comp, GS == gs, RAIN == rain, MAT_STAT == mat_stat, ID==id)
  sumREP_data_subset <- sumREP_data %>% filter(TREATMENT == treatment, SEX == sex, WEIGHT == weight, COMP == comp, GS == gs, RAIN == rain, MAT_STAT == mat_stat, ID==id)
  sumCC_data_subset <- sumCC_data %>% filter(TREATMENT == treatment, SEX == sex, WEIGHT == weight, COMP == comp, GS == gs, RAIN == rain, MAT_STAT == mat_stat, ID==id)
  
  peak_age_REP <- find_peak_age_CI(sumREP_data_subset) 
  peak_age_DIG <- find_peak_age_CI(sumDIG_data_subset)
  peak_age_CC <- find_peak_age_CI(sumCC_data_subset)
  
  peak_age <- rbind(peak_age, data.frame(TREATMENT = treatment,
                                         SEX = sex,
                                         WEIGHT = weight,
                                         COMP = comp,
                                         GS = gs,
                                         RAIN = rain,
                                         MAT_STAT = mat_stat,
                                         REP = peak_age_REP,
                                         DIG = peak_age_DIG,
                                         CC = peak_age_CC,
                                         ID = id))
  
  print(peak_age)
  
  return(peak_age)
}

#### Calculate DIG peak age for all combinations of predictors ####
mean_weight <- mean(PROP_data$WEIGHT_DIFF_PER, na.rm=T)
mean_gs <- mean(PROP_data$GROUPSIZE, na.rm=T)
mean_comp <- mean(PROP_data$COMP_NORM, na.rm=T)
mean_rain <- mean(PROP_data$MonthlyRainfall, na.rm=T)

# create a new dataset with unique individuals and their mean values
FLUT_ID_data <- PROP_data %>%
  group_by(ID) %>%
  summarise(Mean_WEIGHT_DIFF = mean_weight,
            Mean_COMP_NORM = mean_comp,
            Mean_GROUPSIZE = mean_gs,
            Mean_RAIN = mean_rain,
            Mean_WEIGHT_z = 0,
            Mean_COMP_z = 0,
            Mean_GS_z = 0,
            Mean_RAIN_z = 0,
            MaternalRank  = first(MaternalRank),
            LITTER_CODE = first(LITTER_CODE),
            ID = first(ID),
            TREATMENT = first(TREATMENT),
            SEX = first(SEX))

# Initialize an empty data frame to store results
full_data <- data.frame(
  TREATMENT = character(),
  SEX = character(),
  WEIGHT = numeric(),
  COMP = numeric(),
  GS = numeric(),
  RAIN = numeric(),
  MAT_STAT = character(),
  DIG_REP = numeric(),
  CC_DIG = numeric(),
  CC_REP = numeric(),
  ID = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all possible combinations of weight, comp, and gs
for (i in 1:nrow(FLUT_ID_data)) {
  row = FLUT_ID_data[i,]
  peak_age <- PEAK_AGE(row = row)
  
  full_data <- rbind(full_data, peak_age)
  full_data$WEIGHT <- as.numeric(full_data$WEIGHT)
  full_data$COMP <- as.numeric(full_data$COMP)
  full_data$GS <- as.numeric(full_data$GS)
  full_data$RAIN <- as.numeric(full_data$RAIN)
  full_data$ID <- as.numeric(full_data$ID)
}

full_data$WEIGHT <- as.numeric(full_data$WEIGHT)
full_data$COMP <- as.numeric(full_data$COMP)
full_data$GS <- as.numeric(full_data$GS)
full_data$RAIN <- as.numeric(full_data$RAIN)
full_data$ID <- as.numeric(full_data$ID)

################################################################################
### RESULTS in right format ####
################################################################################
peak_data <- full_data

# convert to long format
data_long <- peak_data %>%
  select(ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MAT_STAT, REP.peak_age, DIG.peak_age, CC.peak_age) %>%
  pivot_longer(cols = c(REP.peak_age, DIG.peak_age, CC.peak_age),
               names_to = "Variable", values_to = "Peak_age") %>%
  mutate(Variable = gsub("\\.peak_age", "", Variable)) %>%
  left_join(
    peak_data %>%
      select(ID, REP.ci_lower, DIG.ci_lower, CC.ci_lower) %>%
      pivot_longer(cols = c(REP.ci_lower, DIG.ci_lower, CC.ci_lower),
                   names_to = "Variable", values_to = "ci_lower") %>%
      mutate(Variable = gsub("\\.ci_lower", "", Variable))
  ) %>%
  left_join(
    peak_data %>%
      select(ID, REP.ci_upper, DIG.ci_upper, CC.ci_upper) %>%
      pivot_longer(cols = c(REP.ci_upper, DIG.ci_upper, CC.ci_upper),
                   names_to = "Variable", values_to = "ci_upper") %>%
      mutate(Variable = gsub("\\.ci_upper", "", Variable))
  )

data_long$TREATMENT <- as.factor(data_long$TREATMENT)
data_long$TREATMENT <- factor(data_long$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
data_long$SEX <- as.factor(data_long$SEX)
data_long$MAT_STAT <-as.factor(data_long$MAT_STAT)

data_long <- na.omit(data_long)

