### Calculate transition ages for 2nd Gen FLUT study subjects (by ID!) and include all predictors ####
### XXX, 2024 ########################################################################################

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
find_transition_age_CI <- function(data_higher, data_lower, check_peak_dig = FALSE, n_bootstrap = 1000, confidence_level = 0.95) {
  calculate_transition_age <- function(data_higher, data_lower, check_peak_dig) {
    tryCatch(
      {
        # Create interpolation functions for median estimates using spline interpolation
        interp_higher <- with(data_higher, splinefun(REC_AGE_D, Normalized_prop))
        interp_lower <- with(data_lower, splinefun(REC_AGE_D, Normalized_prop))
        
        # Create a sequence of REC_AGE_D values above the minimum age
        rec_age_seq <- seq(min(min(data_higher$REC_AGE_D), min(data_lower$REC_AGE_D)), max(max(data_higher$REC_AGE_D), max(data_lower$REC_AGE_D)))

        # Calculate Call_prop values
        call_prop_higher <- interp_higher(rec_age_seq)
        call_prop_lower <- interp_lower(rec_age_seq)

        # Determine the peak age of DIG
        peak_dig_age <- rec_age_seq[which.max(call_prop_lower)]
        
        # Find the index where higher > lower and within confidence intervals
        # Only consider ages after the peak DIG age if check_peak_dig is TRUE
        if (check_peak_dig) {
          transition_age_index <- which(call_prop_higher > call_prop_lower &
                                          rec_age_seq > peak_dig_age & # Ensure age is after peak DIG
                                          data_higher$.lower <= call_prop_higher &
                                          data_higher$.upper >= call_prop_higher &
                                          data_lower$.lower <= call_prop_lower &
                                          data_lower$.upper >= call_prop_lower)[1]
        } else {
          transition_age_index <- which(call_prop_higher > call_prop_lower &
                                          data_higher$.lower <= call_prop_higher &
                                          data_higher$.upper >= call_prop_higher &
                                          data_lower$.lower <= call_prop_lower &
                                          data_lower$.upper >= call_prop_lower)[1]
        }
        
        # Determine the point estimate
        transition_age <- rec_age_seq[transition_age_index]

        return(transition_age)
      },
      error = function(e) {
        cat("Error occurred:", conditionMessage(e), "\n")
        return(NA)
      }
    )
  }

  # Bootstrap resampling
  transition_ages <- replicate(n_bootstrap, {
    indices_higher <- sample(nrow(data_higher), replace = TRUE)
    indices_lower <- sample(nrow(data_lower), replace = TRUE)
    calculate_transition_age(data_higher[indices_higher, ], data_lower[indices_lower, ], check_peak_dig = check_peak_dig)
  })
  
  # Calculate point estimate
  point_estimate <- median(transition_ages, na.rm = TRUE)
  
  # Calculate confidence intervals
  ci_lower <- quantile(transition_ages, (1 - confidence_level) / 2, na.rm = TRUE)
  ci_upper <- quantile(transition_ages, 1 - (1 - confidence_level) / 2, na.rm = TRUE)
  
  return(list(transition_age = point_estimate, ci_lower = ci_lower, ci_upper = ci_upper))
}

TRANSITION_AGE <- function(row){
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
                                      MaternalRank = row$MAT_STAT,
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
  PROP_pred$MAT_STAT <- row$MAT_STAT
  PROP_pred$RAIN <- row$Mean_RAIN

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
    group_by(TREATMENT, REC_AGE_D, SEX, WEIGHT, COMP, GS, RAIN, MAT_STAT, ID, Call_type) %>%
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
  
  transition_age <- data.frame(TREATMENT = character(),
                               SEX = character(),
                               WEIGHT = numeric(), 
                               COMP = numeric(),
                               GS = numeric(),
                               RAIN = numeric(),
                               MAT_STAT = character(),
                               DIG_REP = numeric(),
                               CC_DIG = numeric(),
                               ID = numeric(),
                               stringsAsFactors = FALSE)
  
  sumDIG_data_subset <- sumDIG_data %>% filter(TREATMENT == treatment, SEX == sex, WEIGHT == weight, COMP == comp, GS == gs, RAIN == rain, MAT_STAT == mat_stat, ID==id)
  sumREP_data_subset <- sumREP_data %>% filter(TREATMENT == treatment, SEX == sex, WEIGHT == weight, COMP == comp, GS == gs, RAIN == rain, MAT_STAT == mat_stat, ID==id)
  sumCC_data_subset <- sumCC_data %>% filter(TREATMENT == treatment, SEX == sex, WEIGHT == weight, COMP == comp, GS == gs, RAIN == rain, MAT_STAT == mat_stat, ID==id)

  rec_age_DIG_REP <- find_transition_age_CI(sumDIG_data_subset, sumREP_data_subset) 
  rec_age_CC_DIG <- find_transition_age_CI(sumCC_data_subset, sumDIG_data_subset, check_peak_dig = TRUE)
  rec_age_CC_REP <- find_transition_age_CI(sumCC_data_subset, sumREP_data_subset)
  
  transition_age <- rbind(transition_age, data.frame(TREATMENT = treatment,
                                                     SEX = sex,
                                                     WEIGHT = weight,
                                                     COMP = comp,
                                                     GS = gs,
                                                     RAIN = rain,
                                                     MAT_STAT = mat_stat,
                                                     DIG_REP = rec_age_DIG_REP,
                                                     CC_DIG = rec_age_CC_DIG,
                                                     CC_REP = rec_age_CC_REP,
                                                     ID = id))
  
  print(transition_age)
  
  return(transition_age)
}

#### Transition ages: SEMI & FULL ##############################################
# Calculate transition ages for all combinations of predictors #################

mean_weight <- mean(PROP_data$WEIGHT_DIFF_PER)
mean_comp_norm <- mean(PROP_data$COMP_NORM)
mean_gs <- mean(PROP_data$GROUPSIZE)
mean_rain <- mean(PROP_data$MonthlyRainfall)

# create a new dataset with unique individuals and their mean values
FLUT_ID_data <- PROP_data %>%
     group_by(ID) %>%
     summarise(Mean_WEIGHT_DIFF = mean_weight,
               Mean_COMP_NORM = mean_comp_norm,
               Mean_GROUPSIZE = mean_gs,
               Mean_RAIN = mean_rain,
               Mean_WEIGHT_z =0,
               Mean_COMP_z = 0,
               Mean_GS_z = 0,
               Mean_RAIN_z = 0,
               MAT_STAT = first(MaternalRank),
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
  transition_age <- TRANSITION_AGE(row = row)

  full_data <- rbind(full_data, transition_age)
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

rm(FLUT_ID_data)


### RESULTS in right format ####
transition_data <- full_data

# use CC > REP in case this happens AFTER CC > DIG to ensure realistic transition age
transition_data <- transition_data %>%
  mutate(
    CC_DIG.ci_lower = if_else(CC_REP.transition_age > CC_DIG.transition_age, CC_REP.ci_lower, CC_DIG.ci_lower),
    CC_DIG.ci_upper = if_else(CC_REP.transition_age > CC_DIG.transition_age, CC_REP.ci_upper, CC_DIG.ci_upper),
    CC_DIG.transition_age = if_else(CC_REP.transition_age > CC_DIG.transition_age, CC_REP.transition_age, CC_DIG.transition_age)
  )

# rename columns to semi-independence and full independence
transition_data <- transition_data %>%
  rename(SEMI.transition_age = DIG_REP.transition_age,
         SEMI.ci_lower = DIG_REP.ci_lower,
         SEMI.ci_upper = DIG_REP.ci_upper,
         FULL.transition_age = CC_DIG.transition_age,
         FULL.ci_lower = CC_DIG.ci_lower,
         FULL.ci_upper = CC_DIG.ci_upper)

# convert to long format
data_long <- transition_data %>%
  select(ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MAT_STAT, SEMI.transition_age, FULL.transition_age) %>%
  pivot_longer(cols = c(SEMI.transition_age, FULL.transition_age),
               names_to = "Variable", values_to = "Transition_age") %>%
  mutate(Variable = gsub("\\.transition_age", "", Variable)) %>%
  left_join(
    transition_data %>%
      select(ID, SEMI.ci_lower, FULL.ci_lower) %>%
      pivot_longer(cols = c(SEMI.ci_lower, FULL.ci_lower),
                   names_to = "Variable", values_to = "ci_lower") %>%
      mutate(Variable = gsub("\\.ci_lower", "", Variable))
  ) %>%
  left_join(
    transition_data %>%
      select(ID, SEMI.ci_upper, FULL.ci_upper) %>%
      pivot_longer(cols = c(SEMI.ci_upper, FULL.ci_upper),
                   names_to = "Variable", values_to = "ci_upper") %>%
      mutate(Variable = gsub("\\.ci_upper", "", Variable))
  )

data_long$TREATMENT <- as.factor(data_long$TREATMENT)
data_long$TREATMENT <- factor(data_long$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
data_long$SEX <- as.factor(data_long$SEX)
data_long$MAT_STAT <- as.factor(data_long$MAT_STAT)

data_long <- na.omit(data_long)

saveRDS(data_long, 'data/data_long_180_SAMPLE.rds')

### cleanup ###
rm(data_long, B_prop, PROP_priors, transition_data, PROP_data, full_data)
