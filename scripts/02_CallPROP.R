#### Analysis script for 2nd GEN Flutamide study: CALL PROPORTIONS #####################

# Bayesian Multilevel models: logistic regression including: 
# Treatment (T), age (A), sex (S), weight offset (W), 
# also checking for (normalised) competition score (C = H/P ratio & pups & adults),
# group size and rainfall
#
# XXX, 2024 

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(ggplot2) 
library(ggpubr) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # 
library(extrafont)# use font_import() on first use

set.seed(23)

PROP_data <- readRDS('data/PROP_data.rds')

# load priors based on posterior of the first gen:
PROP_priors <- readRDS("data/PROP_priors.rds")

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(2)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(2), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(2)),
      legend.text = element_text(face = 'italic', size = rel(1.5)))
}

get_age_vars <- function(){
  age_vars <- c((30- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D), 
                (45- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D), 
                (60- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D), 
                (75- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D), 
                (90- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D),
                (105- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D), 
                (120- mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D))
  return(age_vars)
}

get_org_value_for_z <- function(z_value, column){
  return ((z_value*sd(column)) + mean(column))   
}

get_z_value_for_org <- function(org_value, column){
  return ((org_value - mean(column))/sd(column))   
}


# y distributions -RAW
plot(density(PROP_data$Sum_BEG/PROP_data$Total_calls), 
     xlab = "REP proportions",
     ylab = "Density")

plot(density(PROP_data$Sum_DIG/PROP_data$Total_calls), 
     xlab = "DIG proportions",
     ylab = "Density")

plot(density(PROP_data$Sum_CC/PROP_data$Total_calls), 
     xlab = "CC proportions",
     ylab = "Density")


########################## BAYES MODELS ########################################
# # FULL MULTIVARIATE MODEL ####
# # model to answer: MAS / TAS / TAM
# # x ~ (M*A*S + T*A*S + T*A*M) + W + C + G + R + (1|LITTER/ID)
bf_REP <- bf(Sum_BEG | trials(Total_calls) ~ (TREATMENT*AGE_z*SEX + MaternalRank*AGE_z*SEX + TREATMENT*AGE_z*MaternalRank) +
                WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) +
               GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID))

bf_DIG <-bf(Sum_DIG | trials(Total_calls) ~ (TREATMENT*AGE_z*SEX + MaternalRank*AGE_z*SEX + TREATMENT*AGE_z*MaternalRank) +
              (TREATMENT*I(AGE_z^2)*SEX + MaternalRank*I(AGE_z^2)*SEX + TREATMENT*I(AGE_z^2)*MaternalRank) +
              WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) + GS_z + RAIN_z + I(RAIN_z^2) + (1|LITTER_CODE/ID))

bf_CC <- bf(Sum_CC | trials(Total_calls) ~ (TREATMENT*AGE_z*SEX + MaternalRank*AGE_z*SEX + TREATMENT*AGE_z*MaternalRank) +
              WEIGHT_z + COMP_NORM_z  + I(COMP_NORM_z^2) + GS_z + I(GS_z^2) + RAIN_z + (1|LITTER_CODE/ID))

multivar_formula <- mvbrmsformula(bf_REP, bf_DIG, bf_CC)

B_prop <- brms::brm(formula = multivar_formula,
                    data = PROP_data, family = zero_inflated_beta_binomial(link='logit'),
                    chains = 4, iter = 10000, warmup = 2500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.999),
                    save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                    threads = threading(4),
                    prior = PROP_priors,
                    file="B_prop")

rm(bf_REP, bf_DIG, bf_CC, multivar_formula, PROP_priors)

B_prop <- readRDS('B_prop.rds')

#### RESULTS: MULTIVAR model ####
summary(B_prop)

plot(B_prop)
pp_check(B_prop, ndraws = 100, resp='SumBEG')
pp_check(B_prop, ndraws = 100, resp='SumDIG')
pp_check(B_prop, ndraws = 100, resp='SumCC')

describe_posterior(
  B_prop$fit,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = c(-0.18, 0.18),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_prop, moment_match=T)


bayes_R2(B_prop)

# coefficient plots ####
posterior_desc <- describe_posterior(
  B_prop$fit,
  effects = "fixed",
  component = "all",
  rope_range = c(-0.18, 0.18),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

REP_desc <- posterior_desc %>%
    filter(str_detect(Parameter, "SumBEG"))
DIG_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumDIG"))
CC_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumCC"))

REP_desc <- REP_desc[(2:26),]
DIG_desc <- DIG_desc[(2:36),]
CC_desc <- CC_desc[c(2:26),]

# REP Coeff
# clean up labels:
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'b_SumBEG_', '')
REP_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", REP_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", REP_desc$Parameter), "DT", "DC"))
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTFLUT', 'DT')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTSUB', 'SC')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'WEIGHT_z', 'Weight offset')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'IGS_zE2', 'Group size^2')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'GS_z', 'Group size')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'AGE_z', 'Age')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'ICOMP_NORM_zE2', 'Competition^2')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'COMP_NORM_z', 'Competition')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'SEXM', 'Male')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'MaternalRankSUB', 'SUB mother')

custom_order <- c('Age:Male:SUB mother',
                  'Male:SUB mother', 
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall', 
                  'DT:Group size^2','SC:Group size^2','Group size^2', 
                  'DT:Group size','SC:Group size','Group size', 
                  'DT:Competition^2','SC:Competition^2','Competition^2',
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Weight offset','SC:Weight offset','Weight offset', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age','SC:Age',"Age", 
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT','SC') 

# Update the order of TREATMENT factor levels
REP_desc$TREATMENT <- factor(REP_desc$TREATMENT, levels = c("DC", "SC", "DT"))
REP_desc$Parameter <- factor(REP_desc$Parameter, levels = custom_order)
# Coeff_REP 700 * 800
ggplot(REP_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  ggtitle('REP proportion')+
  theme_clean()+
  theme(legend.position="none")

# DIG coeff
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'b_SumDIG_', '')
DIG_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", DIG_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", DIG_desc$Parameter), "DT", "DC"))
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTFLUT', 'DT')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTSUB', 'SC')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'WEIGHT_z', 'Weight offset')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IAGE_zE2', 'Age^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'AGE_z', 'Age')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'GS_z', 'Group size')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'ICOMP_NORM_zE2', 'Competition^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'COMP_NORM_z', 'Competition')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'SEXM', 'Male')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IRAIN_zE2', 'Monthly rainfall^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c("Male:SUB mother:Age^2", 
                  "Age:Male:SUB mother" ,
                  'DT:Male:SUB mother','SC:Male:SUB mother',"Male:SUB mother",
                  'DT:SUB mother:Age^2','SC:SUB mother:Age^2',"SUB mother:Age^2", 
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:Male:Age^2','SC:Male:Age^2','Male:Age^2',
                  'DT:Age:Male','SC:Age:Male','Age:Male',
                  'Monthly rainfall^2',
                  'Monthly rainfall',
                  'Group size', 
                  'Competition^2',
                  'Competition', 
                  'Weight offset', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age",
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT','SC') 

DIG_desc$Parameter <- factor(DIG_desc$Parameter, levels = custom_order)
DIG_desc$TREATMENT <- factor(DIG_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG 700*800
ggplot(DIG_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  ggtitle('DIG proportion')+
  theme_clean()+
  theme(legend.position="none")

# CC Coeff
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'b_SumCC_', '')
CC_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", CC_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", CC_desc$Parameter), "DT", "DC"))
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTFLUT', 'DT')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTSUB', 'SC')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'WEIGHT_z', 'Weight offset')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'IGS_zE2', 'Group size^2')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'GS_z', 'Group size')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'AGE_z', 'Age')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'ICOMP_NORM_zE2', 'Competition^2')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'COMP_NORM_z', 'Competition')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'SEXM', 'Male')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'IRAIN_zE2', 'Monthly rainfall^2')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Age:Male:SUB mother',
                  'Male:SUB mother', 
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'Monthly rainfall', 
                  'Group size^2', 
                  'Group size', 
                  'Competition^2',
                  'Competition', 
                  'Weight offset', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age','SC:Age',"Age", 
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT','SC') 

CC_desc$Parameter <- factor(CC_desc$Parameter, levels = custom_order)
CC_desc$TREATMENT <- factor(CC_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_CC 700*800
ggplot(CC_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  ggtitle('CC proportion')+
  theme_clean()+
  theme(legend.position="none")

rm(REP_desc, DIG_desc, CC_desc, posterior_desc, custom_order)

### EMMs TASM ####
### REP ####
REP_ROPE <- c(-0.18, 0.18)
(mat_stat_treat_sex <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank:SEX, var="AGE_z", resp='SumBEG'))

pd(mat_stat_treat_sex)

equivalence_test(mat_stat_treat_sex, range = REP_ROPE)

p_significance(mat_stat_treat_sex, threshold = REP_ROPE)

### DIG ####
DIG_ROPE <- c(-0.18, 0.18)
(mat_stat_treat_sex <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank:SEX, var="AGE_z", resp='SumDIG'))

pd(mat_stat_treat_sex)

equivalence_test(mat_stat_treat_sex, range = DIG_ROPE)

p_significance(mat_stat_treat_sex, threshold = DIG_ROPE)

### CC ####
CC_ROPE <- c(-0.18, 0.18)
(mat_stat_treat_sex <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank:SEX, var="AGE_z", resp='SumCC'))

pd(mat_stat_treat_sex)

equivalence_test(mat_stat_treat_sex, range = CC_ROPE)

p_significance(mat_stat_treat_sex, threshold = CC_ROPE)

rm(mat_stat_treat_sex)

#### MODEL ONTOGENY PLOTS ####
B_prop <- readRDS('B_prop.rds')

# get all needed values
sd_age <- sd(PROP_data$REC_AGE_D)
mean_age <- mean(PROP_data$REC_AGE_D)

range(PROP_data$REC_AGE_D)# 39 130

rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(PROP_data$REC_AGE_D))/sd(PROP_data$REC_AGE_D)
# 1 day steps not needed here 
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

#predictions based on mean values
PROP_pred <- B_prop %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(PROP_data$TREATMENT),
                                    SEX = levels(PROP_data$SEX),
                                    MaternalRank = levels(PROP_data$MaternalRank),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(PROP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(PROP_data$COMP_NORM_z),
                                    GS_z = mean(PROP_data$GS_z),
                                    RAIN_z = mean(PROP_data$RAIN_z),
                                    Total_calls=1), re_formula = NA,  robust = T)

#unscale AGE_z values:
PROP_pred$REC_AGE_D <- PROP_pred$AGE_z * sd_age + mean_age
# ensure right format
PROP_pred$Call_prop <- PROP_pred$.epred
PROP_pred$Call_type <- as.factor(PROP_pred$.category)
PROP_pred$SEX <- as.factor(PROP_pred$SEX)
PROP_pred$MaternalRank <- as.factor(PROP_pred$MaternalRank)
PROP_pred$TREATMENT <- factor(PROP_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
PROP_data$TREATMENT <- factor(PROP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

# normalise predictions to ensure they sum up to 1
PROP_pred_sum <- PROP_pred %>%
  group_by(REC_AGE_D, TREATMENT, SEX, MaternalRank, .draw) %>%
  mutate(Sum_prop = sum(Call_prop))

PROP_pred_normalized <- PROP_pred_sum %>%
  mutate(Normalized_prop = Call_prop / Sum_prop)

rm(PROP_pred_sum)

#### TASM ####
#900*500: TAMS (_re/_na)
ggplot(PROP_pred_normalized, aes(x = REC_AGE_D, y = Normalized_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +
  stat_lineribbon(.width = .95) + 
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), name = 'Call type', labels = c('REP', 'DIG', 'CC'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid(SEX ~ MaternalRank, 
             labeller = labeller(
               MaternalRank = c('DOM' = 'Dominant mother', 'SUB' = 'Subordinate mother'),
               SEX = c('F' = 'Female', 'M' = 'Male')
             ))

# just REP TAMS
REP_pred <- subset(PROP_pred_normalized, Call_type == 'SumBEG')
REP_pred$REP_prop <- REP_pred$Normalized_prop

#900*600: REP_TAMS_
ggplot(REP_pred, aes(x = REC_AGE_D, y = REP_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) + 
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid(SEX ~ MaternalRank, 
             labeller = labeller(
               MaternalRank = c('DOM' = 'Dominant mother', 'SUB' = 'Subordinate mother'),
               SEX = c('F' = 'Female', 'M' = 'Male')
             ))

# different layout
#900*600: REP_TAMS_
ggplot(REP_pred, aes(x = REC_AGE_D, y = REP_prop, color = TREATMENT, fill = TREATMENT, linetype= MaternalRank)) +
  stat_lineribbon(.width = .95) + 
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Maternal\nstatus', labels = c('dominant', 'subordinate'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid( ~ SEX, 
             labeller = labeller(
               SEX = c('F' = 'Female', 'M' = 'Male')
             ))


# just DIG TAMS
DIG_pred <- subset(PROP_pred_normalized, Call_type == 'SumDIG')
DIG_pred$DIG_prop <- DIG_pred$Normalized_prop

#800*500: DIG_TAMS_
ggplot(DIG_pred, aes(x = REC_AGE_D, y = DIG_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) + 
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid(SEX ~ MaternalRank, 
             labeller = labeller(
               MaternalRank = c('DOM' = 'Dominant mother', 'SUB' = 'Subordinate mother'),
               SEX = c('F' = 'Female', 'M' = 'Male')
             ))

# different layout
ggplot(DIG_pred, aes(x = REC_AGE_D, y = DIG_prop, color = TREATMENT, fill = TREATMENT, linetype= MaternalRank)) +
  stat_lineribbon(.width = .95) + 
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Maternal\nstatus', labels = c('dominant', 'subordinate'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid( ~ SEX, 
              labeller = labeller(
                SEX = c('F' = 'Female', 'M' = 'Male')
              ))


# just CC TAMS
CC_pred <- subset(PROP_pred_normalized, Call_type == 'SumCC')
CC_pred$CC_prop <- CC_pred$Normalized_prop

#900*600: CC_TAMS
ggplot(CC_pred, aes(x = REC_AGE_D, y = CC_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) + 
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid(SEX ~ MaternalRank, 
             labeller = labeller(
               MaternalRank = c('DOM' = 'Dominant mother', 'SUB' = 'Subordinate mother'),
               SEX = c('F' = 'Female', 'M' = 'Male')
             ))

# different layout
ggplot(CC_pred, aes(x = REC_AGE_D, y = CC_prop, color = TREATMENT, fill = TREATMENT, linetype= MaternalRank)) +
  stat_lineribbon(.width = .95) + 
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Maternal\nstatus', labels = c('dominant', 'subordinate'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n",  
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_grid( ~ SEX, 
              labeller = labeller(
                SEX = c('F' = 'Female', 'M' = 'Male')
              ))

#clean up
rm(REP_pred, DIG_pred, CC_pred, REP_DIG_pred, DIG_CC_pred, PROP_pred, PROP_data, 
   age_z_vals, rec_age_c, mean_age, sd_age, PROP_pred_normalized, B_prop)

