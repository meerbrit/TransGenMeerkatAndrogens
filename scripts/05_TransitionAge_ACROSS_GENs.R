### Examine the transition ages for FLUT study subjects (F1 & F2) with their sample means ####
### XXX, 2024 #######################################################################

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
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

GEN_data <- readRDS('data/GEN_trans_data.rds')

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

################################################################################
################################################################################

#### ANOVAs ####
## REP -> DIG: SEMI ####
SEMI_data <- GEN_data %>%
  filter(Variable == "SEMI") %>%
  na.omit()
SEMI_data$TREATMENT <- factor(SEMI_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
SEMI_data$GEN <- as.factor(SEMI_data$GEN)
SEMI_data$MAT_STAT <- factor(SEMI_data$MAT_STAT, levels = c('DOM', 'SUB'))
SEMI_data$MAT_STAT_in_F2 <- ifelse(SEMI_data$GEN == 2, as.character(SEMI_data$MAT_STAT), 'F1')
SEMI_data$MAT_STAT_in_F2 <- factor(SEMI_data$MAT_STAT_in_F2, levels = c("F1", "DOM", "SUB"))

SEMI_data$WEIGHT_z <- scale(SEMI_data$WEIGHT)
SEMI_data$COMP_z <- scale(SEMI_data$COMP)
SEMI_data$GS_z <- scale(SEMI_data$GS)
SEMI_data$RAIN_z <- scale(SEMI_data$RAIN)

mean(FULL_data$Transition_age)
sd(FULL_data$Transition_age)

priors_student_horse <- c(set_prior("normal(85,20)", class = "Intercept"), set_prior("horseshoe(1)", class = "b")) #
SEMI_anova_final <- brm(formula = Transition_age| trunc(lb = 1, ub = 180)  ~ (TREATMENT * GEN * SEX) + (MAT_STAT_in_F2 * TREATMENT * SEX) +
                          WEIGHT_z + COMP_z + I(COMP_z^2) + GS_z +I(GS_z^2) +
                          RAIN_z +  I(RAIN_z^2) , 
                        data = SEMI_data,
                        family = student(link = 'identity'),
                        chains = 4, iter = 5000, warmup = 1500, seed = 23425235, control = list(max_treedepth = 20, adapt_delta=0.99),
                        save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                        prior = priors_student_horse,
                        threads = threading(4),
                        file="SEMI_anova_ACROSS_final")

rm(priors_student_horse)


SEMI_anova <- readRDS("models/SEMI_anova_ACROSS_final.rds")

plot(SEMI_anova)
pp_check(SEMI_anova, ndraws=100)

loo_R2(SEMI_anova, moment_match = T)

bayes_R2(SEMI_anova)

summary(SEMI_anova)

describe_posterior(
  SEMI_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(SEMI_anova),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

#### EMMs ####
### Are sex/treatment differences consistent across GENs ####
(treat_gen_sex <- emmeans(SEMI_anova,  ~ TREATMENT:GEN:SEX, nesting = c("MAT_STAT_in_F2 %in% GEN")))

pd(treat_gen_sex)

equivalence_test(treat_gen_sex, range = rope_range(SEMI_anova))

p_significance(treat_gen_sex, threshold = rope_range(SEMI_anova))

# Pairwise contrasts by GEN for each sex
(tgs_contrast <- contrast(treat_gen_sex, interaction = "pairwise", by = c('GEN', 'SEX')))

pd(tgs_contrast)


equivalence_test(tgs_contrast, range = rope_range(SEMI_anova))


p_significance(tgs_contrast, threshold = rope_range(SEMI_anova))


# Compare treatment effects for each sex between GEN1 and GEN2
(tgs_contrast_gen_sex <- contrast(treat_gen_sex, interaction = "pairwise", by = c('TREATMENT', 'SEX')))

pd(tgs_contrast_gen_sex)

equivalence_test(tgs_contrast_gen_sex, range = rope_range(SEMI_anova))

p_significance(tgs_contrast_gen_sex, threshold = rope_range(SEMI_anova))

# Plot TGS
treat_gen_sex <- emmeans(SEMI_anova, ~ TREATMENT:GEN:SEX, nesting = c("MAT_STAT_in_F2 %in% GEN"))
emm_data <- as.data.frame(treat_gen_sex)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", 'FLUT'))
emm_data$GEN <- as.factor(emm_data$GEN)
levels(emm_data$GEN) <- c('F1', 'F2')

# 700* 500: SEMI_TGS
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Sex')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Offspring", y = "Age (days)\n") +
  geom_point(data = SEMI_data, aes(x = TREATMENT, y = Transition_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()+
  facet_wrap(~GEN)

rm(treat_sex)

### Check the influence of the maternal status in F2 on the sexes #### NOT POSS BC OF MISSING INTERACTION IN MODEL?
(F2_only <- emmeans(SEMI_anova,  ~ MAT_STAT_in_F2:TREATMENT:SEX, at = list(GEN = 2, MAT_STAT_in_F2=c('DOM', 'SUB')), nesting = c("MAT_STAT_in_F2 %in% GEN")))

pd(F2_only)

equivalence_test(F2_only, range = rope_range(SEMI_anova))

p_significance(F2_only, threshold = rope_range(SEMI_anova))

# Compare DOM vs. SUB within each treatment and sex 
(f2_contrast <- contrast(F2_only, interaction = "pairwise", by = c("TREATMENT", 'GEN', 'SEX')))

pd(f2_contrast)

equivalence_test(f2_contrast, range = rope_range(SEMI_anova))


p_significance(f2_contrast, threshold = rope_range(SEMI_anova))

# compare the sexes within treatment and mat_stat
(f2_contrast <- contrast(F2_only, interaction = "pairwise", by = c("TREATMENT", 'GEN', 'MAT_STAT_in_F2')))

pd(f2_contrast)

equivalence_test(f2_contrast, range = rope_range(SEMI_anova))

p_significance(f2_contrast, threshold =  rope_range(SEMI_anova))

# Plot TMS
treat_mat_sex <- emmeans(SEMI_anova, ~ MAT_STAT_in_F2:TREATMENT:SEX, at = list(GEN = 2, MAT_STAT_in_F2=c('DOM', 'SUB')), nesting = c("MAT_STAT_in_F2 %in% GEN"))
emm_data <- as.data.frame(treat_mat_sex)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", 'FLUT'))

# no F1
f2_data <- SEMI_data %>%
  filter(GEN == 'F2')

# 700* 500: SEMI_F2_TMS
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Sex')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Offspring", y = "Age (days)\n") +
  geom_point(data = f2_data, aes(x = TREATMENT, y = Transition_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()+
  facet_wrap(~MAT_STAT_in_F2)

rm(treat_mat_sex)

### PLOTS: SEMI ####
# coefficient plots #### 
posterior_desc <- describe_posterior(
  SEMI_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(SEMI_anova),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 32),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'ICOMP_zE2', 'Competition^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IRAIN_zE2', 'Monthly rainfall^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GEN2', 'F2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STAT_in_F2SUB', 'F2 Maternal status = SUB')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STAT_in_F2DOM', 'F2 Maternal status = DOM')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STAT_in_F2F1', 'First Gen')

custom_order <- c("DT:Male:F2 Maternal status = SUB","SC:Male:F2 Maternal status = SUB","Male:F2 Maternal status = SUB",
                  "DT:Male:F2 Maternal status = DOM","SC:Male:F2 Maternal status = DOM","Male:F2 Maternal status = DOM",
                  'DT:F2:Male','SC:F2:Male','F2:Male',
                  'DT:F2','SC:F2','F2',
                  'Monthly rainfall^2', 
                  'Monthly rainfall', 
                  'Group size^2', 
                  'Group size', 
                  'Competition^2', 
                  'Competition',
                  'Weight offset', 
                  "DT:F2 Maternal status = SUB", "SC:F2 Maternal status = SUB", "F2 Maternal status = SUB", 
                  "DT:F2 Maternal status = DOM", "SC:F2 Maternal status = DOM", "F2 Maternal status = DOM" ,
                  'DT:Male','SC:Male',"Male", 
                  'DT','SC')


posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_SEMI 800*900
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1,3), name = "Offspring", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  #ggtitle('Full dependence to semi-dependence')+
  theme_clean()+
  theme(legend.position="none")

## DIG -> CC ####

FULL_data <- GEN_data %>%
  filter(Variable == "FULL") %>%
  na.omit()
FULL_data$TREATMENT <- factor(FULL_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
FULL_data$GEN <- as.factor(FULL_data$GEN)
FULL_data$MAT_STAT <- factor(FULL_data$MAT_STAT, levels = c('DOM', 'SUB'))
FULL_data$MAT_STAT_in_F2 <- ifelse(FULL_data$GEN == 2, as.character(FULL_data$MAT_STAT), 'F1')
FULL_data$MAT_STAT_in_F2 <- factor(FULL_data$MAT_STAT_in_F2, levels = c("F1", "DOM", "SUB"))

FULL_data$WEIGHT_z <- scale(FULL_data$WEIGHT)
FULL_data$COMP_z <- scale(FULL_data$COMP)
FULL_data$GS_z <- scale(FULL_data$GS)
FULL_data$RAIN_z <- scale(FULL_data$RAIN)

mean(FULL_data$Transition_age)
sd(FULL_data$Transition_age)

priors_student_horse <- c(set_prior("normal(120,15)", class = "Intercept"), set_prior("horseshoe(1)", class = "b")) #
FULL_anova_final <- brm(formula = Transition_age| trunc(lb = 1, ub = 180)  ~ (TREATMENT * GEN * SEX) + (MAT_STAT_in_F2 * TREATMENT * SEX) +
                          WEIGHT_z + COMP_z + I(COMP_z^2) + GS_z +I(GS_z^2) +
                          RAIN_z +  I(RAIN_z^2) , 
                        data = FULL_data,
                        family = student(link = 'identity'),
                        chains = 4, iter = 5000, warmup = 1500, seed = 23425235, control = list(max_treedepth = 20, adapt_delta=0.99),
                        save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                        prior = priors_student_horse,
                        threads = threading(4),
                        file="FULL_anova_ACROSS_final")
rm(priors_student_horse)

FULL_anova <- readRDS("models/FULL_anova_ACROSS_final.rds")

plot(FULL_anova)

pp_check(FULL_anova, ndraws=100)

loo_R2(FULL_anova, moment_match=T)

bayes_R2(FULL_anova)

summary(FULL_anova)

describe_posterior(
  FULL_anova,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(FULL_anova),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

#### EMMs ####
### Are sex/treatment differences consistent across GENs ####
(treat_gen_sex <- emmeans(FULL_anova,  ~ TREATMENT:GEN:SEX, nesting = c("MAT_STAT_in_F2 %in% GEN")))

pd(treat_gen_sex)

equivalence_test(treat_gen_sex, range = rope_range(FULL_anova))

p_significance(treat_gen_sex, threshold = rope_range(FULL_anova))

# Pairwise contrasts by GEN for each sex
(tgs_contrast <- contrast(treat_gen_sex, interaction = "pairwise", by = c('GEN', 'SEX')))

pd(tgs_contrast)


equivalence_test(tgs_contrast, range = rope_range(FULL_anova))


p_significance(tgs_contrast, threshold = rope_range(FULL_anova))

# Compare treatment effects for each sex between GEN1 and GEN2
(tgs_contrast_gen_sex <- contrast(treat_gen_sex, interaction = "pairwise", by = c('TREATMENT', 'SEX')))

pd(tgs_contrast_gen_sex)

equivalence_test(tgs_contrast_gen_sex, range = rope_range(FULL_anova))

p_significance(tgs_contrast_gen_sex, threshold = rope_range(FULL_anova))

# Plot TGS
treat_gen_sex <- emmeans(FULL_anova, ~ TREATMENT:GEN:SEX, nesting = c("MAT_STAT_in_F2 %in% GEN"))
emm_data <- as.data.frame(treat_gen_sex)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", 'FLUT'))
emm_data$GEN <- as.factor(emm_data$GEN)
levels(emm_data$GEN) <- c('F1', 'F2')
FULL_data$GEN <- as.factor(FULL_data$GEN)
levels(FULL_data$GEN) <- c('F1', 'F2')

# 700* 500: FULL_TGS
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Sex')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Offspring", y = "Age (days)\n") +
  geom_point(data = FULL_data, aes(x = TREATMENT, y = Transition_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()+
  facet_wrap(~GEN)

rm(treat_gen_sex)

### Check the influence of the maternal status in F2 on the sexes ####
(F2_only <- emmeans(FULL_anova,  ~ MAT_STAT_in_F2:TREATMENT:SEX, at = list(GEN = 2, MAT_STAT_in_F2=c('DOM', 'SUB')), nesting = c("MAT_STAT_in_F2 %in% GEN")))

pd(F2_only)

equivalence_test(F2_only, range = rope_range(FULL_anova))

p_significance(F2_only, threshold = rope_range(FULL_anova))

# Compare DOM vs. SUB within each treatment and sex 
(f2_contrast <- contrast(F2_only, interaction = "pairwise", by = c("TREATMENT", 'GEN', 'SEX')))

pd(f2_contrast)

equivalence_test(f2_contrast, range = rope_range(FULL_anova))

p_significance(f2_contrast, threshold = rope_range(FULL_anova))

# compare the sexes within treatment and mat_stat
(f2_contrast <- contrast(F2_only, interaction = "pairwise", by = c("TREATMENT", 'GEN', 'MAT_STAT_in_F2')))

pd(f2_contrast)

equivalence_test(f2_contrast, range = rope_range(FULL_anova))

p_significance(f2_contrast, threshold =  rope_range(FULL_anova))

# Plot TMS
treat_mat_sex <- emmeans(FULL_anova, ~ MAT_STAT_in_F2:TREATMENT:SEX, at = list(GEN = 2, MAT_STAT_in_F2=c('DOM', 'SUB')), nesting = c("MAT_STAT_in_F2 %in% GEN"))
emm_data <- as.data.frame(treat_mat_sex)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", 'FLUT'))

# no F1
f2_data <- FULL_data %>%
  filter(GEN == 'F2')

# 700* 500: FULL_TMS
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Sex')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Offspring", y = "Age (days)\n") +
  geom_point(data = f2_data, aes(x = TREATMENT, y = Transition_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()+
  facet_wrap(~MAT_STAT_in_F2)

rm(treat_mat_sex)

#### PLOTS: FULL ####
# coefficient plots #### 
posterior_desc <- describe_posterior(
  FULL_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(FULL_anova),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 32),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'ICOMP_zE2', 'Competition^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IRAIN_zE2', 'Monthly rainfall^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GEN2', 'F2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STAT_in_F2SUB', 'F2 Maternal status = SUB')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STAT_in_F2DOM', 'F2 Maternal status = DOM')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STAT_in_F2F1', 'First Gen')

custom_order <- c("DT:Male:F2 Maternal status = SUB","SC:Male:F2 Maternal status = SUB","Male:F2 Maternal status = SUB",
                  "DT:Male:F2 Maternal status = DOM","SC:Male:F2 Maternal status = DOM","Male:F2 Maternal status = DOM",
                  'DT:F2:Male','SC:F2:Male','F2:Male',
                  'DT:F2','SC:F2','F2',
                  'Monthly rainfall^2', 
                  'Monthly rainfall', 
                  'Group size^2', 
                  'Group size', 
                  'Competition^2', 
                  'Competition',
                  'Weight offset', 
                  "DT:F2 Maternal status = SUB", "SC:F2 Maternal status = SUB", "F2 Maternal status = SUB", 
                  "DT:F2 Maternal status = DOM", "SC:F2 Maternal status = DOM", "F2 Maternal status = DOM" ,
                  'DT:Male','SC:Male',"Male", 
                  'DT','SC')


posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_FULL 800*900
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  # ggtitle('Semi-dependence to full independence')+
  theme_clean()+
  theme(legend.position="none")

rm(SEMI_data, FULL_data, SEMI_anova, FULL_anova, cond_data)

### cleanup ###
rm(GEN_data, B_prop, priors, transition_data, GEN_data)
