### Compare DIG peak age for FLUT F1 & F2 study subjects and include all predictors ####
# Check influence of maternal status in F2 #############################################
### XXX, 2024 ##########################################################################

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(writexl)
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

GEN_data <- readRDS('data/DIG_PEAK_GEN_data.rds')

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
# ANOVA ####
################################################################################

DIG_PEAK_data <- GEN_data %>%
  filter(Variable == "DIG") %>%
  na.omit()
DIG_PEAK_data$TREATMENT <- factor(DIG_PEAK_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_PEAK_data$GEN <- as.factor(DIG_PEAK_data$GEN)
DIG_PEAK_data$MAT_STAT <- factor(DIG_PEAK_data$MAT_STAT, levels = c('DOM', 'SUB'))
DIG_PEAK_data$MAT_STAT_in_F2 <- ifelse(DIG_PEAK_data$GEN == 2, as.character(DIG_PEAK_data$MAT_STAT), 'F1')
DIG_PEAK_data$MAT_STAT_in_F2 <- factor(DIG_PEAK_data$MAT_STAT_in_F2, levels = c("F1", "DOM", "SUB"))

DIG_PEAK_data$WEIGHT_z <- scale(DIG_PEAK_data$WEIGHT)
DIG_PEAK_data$COMP_z <- scale(DIG_PEAK_data$COMP)
DIG_PEAK_data$GS_z <- scale(DIG_PEAK_data$GS)
DIG_PEAK_data$RAIN_z <- scale(DIG_PEAK_data$RAIN)

priors_student_horse <- c(set_prior("normal(95,15)", class = "Intercept"), set_prior("horseshoe(1)", class = "b")) #
DIG_PEAK_anova_final <- brm(formula = Peak_age| trunc(lb = 1, ub = 180)  ~ (TREATMENT * GEN * SEX) + (MAT_STAT_in_F2 * TREATMENT * SEX) +
                                  WEIGHT_z + COMP_z + I(COMP_z^2) + GS_z +I(GS_z^2) +
                                  RAIN_z +  I(RAIN_z^2) , 
                                data = DIG_PEAK_data,
                                family = student(link = 'identity'),
                                chains = 4, iter = 5000, warmup = 1500, seed = 23425235, control = list(max_treedepth = 20, adapt_delta=0.99),
                                save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                                prior = priors_student_horse,
                                threads = threading(4),
                                file="DIG_PEAK_anova_ACROSS_final")

DIG_PEAK_anova <- readRDS("DIG_PEAK_anova_ACROSS_final.rds")

pp_check(DIG_PEAK_anova, ndraws=100)
loo_R2(DIG_PEAK_anova, moment_match=T)
bayes_R2(DIG_PEAK_anova)

summary(DIG_PEAK_anova)

describe_posterior(
  DIG_PEAK_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(DIG_PEAK_anova),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

#### EMMs ####
### Are sex/treatment differences consistent across GENs ####
(treat_gen_sex <- emmeans(DIG_PEAK_anova,  ~ TREATMENT:GEN:SEX, nesting = c("MAT_STAT_in_F2 %in% GEN")))

pd(treat_gen_sex)

equivalence_test(treat_gen_sex, range = rope_range(DIG_PEAK_anova))

p_significance(treat_gen_sex, threshold = rope_range(DIG_PEAK_anova))

# Pairwise contrasts by GEN for each sex
(tgs_contrast <- contrast(treat_gen_sex, interaction = "pairwise", by = c('GEN', 'SEX')))

pd(tgs_contrast)

equivalence_test(tgs_contrast, range = rope_range(DIG_PEAK_anova))

p_significance(tgs_contrast, threshold = rope_range(DIG_PEAK_anova))


# Compare treatment effects for each sex between GEN1 and GEN2
(tgs_contrast_gen_sex <- contrast(treat_gen_sex, interaction = "pairwise", by = c('TREATMENT', 'SEX')))

pd(tgs_contrast_gen_sex)

equivalence_test(tgs_contrast_gen_sex, range = rope_range(DIG_PEAK_anova))

p_significance(tgs_contrast_gen_sex, threshold = rope_range(DIG_PEAK_anova))

# Plot TGS
treat_gen_sex <- emmeans(DIG_PEAK_anova, ~ TREATMENT:GEN:SEX, nesting = c("MAT_STAT_in_F2 %in% GEN"))
emm_data <- as.data.frame(treat_gen_sex)
emm_data$PEAK_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", 'FLUT'))
emm_data$GEN <- as.factor(emm_data$GEN)
levels(emm_data$GEN) <- c('F1', 'F2')
levels(DIG_PEAK_data$GEN) <- c('F1', 'F2')

# 700* 500: PEAK_DIG_TGS
ggplot(emm_data, aes(x = TREATMENT, y = PEAK_age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Sex')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Offspring", y = "Age (days)\n") +
  geom_point(data = DIG_PEAK_data, aes(x = TREATMENT, y = Peak_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()+
  facet_wrap(~GEN)

rm(treat_gen_sex)

### Check the influence of the maternal status in F2 on the sexes ####
(F2_only <- emmeans(DIG_PEAK_anova,  ~ MAT_STAT_in_F2:TREATMENT:SEX, at = list(GEN = 2, MAT_STAT_in_F2=c('DOM', 'SUB')), nesting = c("MAT_STAT_in_F2 %in% GEN")))

pd(F2_only)

equivalence_test(F2_only, range = rope_range(DIG_PEAK_anova))

p_significance(F2_only, threshold = rope_range(DIG_PEAK_anova))


# Compare DOM vs. SUB within each treatment and sex << !!!! CANNOT DETECT SEX DIFF BC OF FORMULA ???
(f2_contrast <- contrast(F2_only, interaction = "pairwise", by = c("TREATMENT", 'GEN', 'SEX')))

pd(f2_contrast)

equivalence_test(f2_contrast, range = rope_range(DIG_PEAK_anova))

p_significance(f2_contrast, threshold = rope_range(DIG_PEAK_anova))

# compare the sexes within treatment and mat_stat
(f2_contrast <- contrast(F2_only, interaction = "pairwise", by = c("TREATMENT", 'GEN', 'MAT_STAT_in_F2')))

pd(f2_contrast)

equivalence_test(f2_contrast, range = rope_range(DIG_PEAK_anova))

p_significance(f2_contrast, threshold =  rope_range(DIG_PEAK_anova))

# Plot TMS
treat_mat_sex <- emmeans(DIG_PEAK_anova, ~ MAT_STAT_in_F2:TREATMENT:SEX, at = list(GEN = 2, MAT_STAT_in_F2=c('DOM', 'SUB')), nesting = c("MAT_STAT_in_F2 %in% GEN"))
emm_data <- as.data.frame(treat_mat_sex)
emm_data$PEAK_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", 'FLUT'))

# no F1
f2_data <- DIG_PEAK_data %>%
  filter(GEN == 'F2')

# 700* 500: PEAK_DIG_TMS
ggplot(emm_data, aes(x = TREATMENT, y = PEAK_age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Sex')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Offspring", y = "Age (days)\n") +
  geom_point(data = f2_data, aes(x = TREATMENT, y = Peak_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()+
  facet_wrap(~MAT_STAT_in_F2)

rm(treat_mat_sex)

# coefficient plot ####
posterior_desc <- describe_posterior(
  DIG_PEAK_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(DIG_PEAK_anova),
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

# Coeff_PEAK_DIG 800*900
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1,3), name = "Offspring", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()+
  theme(legend.position="none")


rm(GEN_data, full_data, posterior_desc, GEN_data, PROP_priors,
   DIG_data, DIG_PEAK_anova, emms_peak_FULL, PEAK_full_contrasts_rank, PEAK_full_contrasts_treat,
   PEAK_full_contrasts_stat, peak_data, MIN_AGE, MAX_AGE, theme_clean, custom_order)

