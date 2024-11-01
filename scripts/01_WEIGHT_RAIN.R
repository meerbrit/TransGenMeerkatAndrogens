###### Model weights of 2nd Gen meerkat litters and determine offset to expected weight#####
#### XXX 2024 ####

rm(list = ls(all.names = TRUE)) 
gc() 

# Load required libraries
library(rstan)
library(dplyr)
library(zoo)
library(lubridate)
library(ggplot2)
library(brms)
library(tidyverse)
library(tidybayes)
library(ggokabeito) # colour palette
library(extrafont)# use font_import() on first use
library(bayestestR)

# set seed to duplicate results
set.seed(42)

# half normal prior and weak normal for intercept
priors_halfnormal <- c(set_prior('normal(0,0.5)', class = 'b', lb = 0), set_prior("normal(0,1)", class = "Intercept"))

final_data <- readRDS('data/final_weight_data.rds')

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

#### MODEL ####
WEIGHT_AGE_MOD_FLUT <- brms::brm(formula = scale(Weight) ~ scale(AGE_D) * SEX * scale(Rainfall30D) + (1|ID), 
                                     data = final_data, family = gaussian(link='identity'),
                                     chains = 4, iter = 10000, warmup = 2500, cores = 4, backend = "cmdstanr", 
                                     prior = priors_halfnormal,
                                     control = list(max_treedepth = 25, adapt_delta=0.99999), init=0, 
                                     threads = threading(4),
                                     file ="WEIGHT_AGE_MOD_2NDGEN")
WEIGHT_AGE_MOD_FLUT <- readRDS('models/WEIGHT_AGE_MOD_2NDGEN.rds')

#### Model details ####
summary(WEIGHT_AGE_MOD_FLUT)

plot(WEIGHT_AGE_MOD_FLUT)
pp_check(WEIGHT_AGE_MOD_FLUT, ndraws=100)

# get the rope range  = -0.1 * SDy, 0.1 * SDy
# as its scaled, sd = 1 so -0.1 and 0.1 it is!
ropeRange <- c(-0.1* sd(scale(final_data$Weight)), 0.1 * sd(scale(final_data$Weight)))

describe_posterior(
  WEIGHT_AGE_MOD_FLUT,
  effects = "all", 
  component = "all",
  rope_range = ropeRange,  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(WEIGHT_AGE_MOD_FLUT)

bayes_R2(WEIGHT_AGE_MOD_FLUT)

performance::variance_decomposition(WEIGHT_AGE_MOD_FLUT)


### predict weight for FLUT pups ####
#load predictions 
flut_data <- readRDS("data/flut_weight_data.rds")
weight_pred <- readRDS("data/FLUT_RAIN_weight.rds")

flut_data$WEIGHT_PRED <- as.numeric(weight_pred[,1])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)
flut_data$WEIGHT_CI.lo <- as.numeric(weight_pred[,3])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)
flut_data$WEIGHT_CI.hi <- as.numeric(weight_pred[,4])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)

#plot predictions and actual average weight (jitter)
ggplot(flut_data, aes(x = REC_AGE_D, y = WEIGHT_PRED, color = TREATMENT)) +
  geom_jitter(alpha = 1, aes(shape = "predicted")) +
  geom_smooth(method = 'loess',  alpha = 0.3) +
  geom_jitter(aes(y = AvgWeight, shape = "average")) +
  geom_segment(aes(x = REC_AGE_D, xend = REC_AGE_D, y = AvgWeight, yend = WEIGHT_PRED), linetype = "dashed") +
  scale_fill_okabe_ito() +
  scale_color_okabe_ito(order = c(2, 3, 1), name = "Offspring", labels = c('DC', 'DT', 'SC')) +
  labs(x = "Age in days", y = "Weight") +
  scale_x_continuous(breaks = seq(0, 130, 10)) +
  scale_y_continuous(n.breaks = 10) +
  theme_clean() +
  guides(fill = FALSE, color = guide_legend(title = "Offspring"),
         shape = guide_legend(title = "Weight"))


# expectation of the posterior predictive distribution: epred
prediction130D_data <- WEIGHT_AGE_MOD_FLUT %>% 
  epred_draws(newdata = expand_grid(AGE_D = seq(min(final_data$AGE_D), max(final_data$AGE_D), by=5),
                                    SEX= levels(final_data$SEX),
                                    Rainfall30D = mean(final_data$Rainfall30D)),
              re_formula = NULL, allow_new_levels=T, robust=T)

prediction130D_data$Weight <- prediction130D_data$.epred * sd(final_data$Weight) + mean(final_data$Weight)
prediction130D_data$SEX <- as.factor(prediction130D_data$SEX)

ggplot(prediction130D_data, aes(x = AGE_D, y = Weight)) +  
  stat_lineribbon(.width = .95, alpha = 0.5) +
  geom_segment(data = flut_data, aes(x = REC_AGE_D, xend = REC_AGE_D, y = AvgWeight, yend = WEIGHT_PRED, color = TREATMENT), linetype = "dashed", alpha = 0.7) +
  geom_jitter(data = flut_data, aes(x = REC_AGE_D, y = AvgWeight, colour = TREATMENT, shape = 'average'), alpha = 1, size = 1.5) +
  geom_jitter(data = flut_data, aes(x = REC_AGE_D, y = WEIGHT_PRED, colour = TREATMENT, shape = 'predicted'), alpha = 1, size = 2) +
  scale_color_okabe_ito(order = c(2, 3, 1), name = "Offspring", labels = c('DC', 'DT', 'SC')) +
  scale_fill_okabe_ito(order = c(8), alpha = 0.2, guide = 'none') +
  labs(x = "Age (days)", y = "Weight (g)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_shape_manual(name = "Weight", values = c('average' = 16, 'predicted' = 3)) +
  theme_clean()+
  guides(fill=FALSE)

############## determine WEIGHT OFFSET ####
flut_data$WEIGHT_DIFF <- flut_data$AvgWeight - flut_data$WEIGHT_PRED 

#use percentage
flut_data$WEIGHT_DIFF_PER <-as.numeric(((flut_data$AvgWeight - flut_data$WEIGHT_PRED)/flut_data$WEIGHT_PRED) *100)

ggplot(flut_data, aes(x = REC_AGE_D, y = WEIGHT_DIFF_PER, color = TREATMENT)) +
  geom_point() +
  geom_smooth(method='loess', )+
 # geom_smooth(method='glm' )+
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_color_okabe_ito(order = c(2, 3, 1), name = "Offspring", labels = c('DC', 'DT', 'SC')) +
  labs(x = "Age",
       y = "Weight offset (%)",
       color = "Treatment") +
  theme_clean()


