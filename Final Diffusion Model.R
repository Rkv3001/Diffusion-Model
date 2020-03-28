# Seeting working directory
setwd("C:/Users/h/Desktop/WIP/2")

# Installing necessory libraries
#devtools::install_github("rtdists/rtdists")
install.packages("rtdists")
library("rtdists")
install.packages("dplyr")
library('dplyr')
install.packages("tidyverse")
library("tidyverse")
install.packages("binom")
library("binom")
install.packages("ggplot2")
library("ggplot2")
install.packages("ggpmisc")
library("ggpmisc")


# Read Data
medical <- read.csv("dm.csv")

# Dropping NA values from the data
medical <- medical %>% drop_na()

# Filtering response timebased on conditon
medical <- medical %>%
  filter(rt > 0.025 &  rt < 2.5)

# Creating a response column for upper and lower case where upper = correclty response
# and lower = incorrect response based on the classificcasion and response given
medical$actual <- medical$response
medical$response <- NULL
medical$response <- c(ifelse(medical$classification==medical$actual, "upper", "lower"))
medical$group <- c(ifelse(medical$group == "novice", "novice", "med p"))




#get start values
set.seed(9)

# Function to initialize random values
get_start_values <- function() {
  c(
    a = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5),
    vbe = rnorm(1, 0, 2),
    vbh = rnorm(1, 0, 2),
    vnbe = rnorm(1, 0, 2),
    vnbh= rnorm(1, 0, 2),
    z = runif(1, 0.4, 0.6),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5)
  )
}


#Function to calculate 
ll_diffusion <- function(pars, stimuli)
{
  t1 <- stimuli %>% filter(classification=='blast', difficulty == 'easy')
  t2 <- stimuli %>% filter(classification=='non-blast', difficulty == 'easy')
  t3 <- stimuli %>% filter(classification=='blast', difficulty == 'hard')
  t4 <- stimuli %>% filter(classification=='non-blast', difficulty == 'hard')
  
  
  pdf1=ddiffusion(t1[['rt']],response =t1[["response"]],
                  a=pars["a"],
                  t0=pars["t0"],
                  v=pars["vbe"],
                  z=pars["z"]*pars["a"],
                  sv=pars["sv"],
                  st0=pars["st0"])
  pdf2=ddiffusion(t2[['rt']],response = t2[["response"]],
                  a=pars["a"],
                  t0=pars["t0"],
                  v=pars["vnbe"],
                  z=pars["z"]*pars["a"],
                  sv=pars["sv"],
                  st0=pars["st0"])
  
  pdf3=ddiffusion(t3[['rt']],response = t3[["response"]],
                  a=pars["a"],
                  t0=pars["t0"],
                  v=pars["vbh"],
                  z=pars["z"]*pars["a"],
                  sv=pars["sv"],
                  st0=pars["st0"])
  
  pdf4=ddiffusion(t4[['rt']],response = t4[["response"]],
                  a=pars["a"],
                  t0=pars["t0"],
                  v=pars["vnbh"],
                  z=pars["z"]*pars["a"],
                  sv=pars["sv"],
                  st0=pars["st0"])
  
  
  l=sum(log(pdf1))+sum(log(pdf2))+sum(log(pdf3))+sum(log(pdf4))
  if (l==-Inf)
  {
    return(1e20)
  }
  return(-l)
  
}




# Difusion model
diff_model <- medical %>%
  group_by(id, group) %>%
  nest() %>%
  mutate(fits = map(data, ~rerun(5, nlminb(get_start_values(), ll_diffusion,
                                           stimuli = .,
                                           lower = c(0, 0, -Inf, -Inf,-Inf, -Inf, 0, 0, 0))) %>%
                      map_dfr(~as_tibble(cbind(t(.$par),
                                               logLik = -.$objective,
                                               convergence = .$convergence)))))
# This retunrs the best parameters basis MLE
return_pars <- function(df) {
  nc <- ncol(df)
  mle <- df[which.max(df$logLik), -((nc-1):nc) ]
  which_max <- which(round(max(df$logLik), 3) == round(df$logLik, 3))
  ## exclude actual ML estimate:
  which_max <- which_max[which_max != which.max(df$logLik)]
  ## copy ML estimates
  mle2 <- mle
  ## remove all estimates that are not identifiable:
  mle2[abs(mle - df[which_max[1], -((nc-1):nc) ]) > 0.01] <- NA
  mle2
}

# Saving diffusion model
save(diff_model, file = "Model_fit_diffusion.rda")


View(diff_model %>% mutate(res = map(fits, return_pars)) %>% unnest(res))

#############################################################################################################

load("Model_fit_diffusion.rda")

View(diff_model)

# Pasring nest parameters
pars_fit = diff_model %>% mutate(res = map(fits, return_pars)) %>% unnest(res)

View(pars_fit)

# Creating a single columns for all drift rate cases
pars_fit <- pars_fit %>% gather("drift_rates_case", "v", starts_with("v"))


# dealing missing
pars_fit$a[is.na(pars_fit$a)] = mean(pars_fit$a, na.rm = T)
pars_fit$t0[is.na(pars_fit$t0)] = mean(pars_fit$t0, na.rm = T)
pars_fit$z[is.na(pars_fit$z)] = mean(pars_fit$z, na.rm = T)
pars_fit$sv[is.na(pars_fit$sv)] = mean(pars_fit$sv, na.rm = T)
pars_fit$st0[is.na(pars_fit$st0)] = mean(pars_fit$st0, na.rm = T)
pars_fit$v[is.na(pars_fit$v)] = mean(pars_fit$v, na.rm = T)


# This will add an column of lower and upper predicted proportion based model parameters
pars_fit <- pars_fit  %>% group_by(id, group, drift_rates_case) %>%
  mutate(lower_proportion_pred = pdiffusion(Inf, response="lower", 
                                            a=a, v=v, t0=t0, z=z, sv=sv, st0=st0),
         upper_proportion_pred = pdiffusion(Inf, response="upper", 
                                       a=a, v=v, t0=t0, z=z, sv=sv, st0=st0))

# predicted median upper proportion
pars_fit <- pars_fit  %>% group_by(id, group, drift_rates_case) %>%
  mutate(median_upper_proportion_pred = qdiffusion(c(0.5), response="upper", 
                                                   a=a, v=v, t0=t0, z=z, sv=sv, st0=st0, scale_p = TRUE))

# For upper proportion of the data
pars_fit <- pars_fit %>%
  mutate(vbe = data %>% map(. %>% filter(difficulty == "easy" & classification == 'blast') %>% 
                              summarise(up = mean(response == "upper")) %>% .[["up"]])) %>%
  mutate(vnbe = data %>% map(. %>% filter(difficulty == "easy" & classification == 'non-blast') %>% 
                              summarise(up = mean(response == "upper")) %>% .[["up"]])) %>%
  mutate(vbh = data %>% map(. %>% filter(difficulty == "hard" & classification == 'blast') %>% 
                               summarise(up = mean(response == "upper")) %>% .[["up"]])) %>%
  mutate(vnbh = data %>% map(. %>% filter(difficulty == "hard" & classification == 'non-blast') %>% 
                               summarise(up = mean(response == "upper")) %>% .[["up"]]))
pars_fit <- pars_fit %>% gather("data_case_upper", "upper_proportion_data", c(vbe, vnbe, vbh, vnbh))
pars_fit$condition <- c(pars_fit$drift_rates_case == pars_fit$data_case_upper)
pars_fit <- pars_fit %>%
  filter(condition == TRUE)
pars_fit$condition <- NULL
pars_fit$data_case_upper <- NULL

# for Lower Proportion of the data
pars_fit <- pars_fit %>%
  mutate(lower_proportion_data = 1 - as.numeric(upper_proportion_data))

# Median upper proportion of data
pars_fit <- pars_fit %>%
  mutate(vbe = data %>% map(. %>% filter(difficulty == "easy" & classification == 'blast') %>% 
                              summarise(up = as.numeric(quantile(rt, probs = 0.5))) %>% .[["up"]])) %>%
  mutate(vnbe = data %>% map(. %>% filter(difficulty == "easy" & classification == 'non-blast') %>% 
                               summarise(up = as.numeric(quantile(rt, probs = 0.5))) %>% .[["up"]])) %>%
  mutate(vbh = data %>% map(. %>% filter(difficulty == "hard" & classification == 'blast') %>% 
                              summarise(up = as.numeric(quantile(rt, probs = 0.5))) %>% .[["up"]])) %>%
  mutate(vnbh = data %>% map(. %>% filter(difficulty == "hard" & classification == 'non-blast') %>% 
                               summarise(up = as.numeric(quantile(rt, probs = 0.5))) %>% .[["up"]]))
pars_fit <- pars_fit %>% gather("data_case_upper", "median_upper_proportion_data", c(vbe, vnbe, vbh, vnbh))
pars_fit$condition <- c(pars_fit$drift_rates_case == pars_fit$data_case_upper)
pars_fit <- pars_fit %>%
  filter(condition == TRUE)
pars_fit$condition <- NULL
pars_fit$data_case_upper <- NULL



#############################################################################################################

# ANOVA and plot

pars_fit$upper_proportion_data <- as.numeric(pars_fit$upper_proportion_data)
pars_fit$upper_proportion_pred <- as.numeric(pars_fit$upper_proportion_pred)
pars_fit$median_upper_proportion_pred <- as.numeric(pars_fit$median_upper_proportion_pred)
pars_fit$median_upper_proportion_data <- as.numeric(pars_fit$median_upper_proportion_data)


blast <- pars_fit %>% filter(drift_rates_case=='vbe'||drift_rates_case == 'vbh')
non_blast <- pars_fit %>% filter(drift_rates_case=='vnbe'||drift_rates_case == 'vnbh')

# Plot
#install.packages("ggpmisc")
library("ggpmisc")
my.formula <- y ~ x
# for each group all ids pred and actul proportion plot for blast
# Overall
ggplot(pars_fit, aes(x = upper_proportion_data,y = upper_proportion_pred)) +
  geom_smooth(method = "lm",formula = my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste( ..rr.label..)), 
               parse = TRUE) +
  geom_point()+ ggtitle("Both (Blast & non-blast) - Upper Proportion (actual vs predicted)") +
  facet_wrap(~group)


# for Blast
ggplot(blast, aes(x = upper_proportion_data,y = upper_proportion_pred)) +
  geom_smooth(method = "lm",formula = my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste( ..rr.label..)), 
               parse = TRUE) +
  geom_point()+ ggtitle("Blast - Upper Proportion (actual vs predicted)") +
  facet_wrap(~group)

# For non-blast
ggplot(non_blast, aes(x = upper_proportion_data,y = upper_proportion_pred)) +
  geom_smooth(method = "lm",formula = my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste( ..rr.label..)), 
               parse = TRUE) +
  geom_point()+ ggtitle("Non-Blast - Upper Proportion (actual vs predicted)") +
  facet_wrap(~group)

# For median proportion overall
ggplot(pars_fit, aes(x = median_upper_proportion_data,y = median_upper_proportion_pred)) +
  geom_smooth(method = "lm",formula = my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste( ..rr.label..)), 
               parse = TRUE) +
  geom_point()+ ggtitle("Both (Blast & non-blast) - Median Proportion (actual vs predicted)") +
  facet_wrap(~group)

# For median proportion blast
ggplot(blast, aes(x = median_upper_proportion_data,y = median_upper_proportion_pred)) +
  geom_smooth(method = "lm",formula = my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste( ..rr.label..)), 
               parse = TRUE) +
  geom_point()+ ggtitle("Both (Blast & non-blast) - Median Proportion (actual vs predicted)") +
  facet_wrap(~group)
# For median proportion non-blast
ggplot(non_blast, aes(x = median_upper_proportion_data,y = median_upper_proportion_pred)) +
  geom_smooth(method = "lm",formula = my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste( ..rr.label..)), 
               parse = TRUE) +
  geom_point()+ ggtitle("Both (Blast & non-blast) - Median Proportion (actual vs predicted)") +
  facet_wrap(~group)


##############################################################################################################

# Pasring nest parameters
load("Model_fit_diffusion.rda")

data = diff_model %>% mutate(res = map(fits, return_pars)) %>% unnest(res)

View(data)

# For comparing vbe and vnbe
summary(aov(vbe~vnbe, data = data))
# There is no difference

# For comparing vbe and vbh
summary(aov(vbe~vbh, data = data))
# There is no difference

# For comparing vbe and vnbh
summary(aov(vbe~vnbh, data = data))
# These are signiicantly different

# For comparing vnbe and vbh
summary(aov(vnbe~vbh, data = data))
# No difference

# For comparing vnbe and vnbh
summary(aov(vnbe~vnbh, data = data))
# There is no difference

# For comparing vbh and vnbe
summary(aov(vbh~vnbh, data = data))
# These are signiicantly different
