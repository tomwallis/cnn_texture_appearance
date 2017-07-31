### fit the main model for experiment 12 (spatial texture comparison with powerspec model).

# imports -----------------------------------------------------------------
library(brms)
library(loo)

# load_data ---------------------------------------------------------------
# assuming you're running the script from /code/analysis/:
setwd("../..")
top_dir <- getwd()
# top_dir <- "~/Dropbox/Projects/dnn_texture_appearance"
experiment_num <- 12

# here we load the main data file:
results_dir <- file.path(top_dir, "results", paste0("experiment_", experiment_num))
load(file.path(results_dir, "final_data.RData"))

# setup model variables ---------------------------------------------------

chains <- 6
cores <- 6
seed <- 11242375
iter <- 30000 
warmup <- 10000
thin <- 6

# Full model -----------------------------------------------------------

model_formula <- correct ~ image_model * presentation_cond + 
  (image_model * presentation_cond | subj) + 
  (image_model * presentation_cond | image_code)

# get_prior(model_formula, data = dat,
#           family = bernoulli("logit"))

priors <- c(set_prior("normal(0, 2)", class = "b"), 
            set_prior("cauchy(0, 1)", class = "sd"), 
            set_prior("lkj(2)", class = "cor"))

# make_stancode(model_formula, data = dat,
#               family = bernoulli("logit"),
#               prior = priors)

fit <- brm(model_formula, data = dat,
           family = bernoulli("logit"),
           prior = priors, 
           iter = iter,
           warmup = warmup,
           chains = chains,
           thin = thin,
           cores = cores,
           control = list(adapt_delta = 0.99),           
           save_model = file.path(results_dir, "model_file.stan"))

# do the leave-one-out approximation:
loo <- LOO(fit)

save(fit, loo,
     file = file.path(results_dir, "model_fit.RData"))
