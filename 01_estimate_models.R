################################################################################
# Code by James Bland
#
# This code estimates OLS, instrumental variable regression, as well as pooled 
# and hierarchical Bayesian models for the belief updating model of Grether (1980) 
#
# Keywords: Bayes rule, conservatism, base-rate neglect
#
# If you use any of this code please cite: 
#
#   - Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: 
#     Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN.
#
#   INPUT: data and stan files
#     interval_pooled.stan
#     interval_hierarchical.stan
#     00_scenarios.csv - scenarios used for decisions (i.e., prior and urn composition)
#     00_simulated_data.csv - choice data by simulated agents
#
#   OUTPUT: Three csv files
#     01_fit_reduced_form.rds - OLS and IV Estimates
#     01_fit_pooled.rds - pooled interval model estimates
#     01_fit_hierarchical.rds - hierarchical model estimates
#
################################################################################


# Make sure the working directory is the one that contains this file
# setwd("C:/Users/jbland/Dropbox/LearningRegressionModel/ExperimentAnalysis/ExampleCode")

# load libraries and set options
library(tidyverse)
library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
library(AER) # needed for IV and clustered standard errors
  
  
# load and compile the Stan files

Pooled<-"interval_pooled.stan" |> stan_model()

Hierarchical<-"interval_hierarchical.stan" |> stan_model()
  
################################################################################
# load the data
  
  Scenarios<- "00_scenarios.csv" |> read.csv()

# replace zeros with this and ones with 1-this (only applies for IV estimation)
replacewith<-0.01

  DATA <- "00_simulated_data.csv" |>
    read.csv() |>
    arrange(id,scenario,decision) |>
    mutate(
      
      # determine upper and lower bounds for interval models
      rounded_to_10 = (belief*10-floor(belief*10)) ==0 ,
      rounded_to_05 = (belief*20-floor(belief*20)) ==0 ,
      rounded_to_01 = (belief*100-floor(belief*100)) ==0 ,
      
      belief_lb = ifelse(rounded_to_10,belief-0.05,
                         ifelse(rounded_to_05,belief-0.025,belief-0.005)
      ),
      belief_ub = ifelse(rounded_to_10,belief+0.05,
                         ifelse(rounded_to_05,belief+0.025,belief+0.005)
      ),
      # deal with belief bounds <0 or >1
      belief_lb = ifelse(belief_lb<0,0,belief_lb),
      belief_ub = ifelse(belief_ub>1,1,belief_ub),
      # convert these to bounds for y
      y_lb = logit(belief_lb),
      y_ub = logit(belief_ub),
      
      # signal information
      lambda = ifelse(is.na(signal),0,ifelse(signal,1,-1)*logit(strength))
    ) |>
    # these bits are only needed if you want to estimate IV or OLS as well
    group_by(id,scenario) |>
    mutate(
      # logit belief
      y = ifelse(belief==0,replacewith,ifelse(belief==1,1-replacewith,belief)) |>
        logit(),
      # lagged logit belief
      ly = dplyr::lag(y),
      # lagged Bayesain belief
      lyBayes = logit(prior)+cumsum(lambda)
    ) |>
    ungroup()
  
  
  #-------------------------------------------------------------------------------
  # Estimate OLS and IV models
  #-------------------------------------------------------------------------------
  
  # OLS model
  OLS<- DATA |> 
    lm(formula = y~ly+lambda+0) 
  OLS.se<-OLS|>
    # cluster standard errors by participant
    coeftest(vcov = vcovCL,cluster = ~ id)
  
  # IV model
  # Here we are using lagged Bayesian belief as the instrument 
  # for lagged reported belief
  IV<- ivreg(formula = y~ly+lambda+0 | lyBayes+lambda,data=DATA)
  IV.se<-IV |>
    # cluster standard errors by participant
    coeftest(vcov = vcovCL,cluster = ~ id)
  
  list(
    OLS = OLS,
    OLS.se = OLS.se,
    IV = IV,
    IV.se = IV.se
  ) |>
    saveRDS("01_fit_reduced_form.rds")

  

#-------------------------------------------------------------------------------
# Estimate Bayesian Models
#
# Note: The Stan programs need the data in matrix form, with each column 
# corresponding to a participant-scenario pair, and each row corresponding to 
# a decision
#-------------------------------------------------------------------------------

y_lb <- DATA |>
  pivot_wider(
    id_cols = decision,
    names_from = c(id,scenario),
    values_from = y_lb
  ) |>
  select(-decision)

y_ub <- DATA |>
  pivot_wider(
    id_cols = decision,
    names_from = c(id,scenario),
    values_from = y_ub
  ) |>
  select(-decision)

lambda <- DATA |>
  pivot_wider(
    id_cols = decision,
    names_from = c(id,scenario),
    values_from = lambda
  ) |>
  select(-decision)

id <- DATA |>
  pivot_wider(
    id_cols = decision,
    names_from = c(id,scenario),
    values_from = id
  ) |>
  select(-decision)

id<-id[1,]

#-------------------------------------------------------------------------------
# pooled estimation (takes about 20 min for 30 agents with 4 decisions 
# on a 2017 4-core i7 processor with Win 11)
#-------------------------------------------------------------------------------

  # get the data into a Stan-friendly format 
  
  dStan<-list(
    nScenarios = dim(lambda)[2],
    nDecisions = dim(lambda)[1],
    lambda = lambda,
    prior = (DATA |> filter(decision==1))$prior,
    y_lb = y_lb,
    y_ub = y_ub,
    
    prior_delta = c(1,1),
    prior_beta = c(1,1),
    prior_sigma_eps = 1,
    prior_sigma_eta = 1
  )
  
  Fit_Pooled<- Pooled |>
    sampling(
      
      data=dStan,seed=42,
      
      # Running this on RStan's default 4 chains leads to a bulk ESS warning.
      # An alternative could be to increase the number of iterations
      chains = 8,       
             
      # These options mean we don't store the augmented data, only the parameters
      pars = "ymeas",include=FALSE
      )
  
  # save the STanfit object
  Fit_Pooled |>
    saveRDS("01_fit_pooled.rds")


#-------------------------------------------------------------------------------
# hierarchical estimation (takes about 20 min for 30 agents with 4 decisions 
# on a 2017 4-core i7 processor with Win 11)
#-------------------------------------------------------------------------------

  dStan<-list(
    nScenarios = dim(lambda)[2],
    nDecisions = dim(lambda)[1],
    nParticipants = DATA$id |> unique() |> length(),
    id = id |> unlist() |> as.vector(),
    
    lambda = lambda,
    prior = (DATA |> filter(decision==1))$prior,
    y_lb = y_lb,
    y_ub = y_ub,
    
    prior_MU = list(
      c(1,0.25),
      c(1,0.25),
      c(0,0.25),
      c(0,0.25)
    ),
    prior_TAU = c(1,1,1,1),
    prior_OMEGA = 4
  )
  
  
  Fit_Hierarchical<- Hierarchical |>
    sampling(
      data=dStan,seed=42,
      
      # Running this on RStan's default 4 chains * 2000 iterations leads to a bulk & tail ESS warning.
      chains = 8,

      # These options mean we don't store the augmented and intermediate data, only the parameters
      pars = c("ymeas","z","L_OMEGA"),include=FALSE
      )
  
  Fit_Hierarchical |>
    saveRDS("01_fit_hierarchical.rds")


