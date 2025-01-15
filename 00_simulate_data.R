################################################################################
# Code by James Bland
#
# This code simulates data on responses to belief elicitation questions when faced
# with a two-urn setting common in belief updating literature.
#
# If you use any of this code please cite: 
#
#   - Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: 
#     Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN.
# 
#   INPUT: None
#   OUTPUT: Three csv files
#     00_simulated_parameters.csv - belief updating parameters for each agent
#     00_simulated_scenarios.csv - scenarios used for decisions (i.e., prior and urn composition)
#     00_simulated_data.csv - choice data by simulated agents
################################################################################

# Make sure the working directory is the one that contains this file
# setwd("C:/Users/jbland/Dropbox/LearningRegressionModel/ExperimentAnalysis/ExampleCode")


# load libraries
library(tidyverse)


# the random number seed ensures the same draws each time
set.seed(42)  

#-------------------------------------------------------------------------------
# Set the simulation details
#-------------------------------------------------------------------------------

# rounding behavior. 
# 1 = round to 1 decimal place, 2 = round to 2 decimal places, etc.
rounding <-1

# number of participants
nparticipants<-30
# number of scenarios
nscenarios<-12
# number of decisions per scenario
numdecisions<-4

# consider 2-urn setup up with symmetric compositions
# urn prior (i.e., how likely that urn 1 is selected )
# prior<-runif(nscenarios) #most general case
prior <- sample(seq(0.1, 0.9, by = 0.1), nscenarios, replace = TRUE) #round priors e.g., 50/50

# urn signal strength (i.e., composition of the urn)
# strength<-runif(nscenarios,min=0.5,max = 1) #most general case
strength <- sample(seq(0.6, .9, by = 0.1), nscenarios, replace = TRUE) #suppose there are 10 balls in each urn, so the strength can vary from .6 to .9

# mean and standard deviation
deltaDist <- c(1,0.25) 
betaDist <- c(1,0.25)
# mean and standard deviation of natural logarithm of standard deviation parameters
sigmaEpsDist<-c(log(0.3),0.25)
sigmaEtaDist<-c(log(0.3),0.25)

#-------------------------------------------------------------------------------
# Draw participant-specific parameters
#-------------------------------------------------------------------------------

params<-tibble(
  delta = deltaDist[1]+deltaDist[2]*rnorm(nparticipants),
  beta  = betaDist[1]+betaDist[2]*rnorm(nparticipants),
  sigmaEps = exp(sigmaEpsDist[1]+sigmaEpsDist[2]*rnorm(nparticipants)),
  sigmaEta = exp(sigmaEtaDist[1]+sigmaEtaDist[2]*rnorm(nparticipants))
) |>
  mutate(
    # force the sample means of delta and beta to be exactly one
    delta = delta-mean(delta) + 1,
    beta = beta-mean(beta)+1,
    id = 1:n()
    )

# Save the parameters so that we can compare them
params |>  write.csv("00_agent_parameters.csv")

#-------------------------------------------------------------------------------
# Urn structures
#-------------------------------------------------------------------------------

UrnStructure<-tibble(prior = prior,strength = strength) |>
  mutate(scenario = 1:n())


UrnStructure |>  write.csv("00_scenarios.csv")

#-------------------------------------------------------------------------------
# Simulate data
#-------------------------------------------------------------------------------


# functions we will need to use
logit <- function(x) log(x)-log(1-x)
inv.logit <-function(x) 1/(1+exp(-x))

DATA<-tibble()

d<-expand_grid(params,scenario = 1:nscenarios) |>
  full_join(UrnStructure,by = "scenario") |>
  mutate(
    truth = runif(n())<prior,
    
    # first decision, prior is elicited
    # logit beliefs before classical or rounding measurement error
    ystar = delta*logit(prior)+sigmaEps*rnorm(n()),
    
    decision = 1, # COMMENT: Why is this 1???
    signal = NA # no signal for the first decision
  )

DATA<-rbind(DATA,d)

for (dd in 2:numdecisions) {
  
  d<- d|>
    mutate(
      signal = ifelse(truth,strength,1-strength)<runif(n()),
      
      # logit beliefs before classical or rounding measurement error 
      ystar = delta*ystar+beta*ifelse(signal,1,-1)*logit(strength)+sigmaEps*rnorm(n()),
      decision = dd
    )
  DATA<-rbind(DATA,d)
}

DATA<-DATA |>
  mutate(
    # add in classical measurement error
    ymeas = ystar+sigmaEta*rnorm(n()),
    
    # add in rounding error
    belief = round(inv.logit(ymeas),rounding)
  ) |>
  select(
    id,scenario,decision,prior,strength,truth,signal,belief
  )

# Save the data for estimation
DATA |> write.csv("00_simulated_data.csv")
