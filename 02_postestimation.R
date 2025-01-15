################################################################################
# Code by James Bland
#
# If you use any of this code please cite: 
#
#   - Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: 
#     Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN.
#
# This script demonstrates some ways to report the estimates
# 
#   INPUT: data and stan files
#     interval_pooled.stan
#     interval_hierarchical.stan
#     00_simulated_scenarios.csv - scenarios used for decisions (i.e., prior and urn composition)
#     00_simulated_data.csv - choice data by simulated agents
#
#   OUTPUT: 
#     02_reduced_form_estimates.txt, a regression table showing the OLS and IV estimates
#     02_posterior_moments_summary.html, a summary of some of the posterior moments of both Bayesian models
#     02_hierarchical_population_moments.html, a detailed summary of all of the population-level parameters of the hierarchical model
#     02_hierarchical_individual_estimates.png, a plot of the individual-level "shrinkage" estimates from the hierarchical model
#     02_hierarchical_comparison_to_truth.png, a comparison between the true (simulated) individual-level parameters and their "shrinkage" estimates from the hierarchical model
#
################################################################################

# Make sure the working directory is the one that contains this file
# setwd("C:/Users/jbland/Dropbox/LearningRegressionModel/ExperimentAnalysis/ExampleCode")

# load libraries and set options
library(tidyverse)
library(rstan)
library(kableExtra)
library(stargazer)
library(latex2exp)

################################################################################
# load the data from estimation

ReducedForm<-"01_fit_reduced_form.rds" |> readRDS()
IntervalPooled<-"01_fit_pooled.rds" |>  readRDS()
IntervalHierarchical<-"01_fit_hierarchical.rds" |>  readRDS()

################################################################################
# OLS and IV results summarized in a regression table

stargazer(ReducedForm$OLS,ReducedForm$IV,
          se = list(ReducedForm$OLS.se[,"Std. Error"],
                    ReducedForm$IV.se[,"Std. Error"]
                    ),
          covariate.labels = c("delta","beta")
          ,type="text",
          out = "02_reduced_form_estimates.txt")


################################################################################
# Check the convergence diagnostics of the posterior simulations
#-------------------------------------------------------------------------------

IntervalPooled |>
  check_hmc_diagnostics() |>
  print()

IntervalHierarchical |>
  check_hmc_diagnostics() |>
  print()


#-------------------------------------------------------------------------------
# Getting the Stanfit objects into summary form
# If we are only interested in reporting posterior mean, sd, and percentiles,
# then this is a useful format to have the fitted models in
#-------------------------------------------------------------------------------


IntervalPooledSummary<- summary(IntervalPooled)$summary |>
  data.frame() |>
  rownames_to_column(var = "par") 

IntervalHierarchicalSummary<-summary(IntervalHierarchical)$summary |>
  data.frame() |>
  rownames_to_column(var = "par") 

# For example, we can look at the summary of the posterior distribution from the 
# pooled model like this:

IntervalPooledSummary |>
  print()

# Here the columns are:
#   par: parameter name
#   mean: posterior mean
#   se_mean: standard error of the posterior mean
#         (i.e. the simulation accuracy for the mean)
#   sd: posterior standard deviation
#   X2.5. ... X97.5.: percentiles of the posterior distribution
#   n_eff: effective sample size of simulation
#   Rhat: Gelman-Rubin statistic, a diagnostic for convergence

# Now looking at the hierarchical model summary

IntervalHierarchicalSummary |>
  print()

# Note that we now have indices on the parameters, for example MU[1] corresponds
# to the mean of delta, MU[2] corresponds to the mean of beta, and delta[1]
# corresponds to delta of participant id=1.
# Here we can keep track of these indices as follows:

IntervalHierarchicalSummary<-IntervalHierarchicalSummary |>
  mutate(
    parind = parse_number(par),
    parameter = gsub("[^A-Za-z]","", par)
  )

#-------------------------------------------------------------------------------
# Tabulating the population-level estimates of delta and beta
# Similar to Table 5 of the main text
#-------------------------------------------------------------------------------

# format to display numbers
# %.3f reports to 3 decimal places
fmt<- "%.3f"

dPooled<-IntervalPooledSummary |>
  mutate(
    model = "Pooled",
    mean = sprintf(fmt,mean),
    sd = paste0("(",sprintf(fmt,sd),")")
  ) |>
  filter(par =="delta" | par=="beta") |>
  select(model,par,mean,sd) |>
  pivot_longer(
    cols = c(mean,sd),
    names_to = "quantity",
    values_to = "value"
  )

dHierarchical<-IntervalHierarchicalSummary |>
  # keep only the mean vector
  filter(grepl("MU",par)) |>
  mutate(
    model = "Hierarchical",
    mean = sprintf(fmt,mean),
    sd = paste0("(",sprintf(fmt,sd),")"),
    par = ifelse(par=="MU[1]","delta",ifelse(par=="MU[2]","beta",par))
  ) |>
  filter(par =="delta" | par=="beta") |>
  select(model,par,mean,sd) |>
  pivot_longer(
    cols = c(mean,sd),
    names_to = "quantity",
    values_to = "value"
  )

TAB<-rbind(dPooled,dHierarchical) |>
  pivot_wider(
    id_cols = c(par,quantity),
    names_from = model,
    values_from = value
  ) |>
  mutate(
    par = ifelse(quantity=="sd","",par)
  ) |>
  select(-quantity) |>
  rename(parameter=par) |>
  kbl("html"
    ) |>
  kable_classic(full_width=FALSE) 
TAB |> print()

TAB |>
  readr::write_file("02_posterior_moments_summary.html")

#-------------------------------------------------------------------------------
# Hierarchical model summary of population-level parameters
# Similar to Table C2 of the main text
#-------------------------------------------------------------------------------

parList<-c("$\\delta$","$\\beta$","$\\log\\sigma_\\epsilon$","$\\log\\sigma_\\eta$")
fmt<-"%.3f"

Population<-IntervalHierarchicalSummary |>
  filter(grepl("MU",par) | grepl("TAU",par) | grepl("OMEGA",par)) |>
  mutate(
    parTeX = parList[parind],
    
    # mean and standard deviation in one string
    msd = paste0(sprintf(fmt,mean)," (",sprintf(fmt,sd),")")
  )

omega<-Population |> 
  filter(grepl("OMEGA",par)) |>
  mutate(
    msd = paste0(sprintf(fmt,mean)," (",sprintf(fmt,sd),")"),
    row = floor(parind/10),
    col = parind-floor(parind/10)*10,
    parTeX = parList[col],
    parameter = parList[row],
    msd = ifelse(msd == "1.000 (0.000)","1",msd)
  ) |>
  pivot_wider(
    id_cols = parameter,
    names_from = parTeX,
    values_from = msd
  )

omegaHeader<-cbind("$\\Omega$","","","","")
colnames(omegaHeader)<-names(omega)




TAB<-Population |>
  filter(parameter!="OMEGA")|>
  pivot_wider(
    id_cols = parameter,
    names_from = parTeX,
    values_from = msd
  ) |>
  mutate(parameter = ifelse(parameter=="MU","$\\mu$","$\\tau$")) |>
  rbind(omegaHeader) |>
  rbind(omega) |>
  kbl(caption = "Population-level parameter estimates. Posterior means with posterior standard deviations in parentheses") |>
  row_spec(2, extra_css = "border-bottom: 1px solid") |>
  kable_classic(full_width=FALSE) 

TAB |> print()

TAB |> readr::write_file("02_hierarchical_population_moments.html")

#-------------------------------------------------------------------------------
# Looking at the individual-level estimates from the hierarchical model
# Similar to Figure C3 of the main text
#-------------------------------------------------------------------------------

appender <- function(string) {
  TeX(string) 
}

d<-IntervalHierarchicalSummary |>
  filter(
    parameter =="delta" | parameter == "beta" | parameter =="sigmaeps" | parameter == "sigmaeta"
  ) |>
  group_by(parameter) |>
  mutate(
    parameter = ifelse(
      parameter=="sigmaeps","$\\sigma_\\epsilon$",
      ifelse(parameter=="sigmaeta","$\\sigma_\\eta$",parameter)
    ),
    paramfactor = factor(parameter, levels=c("delta","beta","$\\sigma_\\epsilon$","$\\sigma_\\eta$")),
    cumul = rank(mean)/n()-0.5/n()
  )

plt<-(
  ggplot(d,aes(x=mean,y=cumul))
  +stat_ecdf() 
  +geom_errorbar(aes(xmin = X2.5.,xmax = X97.5.),alpha=0.5)
  +facet_wrap(~paramfactor,scales="free",
              labeller = as_labeller(appender, 
                                     default = label_parsed))
  +theme_bw()
  +ylab("empirical cumulative density of posterior mean")
  +xlab("estimate")
  +labs(caption = "Error bars show 95% Bayesian credible regions (2.5th-97.5th percentile)")
) 
plt |> print()
plt |>
  ggsave(filename="02_hierarchical_individual_estimates.png")

#-------------------------------------------------------------------------------
# Since we have simulated data, we can compare these estimates to the true 
# values that we used to simulate the data
#-------------------------------------------------------------------------------

comparison<-"00_agent_parameters.csv" |>
  read.csv() |>
  # make sure the variable names match the Stan summary
  rename(
    `$\\sigma_\\epsilon$` = sigmaEps,
    `$\\sigma_\\eta$` = sigmaEta,
    parind = id
  ) |>
  select(-X) |>
  pivot_longer(
    cols = c(delta,beta,`$\\sigma_\\epsilon$`,`$\\sigma_\\eta$`),
    names_to = "parameter",
    values_to = "truth"
  ) |>
  # merge this with the Stan summary  of the individual-level parameters
  full_join(
    d,
    by = c("parind","parameter")
  )

plt <- (
  ggplot(comparison,aes(x=truth,y=mean))
  +geom_point(size=0.7)
  +geom_errorbar(aes(ymin = X2.5.,ymax = X97.5.),alpha=0.5)
  +facet_wrap(~paramfactor,scales="free",
             labeller = as_labeller(appender, 
                                    default = label_parsed))
  +theme_bw()
  +xlab("true parameter value")
  +ylab("estimate")
  # 45 degree line
  +geom_abline(slope=1,intercept=0,linetype="dashed")
  +labs(caption = "Vertical coordinates: dots show posterior means, \n error bars show 95% Bayesian credible regions (2.5th-97.5th percentile)")
)
plt |> print()
plt |> ggsave(filename="02_hierarchical_comparison_to_truth.png")
