# RoundingNonBayesianCurveExample

Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN
This repository contains some example code in R and Stan for simulating data then estimating the models described in the above citation. 

Running the following R scripts in order produces all of the outputs:

* 00_simulate_data.R simulates data on responses to belief elicitation questions when faced with a two-urn setting common in belief updating literature
* 01_estimate_models.R estimates OLS, instrumental variable regression, as well as pooled and hierarchical Bayesian models for the belief updating model of Grether (1980)
* 02_postestimation.R demonstrates some ways to report the estimates

01_estimate_models.R requires two .stan files to run:

* interval_pooled.stan for the pooled model
* interval_hierarchical.stan for the hierarchical model


