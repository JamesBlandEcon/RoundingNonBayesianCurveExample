# RoundingNonBayesianCurveExample

This repository contains code in R and Stan for simulating data and then estimating the following versions of the Grether (1980) of belief updating model:
 * OLS
 * IV accounting for classical measurement error
 * Pooled Interval Model accounting for rounding bias Bayesian
 * Hierarchical Interval Model accounting for rounding bias and heterogeneity bias

If you are going to use any of the provided code, please cite:

 * Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN 4996223. http://dx.doi.org/10.2139/ssrn.4996223.


Running the following R scripts in order to produce all of the outputs:

* 00_simulate_data.R simulates data on responses to belief elicitation questions when faced with a two-urn setting common in belief-updating literature
* 01_estimate_models.R estimates OLS, instrumental variable regression, as well as pooled and hierarchical Bayesian models for the belief updating model of Grether (1980)
  * 01_estimate_models.R requires two .stan files to run:
    * interval_pooled.stan for the pooled model
    * interval_hierarchical.stan for the hierarchical model 
* 02_postestimation.R demonstrates some ways to report the estimates




