/*------------------------------------------------------------------------------
 Code by James Bland
 
 This Stan program implements the a pooled Bayesian interval regression model
 for the belief updating model of Grether (1980) 
 
 Keywords: Bayes rule, conservatism, base-rate neglect
 
 if you use any of this code please cite: 

   - Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: 
     Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN.
------------------------------------------------------------------------------*/
data {
  int<lower=0> nScenarios; // number of scenarios
  int<lower=0> nDecisions; // number of decisions
  matrix[nDecisions,nScenarios] lambda; // signal information
  row_vector[nScenarios] prior; // prior for red urn
  matrix[nDecisions,nScenarios] y_lb; // lower bound for y
  matrix[nDecisions,nScenarios] y_ub; // upper bound for y
  
  vector[2] prior_delta; // normal
  vector[2] prior_beta; // normal
  real prior_sigma_eps;// half Cauchy
  real prior_sigma_eta;// half Cauchy
}

transformed data {
  
  // pre-compute the logit of the prior
  row_vector[nScenarios] logit_prior = logit(prior);
  
}

parameters {
  
  // base rate neglect
  real delta;
  // conservatism
  real beta;
  // additive "implementation" error standard deviation
  real<lower=0> sigma_eps;
  // classical measurement error standard deviation
  real<lower=0> sigma_eta;
  
  // augmented logit belief with measurement error
  matrix<lower=y_lb,upper=y_ub>[nDecisions,nScenarios] ymeas;
}

model {
  
  
  // priors
  
  target += normal_lpdf(delta| prior_delta[1],prior_delta[2]);
  target += normal_lpdf(beta| prior_beta[1], prior_beta[2]);
  target += cauchy_lpdf(sigma_eps| 0.0,prior_sigma_eps);
  target += cauchy_lpdf(sigma_eta| 0.0,prior_sigma_eta);
  
  // construct variance-covariance matrix
  // See Appendix A.1 of the main text for more details
  
  matrix[nDecisions,nDecisions] SIGMA = rep_matrix(0.0,nDecisions,nDecisions);
  
  for (tt in 1:nDecisions) {
    for (ss in 1:nDecisions) {
      
      if (tt==ss) {
        SIGMA[tt,ss] += pow(sigma_eta,2.0);
      }
      
      for (pp in 0:(tt-1)) {
        for (qq in 0:(ss-1)) {
          
          if ((tt-pp)==(ss-qq)) {
            SIGMA[tt,ss]+=pow(delta,pp+qq)*pow(sigma_eps,2.0);
          }
        }
      }
      
    }
  }
  
  
   matrix[nDecisions,nScenarios] MU;  
   
   // here I calculate the mean component recursively
   row_vector[nScenarios] lb = logit_prior;
   for (dd in 1:nDecisions) {
     
     lb = delta*lb+beta*lambda[dd,];
     
     MU[dd,] = lb;
     
   }
  
  
  
  for (ss in 1:nScenarios) {
    
  
   // augment the data by logit belief with measurment error
   target += multi_normal_lpdf(ymeas[,ss] | MU[,ss], SIGMA);
    
    
    
  }
  
  

}

