/*------------------------------------------------------------------------------
 Code by James Bland
 
 This Stan program implements the a hierarchical Bayesian interval regression model
 for the belief updating model of Grether (1980) 
 
 Keywords: Bayes rule, conservatism, base-rate neglect
 
 if you use any of this code please cite: 

   - Bland, J. R., & Rosokha, Y. (2024). Rounding the (Non) Bayesian Curve: 
     Unraveling the Effects of Rounding Errors in Belief Updating. Available at SSRN.
------------------------------------------------------------------------------*/
data {
  int<lower=0> nScenarios; // number of scenarios
  int<lower=0> nDecisions; // number of decisions
  int<lower=0> nParticipants; // number of participants
  int id[nScenarios]; // participant unique id
  
  matrix[nDecisions,nScenarios] lambda; // signal information
  row_vector[nScenarios] prior; // prior for red urn
  matrix[nDecisions,nScenarios] y_lb; // lower bound for y
  matrix[nDecisions,nScenarios] y_ub; // upper bound for y
  
  /* prior for the means of individual-level parameters
     In order, the individual-level parameters are:
        delta (base-rate neglect)
        beta (conservatism)
        sigma_eps (additive "implementation" error standard deviation)
        sigma_eta (classical measurement error standard deviation)
    The prior is specified as {(mean1,sd1),(mean2,sd2),(mean3,sd3),(mean4,sd4)}
    where these means and standard deviations are for Normal distributions
  */
  vector[2] prior_MU[4]; 
  /* prior for the standard deviations of individual-level parameters
     The prior for each standard deviation parameter takes on a half-Cauchy
     (truncated to positive numbers) distribution with central tendancy parameter zero
  */
  vector[4] prior_TAU; 
  /* prior for the correlation matrix between individual-level parameters
     The prior for the correlation matrix takes on a LKJ distribution, which 
     means that the prior density of the correlation matrix is proportional 
     to its determinant raised to the power of (prior_OMEGA-1)
  */
  real prior_OMEGA; 
}

transformed data {
  // pre-compute the logit of the prior
  row_vector[nScenarios] logit_prior = logit(prior);
  
}

parameters {
  
  // mean of (delta_i,beta_i,sigma_eps_i,sigma_eta_i)
  vector[4] MU;
  // sd of (delta_i,beta_i,sigma_eps_i,sigma_eta_i)
  vector<lower=0>[4] TAU;
  // Cholesky factor of the correlation matrix of (delta_i,beta_i,sigma_eps_i,sigma_eta_i)
  cholesky_factor_corr[4] L_OMEGA;
  
  // normalized individual-level parameters
  matrix[4,nParticipants] z;
  
  // augmented logit belief with measurement error
  matrix<lower=y_lb,upper=y_ub>[nDecisions,nScenarios] ymeas;
}

transformed parameters {
  
  // individual-level parameters
  vector[nParticipants] delta; // base-rate neglect
  vector[nParticipants] beta; // conservatism
  vector[nParticipants] sigma_eps; // additive error sd
  vector[nParticipants] sigma_eta; // classical measurement error sd
  
  {
    // transform normalized individual-level parameters (z) into the parameters
    // we actually want to estimate
    matrix[4,nParticipants] theta = rep_matrix(MU,nParticipants)+diag_pre_multiply(TAU,L_OMEGA)*z;
    
    delta = theta[1,]';
    beta = theta[2,]';
    sigma_eps = exp(theta[3,]');
    sigma_eta = exp(theta[4,]');
  }
  
}

model {
  
  // hierarchical structure
  
  target += std_normal_lpdf(to_vector(z));
  
  // priors
  for (pp in 1:4) {
    target += normal_lpdf(MU[pp] | prior_MU[pp][1],prior_MU[pp][2]);
    target += cauchy_lpdf(TAU[pp] | 0.0,prior_TAU[pp]);
  }
  target += lkj_corr_cholesky_lpdf(L_OMEGA | prior_OMEGA);
  
  
  /* construct variance-covariance matrix
   See Appendix A.1 of the main text for more details
   Note that this is constant for each participant, so we only need to compute
   it once per participant, not once per scenario.
  */
  matrix[nDecisions,nDecisions] SIGMA[nParticipants];
  
  for (ii in 1:nParticipants) {
  SIGMA[ii] = diag_matrix(rep_vector(pow(sigma_eta[ii],2.0),nDecisions));
  
  for (tt in 1:nDecisions) {
    for (ss in 1:nDecisions) {
      
      for (pp in 0:(tt-1)) {
        for (qq in 0:(ss-1)) {
          
          if ((tt-pp)==(ss-qq)) {
            SIGMA[ii][tt,ss]+=pow(delta[ii],pp+qq)*pow(sigma_eps[ii],2.0);
          }
        }
      }
      
    }
  }
  }
  for (ss in 1:nScenarios) {
   vector[nDecisions] m = rep_vector(0.0,nDecisions);  
   matrix[nDecisions,nDecisions] s = SIGMA[id[ss]];
   real d = delta[id[ss]];
   real b = beta[id[ss]];
   
   // here we calculate the mean component recursively
   real lb = logit_prior[ss];
   for (dd in 1:nDecisions) {
     
     lb = d*lb+b*lambda[dd,ss];
     
     m[dd] = lb;
     
   }
   // augment the data
   target += multi_normal_lpdf(ymeas[,ss] | m, s);
    
    
    
  }
  
  

}

generated quantities {
  /* Convert the Cholesky-factored correlation matrix back into the actual
  correlation matrix
  */
  matrix[4,4] OMEGA = L_OMEGA*L_OMEGA';
  
}

