data{
  int<lower=1> N; //number of trees
  int<lower=1> P; //number of plots
  int<lower=1> Plot[N]; //plot index per tree
  int TotNewStems[N];//the response variable, post-fire stem recruitment
  
}
parameters{
  vector[P] eta_a_plot; //non-centered parameterization to improve sampling
  real a; //global intercept
  real<lower=0> sigma_plot;  //variance
  real phi_reciprocal;
}
transformed parameters {
  vector<lower=0>[N] mu;  
  vector[P] vary_a;   //varying intercept component
  real<lower = 0> phi;
  vary_a = sigma_plot*eta_a_plot;
  for ( i in 1:N ) {
    mu[i] = exp(a + vary_a[Plot[i]]);
  }
  phi = 1. / phi_reciprocal;
}
model{
  a ~ normal( 0 , 10 );
  eta_a_plot ~ normal( 0 , 1 );
  sigma_plot ~ cauchy( 0 , 1 );
  phi_reciprocal ~ cauchy(0, 5);
  TotNewStems ~ neg_binomial_2(mu, phi); //final likelihood part of model
}
generated quantities{
  vector[N] y_rep;
  vector[N] log_lik;
  for ( i in 1:N ) {
    y_rep[i] = neg_binomial_2_rng(mu[i], phi);
    log_lik[i] = neg_binomial_2_lpmf(TotNewStems[i] | mu[i], phi);
  }
}