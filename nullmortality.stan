data{
  int<lower=1> N; //number of trees
  int<lower=1> P; //number of plots
  int<lower=1> Plot[N]; //plot index per tree
  int TotKilled[N]; //the response variable, belowground mortality
}
parameters{
  real a;
  vector[P] eta_a_plot; //non-centered parameterization to improve sampling
  real<lower=0> sigma_plot;  //variance
}
transformed parameters {
  vector[N] p;  //probability of belowground mort
  vector[P] vary_a;   //varying intercept component for plot
  vary_a = sigma_plot*eta_a_plot;
  for ( i in 1:N ) {
    p[i] = a + vary_a[Plot[i]];
  }   // individual scale part of the model, moved as part of N-C parameterization.
}
model{
  a ~ normal( 0 , 5 );
  eta_a_plot ~ normal( 0 , 1 );
  sigma_plot ~ normal( 0 , 1 );
  TotKilled ~ binomial_logit( 1 , p ); //likelihood
}
generated quantities{
  vector[N] pq; 
  vector[N] log_lik; //for calculating loo
  for ( i in 1:N ) {
    pq[i] = a + vary_a[Plot[i]];
    log_lik[i] = binomial_logit_log(TotKilled[i] , 1 , pq[i] );
  }
}