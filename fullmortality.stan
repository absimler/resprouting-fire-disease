data{
  int<lower=1> N; //number of trees
  int<lower=1> P; //number of plots
  int<lower=1> Nsims; //number of simulations
  int<lower=1> Plot[N]; //plot index per tree
  int TotKilled[N]; //the response variable, belowground mortality
  real hostCWD[P]; //plot level predictor - host CWD volume
  real standingA[P]; //plot level standing dead host BA for stage A (recent dead w/fine fuels)
  real standingB[P]; //plot level standing dead host BA for stage B (w/o fine fuels)
  real climate[P]; //plot-level annual mean temperature
  real LiveBAc[N]; //individual level predictor, prefire tree size
  real cwdseq[Nsims]; //simulation data for plots
  real Aseq[Nsims]; //simulation data
  real Bseq[Nsims]; //simulation data
  real dbhseq[Nsims]; //simulation data
}
parameters{
  vector[P] eta_a_plot; //non-centered parameterization to improve sampling
  real a; //global intercept
  real bdbh; //parameter for tree size;
  real bcwd; //parameter for CWD fuels volume;
  real ba;   // parameter for standing fuels A;
  real bb; //parameter for standing fuels B;
  real bclim; //parameter for climate/temperature
  real<lower=0> sigma_plot;  //variance
}
transformed parameters {
  vector[N] p;  //probability of belowground mort
  vector[P] vary_a;   //varying intercept component for plot
  vector[P] pred_a;   // average intercept + plot level effects component
  for (j in 1:P){
    pred_a[j] = bcwd*hostCWD[j] + ba*standingA[j] + bb*standingB[j] + bclim*climate[j];
  }   // plot-level part of the model, moved as part of non-centered parameterization.
  vary_a = pred_a + sigma_plot*eta_a_plot;
  for ( i in 1:N ) {
    p[i] = a + vary_a[Plot[i]] + bdbh * LiveBAc[i];
  }   // individual scale part of the model, moved as part of N-C parameterization.
}
model{
  bdbh ~ normal( 0 , 5 ); //priors
  bcwd ~ normal( 0 , 5 );
  ba ~ normal( 0 , 5 );
  bb ~ normal( 0 , 5 );
  bclim~ normal( 0 , 5 );
  a ~ normal( 0 , 5 );
  eta_a_plot ~ normal( 0 , 1 );
  sigma_plot ~ normal( 0 , 1 );
  TotKilled ~ binomial_logit( 1 , p ); //likelihood
}
generated quantities{
  vector[Nsims] p_sims; //simulated probability for plots
  vector[Nsims] y_sims; //simulated outcomes for plotting
  vector[N] pq; 
  vector[N] log_lik; //for calculating loo
  for ( i in 1:N ) {
    pq[i] = a + vary_a[Plot[i]] + bdbh*LiveBAc[i];
    log_lik[i] = binomial_logit_log(TotKilled[i] , 1 , pq[i] );
  }
  for (q in 1:Nsims){
    p_sims[q] = inv_logit(a + bcwd*cwdseq[q] + ba*Aseq[q] + bdbh*dbhseq[q] + bb*Bseq[q]);
    y_sims[q] = bernoulli_rng(p_sims[q]);
  }
}
