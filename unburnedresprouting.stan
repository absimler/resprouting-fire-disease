data{
    int<lower=1> N; //number of trees
    int<lower=1> P; //number of plots
    int<lower=1> Plot[N]; //plot index per tree
    int TotNewStems[N];//the response variable, post-fire stem recruitment
    real LiveBAc[N]; //prefire tree size
    real TopkillBA[N]; //extent of individual damage by fire (sum dead basal area)
    real SODdead[P]; //plot-level pre-fire mortality of SOD hosts
    real climate[P]; //plot-level mean annual temp, 5 years post fire
}
parameters{
    vector[P] eta_a_plot; //non-centered parameterization to improve sampling
    real a; //global intercept
    real bdbh; //parameter for treesize
    real btop; //parameter for topkill extent
    real bSOD; //parameter for SOD related mortality pre-fire
    real bclim; //parameter for mean annual temperature
    real<lower=0> sigma_plot;  //variance
    real phi_reciprocal;
}
transformed parameters {
vector<lower=0>[N] mu;  
vector[P] vary_a;   //varying intercept component
vector[P] pred_a;
real<lower = 0> phi;
for (j in 1:P){
pred_a[j] = bSOD*SODdead[j] + bclim*climate[j];
}   // this is the plot-level part of the model, moved as part of N-C parameterization.
vary_a = pred_a + sigma_plot*eta_a_plot;
for ( i in 1:N ) {
mu[i] = exp(a + vary_a[Plot[i]] + bdbh*LiveBAc[i] + btop*TopkillBA[i]);
}
phi = 1. / phi_reciprocal;
}
model{
a ~ normal( 0 , 5 );
bdbh ~ normal(0,5);
bclim ~ normal(0,5);
btop ~ normal(0,5);
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