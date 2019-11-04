data {
  int<lower=0> nIV;
  int<lower=0> nOral;
  real<lower=0> doseIV;
  real<lower=0> doseOral;
  vector<lower=0>[nIV] timeIV;
  vector<lower=0>[nIV] concIV;
  vector<lower=0>[nOral] timeOral;
  vector<lower=0>[nOral] concOral;
}
parameters {
  real<lower=0> CL;
  real<lower=0> V;
  real<lower=0,upper=1> F;
  real<lower=0> ka; 
  real<lower=0> kIV; 
  real<lower=0> cIV; 
  real<lower=0,upper=100> sigma;
}
transformed parameters {
  real k;
  real c0;
  real AUCIV;
  real c0star;
  real AUCOral;
  real ta0_5;
  real t0_5;
  vector[nIV] predIV;
  vector[nOral] predOral;
  
  k = CL / V;
  c0 = doseIV / V;
  AUCIV = doseIV / CL + cIV / kIV;
  c0star = doseOral * (ka / (ka - k)) * F / V;
  AUCOral = c0star / k;
  ta0_5 = log(2) / ka;
  t0_5 = log(2) / k;
  predIV = c0 * exp(-k * timeIV) + cIV * exp(-kIV * timeIV);
  predOral = c0star * (exp(-k * timeOral) - exp(-ka * timeOral));
}
model {
  // IV component
  kIV ~ normal(0.4, 1);
  cIV ~ lognormal(1,10);
  V ~ lognormal(2,10);
  CL ~ lognormal(1,10);
  concIV ~ normal(predIV, sigma);
  
  // oral component
  ka ~ normal(0.4,1);
  concOral ~ normal(predOral, sigma);
}