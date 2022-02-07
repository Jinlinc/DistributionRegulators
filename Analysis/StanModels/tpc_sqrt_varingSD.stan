// stan assuming normal distribution with sqrt transformed input, varing sd among temperature treatment
functions{
  real briere_tpc(real Ti, real a, real b, real RTmin, real RTmax){
    real fec;
    if (Ti <= RTmin || Ti >= RTmax){
      fec = 0; // 0
    }else {
      fec = a*Ti*(Ti - RTmin)*(RTmax - Ti)^(1/b); 
    }
    return fec;
  }
}
data {
  int<lower=0> N;  // the totol number of rows of the data
  int<lower=0> Nsp;  // number of species = 9, including Melanogaster
  int<lower=0> Ntemps; // number of temperature treatment = 7
  vector<lower=0>[N] Temp;  // temperature treatment, length = nrow(data), real type
  matrix[N,Nsp] sp;  // species treatment represented by a matrix
  matrix[N,Ntemps] temps;  // temperatre treatment represented by a matrix
  real<lower=0> y[N];  // observation: offspring counts - integer
}
parameters {
  vector<lower=0>[Nsp] a;  // themal performance function parameters - vector!
  vector<lower=0>[Nsp] b;  
  vector<upper=17>[Nsp] RTmin;
  vector<lower=26, upper=35>[Nsp] RTmax;
  vector<lower=0>[Ntemps] sigma_y;  // assuming the variation will change among temperature 
  // super parameters
  real<lower=0> mu_a;  
  real<lower=0> mu_b;
  real<upper=17> mu_RTmin;
  real<lower=26, upper=35> mu_RTmax;
  real<lower=0> mu_sigma_y;
  real<lower=0> sigmasq_sigma_y;
  real<lower=0> sigmasq_a;
  real<lower=0> sigmasq_b;
  real<lower=0> sigmasq_RTmin;
  real<lower=0> sigmasq_RTmax;
  
}
transformed parameters {
  real<lower=0> sigma_sigma_y = sqrt(sigmasq_sigma_y);
  real<lower=0> sigma_a = sqrt(sigmasq_a);
  real<lower=0> sigma_b = sqrt(sigmasq_b);
  real<lower=0> sigma_RTmin = sqrt(sigmasq_RTmin);
  real<lower=0> sigma_RTmax = sqrt(sigmasq_RTmax);
  // calculating lambda (expected value) for each data entry
  vector[N] sigma_y_temps = temps*sigma_y;
  vector[N] a_sp = sp*a;
  vector[N] b_sp = sp*b;
  vector[N] RTmin_sp = sp*RTmin;
  vector[N] RTmax_sp = sp*RTmax;
  vector<lower=0>[N] theta;
  // these a,b,RTmax,RTmin is for estimating sqrt(offspring count) rather than offspring count; if theta = sqrt(est), the fitting doesn't converge well
  for (i in 1:N){theta[i] = briere_tpc(Temp[i], a_sp[i], b_sp[i], RTmin_sp[i], RTmax_sp[i]);}
}
model {
  // defining prior 
  mu_a ~ normal(0, 1);  // non-informative prior
  mu_b ~ normal(0, 10);
  mu_RTmin ~ normal(15, 10);
  mu_RTmax ~ normal(30, 10);
  mu_sigma_y ~ normal(0,10);
  sigmasq_sigma_y ~ inv_gamma(0.001, 0.001); // this is a non-informative prior, seems widely used to model sigmasq
  sigmasq_a ~ inv_gamma(0.001, 0.001);  
  sigmasq_b ~ inv_gamma(0.001, 0.001);
  sigmasq_RTmin ~ inv_gamma(0.001, 0.001);
  sigmasq_RTmax ~ inv_gamma(0.001, 0.001);
  
  // defining super model for estimating parameters
  a ~ normal(mu_a, sigma_a); // vectorized
  b ~ normal(mu_b, sigma_b);  // vectorized
  RTmin ~ normal(mu_RTmin, sigma_RTmin);  // vectorized
  RTmax ~ normal(mu_RTmax, sigma_RTmax);  // vectorized
  sigma_y ~ normal(mu_sigma_y, sigma_sigma_y);
  
  // defining model for estimating observation
  for (i in 1:N){
    y[i] ~ normal(theta[i], sigma_y_temps[i]);} 
}
generated quantities {
  real y_new[N];
  real RTopt[Nsp];
  vector[N] log_lik;
  for (i in 1:N) {
    y_new[i] = normal_rng(theta[i], sigma_y_temps[i]);
    log_lik[i] = normal_lpdf(y[i] | theta[i], sigma_y_temps[i]);
  }
  for (j in 1:Nsp){
    RTopt[j] = (2*b[j]*RTmax[j] + (b[j]+1)*RTmin[j] + sqrt((2*b[j]*RTmax[j])^2 + ((b[j]+1)*RTmin[j])^2 - 4*b[j]*b[j]*RTmin[j]*RTmax[j]))/(4*b[j] + 2);
  }
}
