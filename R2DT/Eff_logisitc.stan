
data {
  // Hyperparameters
  real gamma_mean;
  real<lower=0> gamma_sd;
  real zeta_mean;
  real<lower=0> zeta_sd;
  real eta_mean;
  real<lower=0> eta_sd;
  // Fixed trial parameters
  int<lower=1> num_doses;
  real coded_doses[num_doses]; 
// data
  int<lower=1> num_patients;
  int<lower=1, upper=num_doses> num_used_doses; //number of doses with data 
  int<lower=0, upper=num_patients> counts[num_used_doses]; //reduced x_e to exclude doses not explored
  int<lower=1, upper=num_patients> dose_counts[num_used_doses]; //reduced n to excluding doses not explored
  vector[num_used_doses] X;
  vector[num_used_doses] X2;
 }

parameters {
  // Coefficients in efficacy logit model:
  real gamma;
  real zeta;
  real eta;
}

model {
//priors 
gamma ~ normal(gamma_mean, gamma_sd);
zeta ~ normal(zeta_mean, zeta_sd);
eta ~ normal(eta_mean, eta_sd);

//liklihood
counts ~ binomial_logit(dose_counts, gamma + zeta*X + eta*X2);
}
 
generated quantities{
  real<lower=0, upper=1>p_e[num_doses];
for (i in 1:num_doses) {
p_e[i] = inv_logit(gamma + zeta * coded_doses[i] + eta * coded_doses[i] * coded_doses[i] );
}   

}














