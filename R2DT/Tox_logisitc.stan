
data {
  // Hyperparameters
  real alpha_mean;
  real<lower=0> alpha_sd;
  real beta_mean;
  real<lower=0> beta_sd;

  // Fixed trial parameters
  int<lower=1> num_doses;
  real coded_doses[num_doses]; 
// data
  int<lower=1> num_patients;
  int<lower=1, upper=num_doses> num_used_doses; //number of doses with data 
  int<lower=0, upper=num_patients> counts[num_used_doses]; //reduced x_t to excluding doses not explored
  int<lower=1, upper=num_patients> dose_counts[num_used_doses]; //reduced n to excluding doses not explored
  vector[num_used_doses] X;
 }

parameters {
  // Coefficients in toxicity logit model:
  real alpha;
  real beta;

}

model {
//priors 
alpha ~ normal(alpha_mean, alpha_sd);
beta ~ normal(beta_mean, beta_sd);


//liklihood
counts ~ binomial_logit(dose_counts, alpha + beta*X);
}
 
generated quantities{
  real<lower=0, upper=1>p_t[num_doses];
for (i in 1:num_doses) {
p_t[i] = inv_logit(alpha + beta * coded_doses[i]);
}   

}














