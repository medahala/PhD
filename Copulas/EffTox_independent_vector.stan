
data {
  // Hyperparameters
  real alpha_mean;
  real<lower=0> alpha_sd;
  real beta_mean;
  real<lower=0> beta_sd;
  real gamma_mean;
  real<lower=0> gamma_sd;
  real zeta_mean;
  real<lower=0> zeta_sd;
  real eta_mean;
  real<lower=0> eta_sd;
  // Fixed trial parameters
  int<lower=1> num_doses;
  int<lower=0> x_e[num_doses]; // summed binary data for efficacious event
  int<lower=0> x_t[num_doses]; // summed toxicity event 
  int<lower=0> n[num_doses]; //number of patients treated at each dose
  vector[num_doses] X; //vector for X
  vector[num_doses] X2; //vector for X**2  
  
}


parameters {
  // Coefficients in toxicity logit model:
  real alpha;
  real beta;
  // Coefficients in efficacy logit model:
  real gamma;
  real zeta;
  real eta;
}

model {
//priors 
alpha ~ normal(alpha_mean, alpha_sd);
beta ~ normal(beta_mean , beta_sd);
gamma ~ normal(gamma_mean, gamma_sd);
zeta ~ normal(zeta_mean , zeta_sd);
eta ~ normal(eta_mean, eta_sd);

//liklihood
x_t ~ binomial_logit(n, alpha+beta*X);
x_e ~ binomial_logit(n, gamma+zeta*X+eta*X2);
}

generated quantities{
  real<lower=0, upper=1>p_e[num_doses];
  real<lower=0, upper=1>p_t[num_doses];
for (i in 1:num_doses) {
p_t[i] = inv_logit(alpha + beta * X[i]);
p_e[i] = inv_logit(gamma + zeta * X[i] + eta * X2[i]);
}   
  
}







