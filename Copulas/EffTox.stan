// Dose-Finding Based on Efficacy-Toxicity Trade-Offs, by Thall, Cook
// Effective sample size for computing prior hyperparameters in Bayesian phase I-II dose-finding,
//  by Thall, Herrick, Nguyen, Venier, Norris

functions {
  real log_joint_pdf(real[] coded_doses, real[] coded_doses_squ,
                     int num_patients, int[] eff, int[] tox, int[] doses,
                     real alpha, real beta, real gamma, real zeta, real eta, real psi) {
    real p;
    p = 0;
    for(j in 1:num_patients) {
      real prob_eff;
      real prob_tox;
      real p_j;
      prob_eff = inv_logit(gamma + zeta * coded_doses[doses[j]] + eta * coded_doses_squ[doses[j]]);
      prob_tox = inv_logit(alpha + beta * coded_doses[doses[j]]);
      p_j = prob_eff^eff[j] * (1. - prob_eff)^(1. - eff[j]) * prob_tox^tox[j] *
              (1. - prob_tox)^(1. - tox[j]) + (-1.)^(eff[j] + tox[j]) * prob_eff *
              prob_tox * (1. - prob_eff) * (1. - prob_tox) *
              psi;
      p = p + log(p_j);
    }
    return p;
  }
}

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
  real  psi_l;
  real  psi_u;
  // Fixed trial parameters
  int<lower=1> num_doses;

  // Observed trial outcomes
  int<lower=0> num_patients;
  int<lower=0, upper=1> eff[num_patients]; // Binary efficacy event for patients j=1,..,num_patients
  int<lower=0, upper=1> tox[num_patients]; // Binary toxicity event for patients j=1,..,num_patients
  int<lower=1, upper=num_doses> doses[num_patients]; 
  real coded_doses[num_doses];
  real coded_doses_squ[num_doses]; // The square of coded_doses
}


parameters {
  // Coefficients in toxicity logit model:
  real alpha;
  //real beta;
  real  beta;
  // Coefficients in efficacy logit model:
  // Why not delta and epsilon? They have such common, alternative usage that I skipped them.
  real gamma;
  real zeta;
  real eta;
  // Association:
  real<lower=-1, upper=1> psi;
}

transformed parameters {
  real<lower=0, upper=1> prob_eff[num_doses]; // Posterior probability of efficacy at doses i=1,...,num_doses
  real<lower=0, upper=1> prob_tox[num_doses]; // Posterior probability of toxicity at doses i=1,...,num_doses

  for(i in 1:num_doses)
  {
    prob_tox[i] = inv_logit(alpha + beta * coded_doses[i]);
    prob_eff[i] = inv_logit(gamma + zeta * coded_doses[i] + eta * coded_doses_squ[i]);

  }
}

model {
  target += normal_lpdf(alpha | alpha_mean, alpha_sd);
   target += normal_lpdf(beta | beta_mean, beta_sd);
   target += normal_lpdf(gamma | gamma_mean, gamma_sd);
  target += normal_lpdf(zeta | zeta_mean, zeta_sd);
  target += normal_lpdf(eta | eta_mean, eta_sd);
  target += uniform_lpdf(psi | psi_l, psi_u);
  
  target += log_joint_pdf(coded_doses, coded_doses_squ, num_patients, eff, tox, doses,
                          alpha, beta, gamma, zeta, eta, psi);
}

