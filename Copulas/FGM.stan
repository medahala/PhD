
functions {
  real log_joint_pdf(int num_patients, int[] eff, int[] tox, real prob_eff, real prob_tox, real psi ){
    real p;
    p = 0;
    for(j in 1:num_patients) {
      real p_j;
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
  real<lower=0> peshape1;
  real<lower=0> peshape2;
  real<lower=0> ptshape1;
  real<lower=0> ptshape2;
  real psi_l;
  real psi_u;

  // Observed trial outcomes
  int<lower=0> num_patients;
  int<lower=0, upper=1> eff[num_patients]; // Binary efficacy event for patients j=1,..,num_patients
  int<lower=0, upper=1> tox[num_patients]; // Binary toxicity event for patients j=1,..,num_patients
}


parameters {
  // Coefficients in toxicity logit model:
  real <lower=0, upper=1> pe;
  real <lower=0, upper=1> pt;
  real <lower=-1, upper=1>psi;
}


model {
  
  pe ~ beta(peshape1,peshape2);
  pt ~ beta(ptshape1,ptshape2);
  psi ~ uniform(psi_l,psi_u)  ;

  target += log_joint_pdf(num_patients, eff, tox, pe, pt, psi);
  

}

generated quantities {
 real<lower=-1, upper=1> tau;
  tau =  (pe * pt * (1. - pe) * (1. - pt) * psi)/sqrt(pe * pt * (1. - pe) * (1. - pt)) ;
}


