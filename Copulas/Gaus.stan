functions{

real binormalcdf(real z1, real z2, real rho) {
    if (z1 != 0 || z2 != 0) {
    real denom = fabs(rho) < 1.0 ? sqrt((1 + rho) * (1 - rho)) : not_a_number();
    real a1 = (z2 / z1 - rho) / denom;
    real a2 = (z1 / z2 - rho) / denom;
    real product = z1 * z2;
    real delta = product < 0 || (product == 0 && (z1 + z2) < 0);
    return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2);
  }
  return 0.25 + asin(rho) / (2 * pi());
}

  
    real gaus_cop(int y1, int y2, real p1, real p2, real theta){
    if (y1>=0 && y2>=0){
    real u1;
    real u2;
    u1=bernoulli_cdf(y1, p1);
    u2=bernoulli_cdf(y2, p2);
    if (u1==1 || u2==1) return(u1*u2);
    else return (binormalcdf(inv_Phi(u1),inv_Phi(u2),theta)); 
    }
    else return(0);
  }
  real joint(int eff, int tox, real prob_eff, real prob_tox, real theta){
        return (gaus_cop(eff,tox,prob_eff,prob_tox,theta) - gaus_cop((eff-1),tox,prob_eff,prob_tox,theta) 
            -gaus_cop(eff,(tox-1),prob_eff,prob_tox,theta) + gaus_cop((eff-1),(tox-1),prob_eff,prob_tox,theta)); 
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
  real <lower=-1, upper=1> psi;
}


model {
vector[num_patients] Lik;   
  
  pe ~ beta(peshape1, peshape2);
  pt ~ beta(ptshape1, ptshape2);
  psi ~ uniform(psi_l,psi_u)  ;

   for (i in 1:num_patients)
  Lik[i] = log(joint(eff[i],tox[i], pe, pt, psi));

  target += sum(Lik);
  

}

generated quantities {
 real<lower=-1, upper=1> tau;
  tau = ( joint(0, 0, pe, pt, psi) -  (1. - pe)*(1. - pt))/sqrt(pe * pt * (1. - pe) * (1. - pt));
  
}

