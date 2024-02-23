
data {
  // Hyperparameters
  real<lower=0> peshape1;
  real<lower=0> peshape2;
  real<lower=0> ptshape1;
  real<lower=0> ptshape2;

  // Observed trial outcomes
  int<lower=0> num_patients;
  int<lower=0, upper=1> eff[num_patients]; // Binary efficacy event for patients j=1,..,num_patients
  int<lower=0, upper=1> tox[num_patients]; // Binary toxicity event for patients j=1,..,num_patients
}

transformed data{
  int<lower=0, upper=num_patients> x_e;
  int<lower=0, upper=num_patients> x_t;
  x_e=  sum(eff);
  x_t=  sum(tox);
}

parameters {
  real<lower=0, upper=1> pe;
  real<lower=0, upper=1> pt; 
}

model {
 pe ~ beta(peshape1 , peshape2 ) ; 
 pt ~ beta(ptshape1 , ptshape2 ) ; 
 
 x_e ~ binomial(num_patients , pe);
 x_t ~ binomial(num_patients , pt);
 
}


