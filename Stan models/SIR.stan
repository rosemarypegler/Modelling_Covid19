// Reference: https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
               
      // declare the components of the model
      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];
      
      // the parameter that is going to be estimated
      real beta = theta[1];
      
      // declare the differential equations
      // this model uses a fixed recovery rate of 5 days (gamma = 0.2) due to previous literature results
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - 0.2 * I;
      real dR_dt =  0.2 * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days; // number of days to estimate for
  real y0[3]; // initial conditions
  real t0; // start day
  real ts[n_days]; // the remaining days
  int N; // total population
  int infected[n_days]; // infected data
  int removed[n_days]; // removed data
}

transformed data {
  real x_r[0];
  int x_i[1]={N};
}

parameters {
  real<lower=0> beta; // parameter to be estimated, beta
  real<lower=0> phi_invI; // defines the inverse of the overdispersion parameter for infected
  real<lower=0> phi_invR; // defines the inverse of the overdispersion parameter for removed
}

transformed parameters{
  real y[n_days, 3]; // creates an array to store the integration results
  real phiI = 1. / phi_invI; // defines the overdispersion parameter for infected
  real phiR = 1. / phi_invR; // defines the overdispersion parameter for removed
  {
    real theta[1];
    theta[1] = beta;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);// integrates the ODEs and stores the result
  }
}

model {
  //priors
  beta ~ normal(2, 1); 
  phi_invI ~ exponential(5);
  phi_invR ~ exponential(5);
  
  //sampling distribution
  // chosen a negative binomial distribution due to the number of cases being count data
  infected ~ neg_binomial_2(col(to_matrix(y), 2), phiI);
  removed ~ neg_binomial_2(col(to_matrix(y), 3), phiR);
}

generated quantities {
  real R0 = beta / 0.2; // basic reproductive ratio
  real pred_infected[n_days]; // predicted infected
  real pred_removed[n_days]; // predicted removed
  pred_infected = neg_binomial_2_rng(col(to_matrix(y), 2) + 1e-6, phiI);
  pred_removed = neg_binomial_2_rng(col(to_matrix(y), 3) + 1e-6, phiR);
}
