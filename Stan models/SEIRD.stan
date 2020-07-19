// Reference: https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real E = y[2];
      real I = y[3];
      real R = y[4];
      real D = y[5];
      real N = x_i[1];
      
      real beta = theta[1];
      real delta = theta[2];
      
      real dS_dt = -beta * I * S / N;
      real dE_dt = beta * I * S / N - 0.2 * E;
      real dI_dt =  0.2 * E - 0.2 * I-delta*I;
      real dR_dt =  0.2 * I;
      real dD_dt = delta*I;
      
      return {dS_dt, dE_dt, dI_dt, dR_dt, dD_dt};
  }
}

data {
  int<lower=1> n_days;
  real y0[5];
  real t0;
  real ts[n_days];
  int N;
  int infected[n_days];
  int recovered[n_days];
  int deaths[n_days];
}

transformed data {
  real x_r[0];
  int x_i[1]={N};
}

parameters {
  real<lower=0> beta;
  real<lower=0> delta;
  real<lower=0> phi_invI;
  real<lower=0> phi_invR;
  real<lower=0> phi_invD;
}

transformed parameters{
  real y[n_days, 5];
  real phiI = 1. / phi_invI;
  real phiR = 1. / phi_invR;
  real phiD = 1. / phi_invD;
  {
    real theta[2];
    theta[1] = beta;
    theta[2]=delta;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}

model {
  //priors
  beta ~ normal(2, 1);
  delta ~ uniform(0,1);
  phi_invI ~ exponential(5);
  phi_invR ~ exponential(5);
  phi_invD ~ exponential(5);
  
  //sampling distribution
  infected ~ neg_binomial_2(col(to_matrix(y), 3), phiI);
  recovered ~ neg_binomial_2(col(to_matrix(y), 4), phiR);
  deaths ~ neg_binomial_2(col(to_matrix(y), 5), phiD);
}

generated quantities {
  real R0 = beta / 0.2;
  real pred_cases[n_days];
  real pred_recovered[n_days];
  real pred_deaths[n_days];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 3) + 1e-5, phiI);
  pred_recovered = neg_binomial_2_rng(col(to_matrix(y), 4) + 1e-5, phiR);
  pred_deaths = neg_binomial_2_rng(col(to_matrix(y), 5) + 1e-5, phiD);
}
