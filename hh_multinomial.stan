//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real Sk_Sa = y[1];
      real Sk_Ia = y[2];
      real Ik_Sa = y[3];
      real Ik_Ia = y[4];
      real N = x_i[1];
      
      real lambda_k = theta[1];
      real lambda_a = theta[2];
      real mu_k = theta[3];
      real mu_a = theta[4];
      real delta_a_k = theta[5];
      real delta_k_a = theta[6];
      
      real dSk_Sa_dt = -lambda_k*Sk_Sa + mu_a*Sk_Ia + mu_k*Ik_Sa  ;
      real dSk_Ia_dt =  lambda_a*Sk_Sa - mu_a*Sk_Ia - (lambda_k + delta_a_k)*Sk_Ia  ;
      real dIk_Sa_dt =  lambda_k*Sk_Sa - mu_k*Ik_Sa - (lambda_a + delta_k_a)*Ik_Sa  ;
      real dIk_Ia =  (lambda_a + delta_k_a)*Ik_Sa + (lambda_k + delta_a_k)*Sk_Ia - Ik_Ia*(mu_k + mu_a) ;
      
      return {dSk_Sa_dt,dSk_Ia_dt,dIk_Sa_dt,dIk_Ia};
  }
}


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> n_days;
  real y0[4];
  real t0;
  real ts[n_days];
  int N;
  int cases[4];
}

transformed data {
  real x_r[0];
  int x_i[1] = { N };
}

parameters {
  real<lower=0> lambda_k;
  real<lower=0> lambda_a;
  real<lower=0> mu_k;
  real<lower=0> mu_a;
  real<lower=0> delta_a_k;
  real<lower=0> delta_k_a;
}

transformed parameters{
  real y[n_days, 4];
  vector[4] y_last;  
  real sum_last;
  vector[4] y_prob_last; 
  {
    real theta[6];
    theta[1] = lambda_k;
    theta[2] = lambda_a;
    theta[3] = mu_k;
    theta[4] = mu_a;
    theta[5] = delta_a_k;
    theta[6] = delta_k_a;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
    y_last = to_vector(y[n_days,:]);
    sum_last = sum(y_last);
    y_prob_last = y_last/sum_last;
  }
}

model {
  //priors
  lambda_k ~ normal(0.5, 0.5); //truncated at 0, close to uniform on [0,1]
  lambda_a ~ normal(0.5, 0.5); //truncated at 0, close to uniform on [0,1]
  mu_k ~ normal(0.016667, 0.01); //truncated at 0, close to uniform on [0,1]
  mu_a ~ normal(0.07142, 0.5); //truncated at 0, close to uniform on [0,1]
  delta_a_k ~ normal(0.5, 0.5); //truncated at 0, close to uniform on [0,1]
  delta_k_a ~ normal(0.5, 0.5); //truncated at 0, close to uniform on [0,1]

  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people
  cases ~ multinomial(y_prob_last);
}
