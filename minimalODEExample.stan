functions {
  vector odefunc(real t, vector y, real theta) {
    vector[2] dydt;

    dydt[1] = y[2];
    dydt[2] = -theta * y[1];
    return dydt;
  }
}

data {
  int<lower=1> T; //num timesteps
  real t0; //initial_time
  real theta; 
  array[T-1] real times;
  vector[2] y0_prior_mu;
  vector[2] y0_prior_std;
  real sigma; //measurement noise
  vector[T] z;
}

parameters {
  vector[2] y0_raw;
}

transformed parameters {
  vector[2] y0;
  y0 = y0_raw .* y0_prior_std + y0_prior_mu;
}

model {
  array[T] vector[2] y;
  
  y[1] = y0;
  y[2:T] = ode_rk45(odefunc, y0, t0, times, theta);

  for (t in 1:T) {
     z[t] ~ normal(y[t][1],sigma);
  }
}
