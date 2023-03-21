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
  real sigma; //measurement noise
  vector[2] y0;
}

generated quantities {
  array[T] vector[2] y;
  vector[T] z;
  y[1] = y0;
  y[2:T] = ode_rk45(odefunc, y0, t0, times, theta);

  for (t in 1:T) {
     z[t] = normal_rng(y[t][1],sigma);
  }
}
