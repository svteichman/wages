// copied from 12-15
// adding intercept for probability multinomial logit 
// adding coefficient to d_ij
// make all params indexed by j, then setting Beta^j to 0 for identifiaibility of multinomial logit 
// constrain latent space by fixing 2 positions and constraining the 3rd to have a positive x value 
data {
  int<lower=0> N; 
  int<lower=0> T;
  int<lower=0> K;
  real y[N,T];
  int a[N,T,K];
  int<lower=0,upper=1> r[T-1,N];
  real mean_log_wage;
}
transformed data {
  int job[N,T];
  real ly[N,T];
  for (n in 1:N) {
    for (t in 1:T) {
      for (k in 1:K) {
        if (a[n,t,k]==1)
          job[n,t] = k;
      }
      ly[n,t] = log(y[n,t]);
    }
  }
}
parameters{
  real<lower=0> tau_w0;
  real<lower=0> tau_w;
  real<lower=0> tau_z;
  real free_w[K,T-1,2];
  real free_w1[K-3,1,2];
  real<lower=0> free_w_3_1_1;
  real free_w_3_1_2;
  real z[N,2];
 // real<lower=0> tau_b;
  real free_beta0[K-1];
  real free_beta1;
  real free_beta2;
  real<lower=0> tau;
  real<lower=0> tau_m;
  real mu[N];
  real nu[K];
}
transformed parameters {
  real w[K,T,2];
  real beta0[K];
  real beta1[K];
  real beta2[K];
  real p[N,T,K];
  real<lower=0> sigmasq_i[N,T];
  real<lower=0> sigma_i[N,T];
  w[1,1,1] = 0;
  w[1,1,2] = 0;
  w[2,1,1] = 1;
  w[2,1,2] = 0;
  w[3,1,1] = free_w_3_1_1;
  w[3,1,2] = free_w_3_1_2;
  beta0[1] = 0;
  beta1[1] = 0;
  beta2[1] = 0;
  for (k in 4:K) {
    w[k,1,1] = free_w1[(k-3),1,1];
    w[k,1,2] = free_w1[(k-3),1,2];
  }
  for (k in 1:K) {
    for (t in 2:T) {
      w[k,t,1] = free_w[k,t-1,1];
      w[k,t,2] = free_w[k,t-1,2];
    }
  }
  for (k in 2:K) {
    beta0[k] = free_beta0[k-1];
    beta1[k] = free_beta1;
    beta2[k] = free_beta2;
  }
  for (n in 1:N) {
    for (k in 1:K) {
      p[n,1,k] =  -distance(to_vector(z[n]),to_vector(w[k,1]));
      for (t in 2:T) {
        p[n,t,k] = beta0[k] + a[n,t-1,k]*beta1[k] - distance(to_vector(z[n]),to_vector(w[k,t]))*beta2[k];
      }
    }
  }
  for (n in 1:N) {
    sigmasq_i[n,1] = 1/tau;
    for (t in 2:T) {
      sigmasq_i[n,t] = 1/tau*(1 + (1-r[t-1,n])/tau_m);
    }
  }
  sigma_i = sqrt(sigmasq_i);
}
model{
  tau_w0 ~ gamma(.5, 1);
  tau_z ~ gamma(.5, 1);
  tau ~ gamma(1, 1);
  tau_w ~ gamma(.5, 1); 
//  tau_b ~ gamma(.5, 1);
  tau_m ~ gamma(1, 1);
  // latent positions 
  free_w_3_1_1 ~ gamma(.5, 1);
  free_w_3_1_2 ~ normal(0.0, 1/sqrt(tau_w0));
  for (k in 1:(K-3)) {
    free_w1[k,1] ~ normal(0.0, 1/sqrt(tau_w0));
  }
  for (k in 1:K) {
    for (t in 2:T) {
      free_w[k,t-1] ~ normal(free_w[k,t-1], 1/sqrt(tau_w));
    }
  }
  for (n in 1:N) {
    z[n] ~ normal(0.0,1/sqrt(tau_z)); 
    mu[n] ~ normal(mean_log_wage, 5); // large variance for now 
  }
  for (k in 1:K) {
    nu[k] ~ normal(0, 5); 
  }
  //beta ~ normal(0.0,1/sqrt(tau_b)); 
  free_beta0 ~ normal(1, sqrt(0.1));
  free_beta1 ~ normal(1, sqrt(0.1));
  free_beta2 ~ normal(1, sqrt(0.1));
  for (t in 1:T) {
    for (n in 1:N) {
      job[n,t] ~ categorical_logit(to_vector(p[n,t]));
    }
  }
  for (t in 1:T) {
    for (n in 1:N) {
      ly[n,t] ~ normal(mu[n] + nu[job[n,t]],sigma_i[n,t]);
    }
  }
}
