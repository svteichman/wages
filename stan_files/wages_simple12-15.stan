// copied from 9-30
// fixing 3 sectors in the first time point 
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
  real z[N,2];
 // real<lower=0> tau_b;
  real beta;
  real<lower=0> tau;
  real<lower=0> tau_m;
  real mu[N];
  real nu[K];
}
transformed parameters {
  real w[K,T,2];
  real p[N,T,K];
  real<lower=0> sigmasq_i[N,T];
  real<lower=0> sigma_i[N,T];
  w[1,1,1] = 0;
  w[1,1,2] = 0;
  w[2,1,1] = 1;
  w[2,1,2] = 0;
  w[3,1,1] = 1.0/2.0;
  w[3,1,2] = sqrt(3.0)/2.0;
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
  for (n in 1:N) {
    for (k in 1:K) {
      p[n,1,k] =  -distance(to_vector(z[n]),to_vector(w[k,1]));
      for (t in 2:T) {
        p[n,t,k] = a[n,t-1,k]*(beta) - distance(to_vector(z[n]),to_vector(w[k,t]));
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
  for (k in 1:(K-3)) {
    free_w1[k,1] ~ normal(0.0,1/sqrt(tau_w0));
  }
  for (k in 1:K) {
    for (t in 2:T) {
      free_w[k,t-1] ~ normal(free_w[k,t-1],1/sqrt(tau_w));
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
  beta ~ normal(1, sqrt(0.1));
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
