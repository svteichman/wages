// copied from 7-3
// adding a variance component for switching jobs 
data {
  int<lower=0> N; 
  int<lower=0> T;
  int<lower=0> K;
  int<lower=0> Xnum;
  real y[N,T];
  real x[N,T,Xnum];
  int a[N,T,K];
  int<lower=0,upper=1> r[T-1,N];
  int<lower=0,upper=1> s[T-1,N];
  real<lower=0> first_class_means[K];
}
transformed data {
  int job[N,T];
  real ly[N,T];
  real adj_x[Xnum, N*T];
  real x_norm[N,T,Xnum];
  for (n in 1:N) {
    for (t in 1:T) {
      for (k in 1:K) {
        if (a[n,t,k]==1)
          job[n,t] = k;
      }
      ly[n,t] = log(y[n,t]);
      adj_x[1,t+(n-1)*T] = x[n,t,1];
      adj_x[2,t+(n-1)*T] = x[n,t,2];
    }
  }
  // normalize X by subtracting the mean and dividing by the variance 
  for (n in 1:N) {
    for (t in 1:T) {
      x_norm[n,t,1] = (x[n,t,1] - mean(adj_x[1])/sd(adj_x[1]));
      x_norm[n,t,2] = (x[n,t,2] - mean(adj_x[2])/sd(adj_x[2]));
    }
  }
}
parameters{
  real<lower=0> tau_w0;
  real<lower=0> tau_w;
  real w[K,T,2];
  real z[N,2];
  real<lower=0> tau_bg;
  real beta;
  real<lower=0> tau_k[K];
  real<lower=0> tau_m;
}
transformed parameters {
  real<lower=0> sigmasq_k[K];
  real<lower=0> sigma_k[K];
  real p[N,T,K];
  real<lower=0> sigmasq_i[N,T];
  real<lower=0> sigma_i[N,T];
  real<lower=0> sigmasq_m;
  real<lower=0> sigma_m;

  sigmasq_m = 1/tau_m;
  sigma_m = sqrt(sigmasq_m);
  for (k in 1:K) {
    sigmasq_k[k] = 1/tau_k[k];
    sigma_k[k] = sqrt(sigmasq_k[k]);
  }

  for (n in 1:N) {
    for (k in 1:K) {
      p[n,1,k] =  -distance(to_vector(z[n]),to_vector(w[k,1]));
      for (t in 2:T) {
        p[n,t,k] = a[n,t-1,k]*(-beta) + (1-a[n,t-1,k])*beta - distance(to_vector(z[n]),to_vector(w[k,t]));
      }
    }
  }
  for (n in 1:N) {
    sigmasq_i[n,1] = sigmasq_k[job[n,1]];
    for (t in 2:T) {
      sigmasq_i[n,t] = sigmasq_k[job[n,t]] + (1-r[t-1,n])*sigmasq_m;
    }
  }
  sigma_i = sqrt(sigmasq_i);
}
model{
  tau_w0 ~ gamma(.5, 1);
  tau_w ~ gamma(.5, 1);
  tau_bg ~ gamma(.5, 1);
  tau_m ~ gamma(1, 1);
  // latent positions 
  for (k in 1:K) {
    w[k,1] ~ normal(0.0,1/sqrt(tau_w0));
    for (t in 2:T) {
      w[k,t] ~ normal(w[k,t-1],1/sqrt(tau_w));
    }
  }
  for (n in 1:N) {
    z[n] ~ normal(0.0,1/sqrt(0.5)); // replaced tau_z with 0.5 here 
  }
  beta ~ normal(0.0,1/sqrt(tau_bg)); 
  for (t in 1:T) {
    for (n in 1:N) {
      job[n,t] ~ categorical_logit(to_vector(p[n,t]));
    }
  }
  tau_k ~ gamma(0.71, 1); 
  for (n in 1:N) {
    ly[n,1] ~ normal(first_class_means[job[n,1]],sigma_i[n,1]);
  }
  for (t in 2:T) {
    for (n in 1:N) {
      ly[n,t] ~ normal(ly[n,t-1],sigma_i[n,t]);
    }
  }
}