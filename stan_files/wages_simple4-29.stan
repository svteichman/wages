// copied from 4-25
// adding covariates to persistence parameter
data {
	int<lower=0> N; 
	int<lower=0> T;
	int<lower=0> K;
	int<lower=0> Xnum;
	real y[N,T];
	real x[N,T,Xnum];
	int a[N,T,K];
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
	real<lower=0> sigmasq_w0;
	real<lower=0> sigmasq_w;
	real w[K,T,2];
	real z[N,2];
	real<lower=0> sigmasq_b;
	real<lower=0> sigmasq_g0;
	real<lower=0> sigmasq_g;
	real<lower=0> sigmasq_a;
	real beta;
	real gamma;
	real alpha[Xnum];
	real<lower=0> sigmasq_k[K];
	real gamma_it[N,T-1];
}
transformed parameters {
	real<lower=0> sigma_w0;
	real<lower=0> sigma_w;
	real<lower=0> sigma_b;
	real<lower=0> sigma_g0;
	real<lower=0> sigma_g;
	real<lower=0> sigma_a;
	real<lower=0> sigma_k[K];
	real p[N,T,K];
	real<lower=0> sigmasq_i[N,T];
	real<lower=0> sigma_i[N,T];

	sigma_w0 = sqrt(sigmasq_w0);
	sigma_w = sqrt(sigmasq_w);
	sigma_b = sqrt(sigmasq_b);
	sigma_g0 = sqrt(sigmasq_g0);
	sigma_g = sqrt(sigmasq_g);
	sigma_a = sqrt(sigmasq_a);
	sigma_k = sqrt(sigmasq_k);

	for (n in 1:N) {
		for (k in 1:K) {
			p[n,1,k] =  -distance(to_vector(z[n]),to_vector(w[k,1]));
			for (t in 2:T) {
				p[n,t,k] = a[n,t-1,k]*gamma_it[n,t-1] + (1-a[n,t-1,k])*beta - distance(to_vector(z[n]),to_vector(w[k,t]));
			}
		}
	}
	for (n in 1:N) {
		for (t in 1:T) {
			sigmasq_i[n,t] = sigmasq_k[job[n,t]];
		}
	}
	sigma_i = sqrt(sigmasq_i);
}
model{
	sigmasq_w0 ~ scaled_inv_chi_square(1,5);
	sigmasq_w ~ scaled_inv_chi_square(1,25);
	sigmasq_b ~ scaled_inv_chi_square(1,5);
	sigmasq_b ~ scaled_inv_chi_square(1,5);
	sigmasq_g ~ scaled_inv_chi_square(1,5);
	sigmasq_a ~ scaled_inv_chi_square(1,5);
	sigmasq_k ~ scaled_inv_chi_square(1,1); 
	// latent positions 
	for (k in 1:K) {
		w[k,1] ~ normal(0.0,sigma_w0);
		for (t in 2:T) {
			w[k,t] ~ normal(w[k,t-1],sigma_w);
		}
	}
	for (n in 1:N) {
		z[n] ~ normal(0.0,sqrt(20)); // replaced sigma^2_z with 20 here 
	}
	beta ~ normal(0.0,sigma_b); 
	gamma ~ normal(0.0,sigma_g0);
	alpha ~ normal(0.0,sigma_a);
	for (n in 1:N) {
		for (t in 2:T) {
			gamma_it[n,t-1] ~ normal(gamma + dot_product(to_vector(alpha)',to_vector(x_norm[n,t])), sigma_g);
		}
	}
	for (t in 1:T) {
		for (n in 1:N) {
			job[n,t] ~ categorical_logit(to_vector(p[n,t]));
		}
	}
	for (n in 1:N) {
		ly[n,1] ~ normal(first_class_means[job[n,1]],sigma_i[n,1]);
	}
	for (t in 2:T) {
		for (n in 1:N) {
			ly[n,t] ~ normal(ly[n,t-1],sigma_i[n,t]);
		}
	}
}
