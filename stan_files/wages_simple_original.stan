data {
	int<lower=0> N; 
	int<lower=0> T;
	int<lower=0> K;
	real y[N,T];
	real x[N,T];
	int a[N,T,K];
	int<lower=0,upper=1> r[T-1,N];
	int<lower=0,upper=1> s[T-1,N];
	real<lower=0> first_class_means[K];
}
transformed data {
	int job[N,T];

	for (n in 1:N) {
		for (t in 1:T) {
			for (k in 1:K) {
				if (a[n,t,k]==1)
					job[n,t] = k;
			}
		}
	}
}
parameters{
	real<lower=0> sigmasq_w;
	real<lower=0> sigmasq_z;
	real w[K,T,2];
	real z[N,2];
	real<lower=0> xi;
	real<lower=0> sigmasq_bk;
	real beta_k[K];
	real<lower=0> sigmasq_a;
	real alpha;
	real<lower=0> sigmasq_d;
	real delta;
	real<lower=0> sigmasq_b;
	real beta;
	real<lower=0> sigmasq_u;
	real<lower=0> sigmasq_v;
	real u[K];
	real v[K];
//	real<lower=0,upper=1> eta[N,T-1];
//	real<lower=0,upper=1> lambda[N,T-1];
	real<lower=0> sigmasq_k[K];
	real<lower=0> sigmasq_t;
}
transformed parameters {
	real<lower=0> sigma_r;
	real<lower=0> sigmasq_r;
	real<lower=0> sigma_w;
	real<lower=0> sigma_z;
	real<lower=0> sigma_bk;
	real<lower=0> sigma_a;
	real<lower=0> sigma_d;
	real<lower=0> sigma_b;
	real<lower=0> sigma_u;
	real<lower=0> sigma_v;
	real<lower=0> sigma_k[K];
	real<lower=0> sigma_t;
	real<lower=0,upper=1> p[N,T,K];
	real<lower=0> sigmasq_ik[N,T];
	real<lower=0> sigma_ik[N,T];
	vector[K] p_vec;
	vector[K] new_p_vec;

	sigma_r = fabs(sqrt(min(sigmasq_k)) - 0.1);
	sigmasq_r = sigma_r^2;
	sigma_w = sqrt(sigmasq_w);
	sigma_z = sqrt(sigmasq_z);
	sigma_bk = sqrt(sigmasq_bk);
	sigma_a = sqrt(sigmasq_a);
	sigma_d = sqrt(sigmasq_d);
	sigma_b = sqrt(sigmasq_b);
	sigma_u = sqrt(sigmasq_u);
	sigma_v = sqrt(sigmasq_v);
	sigma_k = sqrt(sigmasq_k);
	sigma_t = sqrt(sigmasq_t);
	for (n in 1:N) {
		for (t in 1:T) {
			for (k in 1:K) {
				p[n,t,k] = inv_logit(beta_k[k] - xi*distance(to_vector(z[n]),to_vector(w[k,t])) + alpha*x[n,t]);
			}
			p_vec = to_vector(p[n,t]);
			new_p_vec = p_vec/sum(p[n,t]);
			p[n,t] = to_array_1d(new_p_vec);
		}
	}
	for (n in 1:N) {
		sigmasq_ik[n,1] = sigmasq_k[job[n,1]] - sigmasq_r;
		for (t in 2:T) {
			sigmasq_ik[n,t] = sigmasq_k[job[n,t]] - r[t-1,n]*sigmasq_r + (1-r[t-1,n])*(1-s[t-1,n])*sigmasq_t;
		}
	}
	sigma_ik = sqrt(sigmasq_ik);

}
model{
	vector[K] probs;

	sigmasq_w ~ scaled_inv_chi_square(1,5);
	sigmasq_z ~ scaled_inv_chi_square(1,5);
	sigmasq_bk ~ scaled_inv_chi_square(1,5);
	sigmasq_a ~ scaled_inv_chi_square(1,5);
	sigmasq_d ~ scaled_inv_chi_square(1,13);
	sigmasq_b ~ scaled_inv_chi_square(1,13);
	sigmasq_u ~ scaled_inv_chi_square(1,3);
	sigmasq_v ~ scaled_inv_chi_square(1,3);
	sigmasq_k ~ scaled_inv_chi_square(1,1); 
	sigmasq_t ~ scaled_inv_chi_square(1,4);
	for (k in 1:K) {
		for (t in 1:T) {
			w[k,t] ~ normal(0.0,sigma_w);
		}
	}
	for (n in 1:N) {
		z[n] ~ normal(0.0,sigma_z);
	}
	xi ~ gamma(1,1);
	beta_k ~ normal(0.0,sigma_bk); 
	alpha ~ normal(0.0,sigma_a);
	delta ~ normal(0.0,sigma_d);
	beta ~ normal(0.0001,sigma_b);
	u ~ normal(0.0, sigma_u);
	v ~ normal(0.0, sigma_v);
	for (n in 1:N) {
		a[n,1] ~ multinomial(to_vector(p[n,1]));
	}
	for (t in 2:T) {
		for (n in 1:N) {
			r[t-1,n] ~ bernoulli(inv_logit(beta*x[n,t] + v[job[n,t-1]]));
			if (r[t-1,n]==0) {
				s[t-1,n] ~ bernoulli(inv_logit(delta*x[n,t] + u[job[n,t-1]]));
			}
			else {
				s[t-1,n] ~ bernoulli(0);
//				s[t-1,n] ~ bernoulli(0.5);
			}
			if (1-r[t-1,n]*(1-s[t-1,n])==0) {
				probs = rep_vector(0.0,K);
				probs[job[n,t-1]] = 1.0;
				a[n,t] ~ multinomial(probs);
//				a[n,t] ~ multinomial(to_vector(p[n,t]));
			}
			else {
				a[n,t] ~ multinomial(to_vector(p[n,t]));
			}
		}
	}
	for (n in 1:N) {
		y[n,1] ~ lognormal(first_class_means[job[n,1]],sigma_ik[n,1]);
	}
	for (t in 2:T) {
		for (n in 1:N) {
			y[n,t] ~ lognormal(log(y[n,t-1]),sigma_ik[n,t]);
		}
	}
}
