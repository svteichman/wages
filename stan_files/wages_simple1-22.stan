// copied from 11-5
// old version, uses s|r and t|s,r instead of r,s,t drawn from multinomial
// runs, but occasionally initial values produce likelihood of 0 
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

	for (n in 1:N) {
		for (t in 1:T) {
			for (k in 1:K) {
				if (a[n,t,k]==1)
					job[n,t] = k;
			}
		}
	}
	for (n in 1:N) {
		for (t in 1:T) {
			ly[n,t] = log(y[n,t]);
		}
	}
	for (n in 1:N) {
		for (t in 1:T) {
			adj_x[1,t+(n-1)*T] = x[n,t,1];
			adj_x[2,t+(n-1)*T] = x[n,t,2];
		}
	}
}
parameters{
	real<lower=0> sigmasq_w0;
	real<lower=0> sigmasq_w;
	real<lower=0> sigmasq_z;
	real w[K,T,2];
	real z[N,2];
	real<lower=0> sigmasq_bk;
	real beta_k[K];
	real<lower=0> sigmasq_a;
	vector[Xnum] alpha;
	real<lower=0> sigmasq_d;
	vector[Xnum] delta;
	real<lower=0> sigmasq_b;
	vector[Xnum] beta;
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
	real<lower=0> sigma_w0;
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
	vector[Xnum] adj_alpha;
	vector[Xnum] adj_beta;
	vector[Xnum] adj_delta;
	vector[K] p_vec;
	vector[K] new_p_vec;

	sigma_r = fabs(sqrt(min(sigmasq_k)) - 0.1);
	sigmasq_r = sigma_r^2;
	sigma_w0 = sqrt(sigmasq_w0);
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
	adj_alpha[1] = alpha[1]/mean(adj_x[1]);
	adj_alpha[2] = alpha[2]/mean(adj_x[2]);
	adj_beta[1] = beta[1]/mean(adj_x[1]);
	adj_beta[2] = beta[2]/mean(adj_x[2]);
	adj_delta[1] = delta[1]/mean(adj_x[1]);
	adj_delta[2] = delta[2]/mean(adj_x[2]);
	for (n in 1:N) {
		for (t in 1:T) {
			for (k in 1:K) {
				p[n,t,k] = inv_logit(beta_k[k] - distance(to_vector(z[n]),to_vector(w[k,t])) + dot_product(adj_alpha',to_vector(x[n,t])));
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

	sigmasq_w0 ~ scaled_inv_chi_square(1,5);
	sigmasq_w ~ scaled_inv_chi_square(1,25);
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
		w[k,1] ~ normal(0.0,sigma_w0);
		for (t in 2:T) {
			w[k,t] ~ normal(w[k,t-1],sigma_w);
		}
	}
	for (n in 1:N) {
		z[n] ~ normal(0.0,sigma_z);
	}
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
			r[t-1,n] ~ bernoulli(inv_logit(dot_product(adj_beta',to_vector(x[n,t])) + v[job[n,t-1]]));
			if (r[t-1,n]==0) {
				s[t-1,n] ~ bernoulli(inv_logit(dot_product(adj_delta',to_vector(x[n,t])) + u[job[n,t-1]]));
			}
			else {
				s[t-1,n] ~ bernoulli(0.0000001);
			}
			if (1-r[t-1,n]*(1-s[t-1,n])==0) {
				probs = rep_vector(0.0000001,K);
				probs[job[n,t-1]] = 1.0 - (K-1)*0.0000001;
				a[n,t] ~ multinomial(probs);
			}
			else {
				a[n,t] ~ multinomial(to_vector(p[n,t]));
			}
		}
	}
	for (n in 1:N) {
		ly[n,1] ~ normal(first_class_means[job[n,1]],sigma_ik[n,1]);
	}
	for (t in 2:T) {
		for (n in 1:N) {
			ly[n,t] ~ normal(ly[n,t-1],sigma_ik[n,t]);
		}
	}
}
