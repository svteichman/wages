data {
	int<lower=0> N; 
	int<lower=0> T;
	int<lower=0> K;
	real y[N,T];
	real x[N,T];
	int a[N,T];
	int<lower=0,upper=1> r[T-1,N];
	int<lower=0,upper=1> s[T-1,N];
	int<lower=0,upper=1> l[T-1,N];
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
//	real<lower=0,upper=1> p[N,T,K];
	real<lower=0> sigmasq_d;
	real delta;
	real<lower=0> sigmasq_b;
	real beta;
	real<lower=0> sigmasq_u;
	real<lower=0> sigmasq_v;
	real u[K];
	real v[K];
	real<lower=0,upper=1> eta[N,T];
	real<lower=0,upper=1> lambda[N,T];
	real<lower=0> sigmasq_k[K];
	real<lower=0> sigmasq_t;
	real<lower=0> sigmasq_r;
	real<lower=0> sigmasq_ik[N,T];
	real mu_k[T,K];
//	real<lower=0> sigmasq_theta;
//	real<lower=0> sigmasq_psi;
//	real theta
//	real psi
	real mu_ik[N,T];
}
transformed parameters {
	real<lower=0> sigma_r;
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
	real<lower=0> sigma_ik[N,T];
	real<lower=0,upper=1> p[N,T,K];

	sigma_r = sqrt(min(sigmasq_k)) - 0.1;
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
	sigma_ik = sqrt(sigmasq_ik);
	for (n in 1:N) {
		for (t in 1:T) {
			for (k in 1:K) {
				p[n,t,k] = logit(beta_k[k] - xi*distance(to_vector(z[n]),to_vector(w[k,t])) + alpha*x[n,t]);
			}
			p[n,t] = as.vector(p[n,t])/sum(p[n,t];
		}
	}

}
model{
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
//	sigmasq_theta ~ scaled_inv_chi_square(1,1);
//	sigmasq_psi ~ scaled_inv_chi_square(1,1);
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
		a[n,1] ~ multinomial(p[n,t]);
	}
	for (t in 2:T) {
		for (n in 1:N) {
			eta[n,t] = logit(beta*x[n,t] + v[a[n,t-1]]);
			lambda[n,t] = logit(delta*x[n,t] + u[a[n,t-1]]);
			r[t-1,n] ~ bernoulli(eta[n,t]);
			temp ~ bermoulli(lambda[n,t])
			s[t-1,n] ~ (1-r[t-1,n])*temp;
			l[t-1,n] = (1-r[t-1,n])*(1-s[t-1,n]);
			temp ~ multinomial(p[n,t]);
			a[n,t] ~ (r[t-1,n]+s[t-1,n])*a[n,t] + l[t-1,n]*temp;
		}
	}
	for (n in 1:N) {
		for (t in 1:T) {
			sigmasq_ik[n,t] = sigmasq_k[a[n,t]] - r[t-1,n]*sigmasq_r + l[t-1,n]*sigmasq_t;
		}
	}
	for (k in 1:K) {
		mu_k[k,1] = mean(log(y[,1]));
	}
	for (n in 1:N) {
		mu_ik[n,1] = mu_k[a[n,1],1];
		log(y[n,1]) ~ normal(mu_ik[n,1],sigma_ik[n,1]);
	}
//	psi ~ normal(0.0,sigmasq_psi);
//	theta ~ normal(0.0,sigmasq_theta);
	for (t in 2:T) {
//		for (k in 1:K) {
//			sum = 0
//			count = 0
//			for (n in 1:N)
//				if (a[n,t-1] == k)
//					sum = sum + y[n,t-1];
//					count = count + 1;
//			mu_k[t,k] = sum/count;
		}
		for (n in 1:N) {
			mu_ik[n,t] = log(y[n,t-1]);
			log(y[n,t]) ~ normal(mu_ik[n,t],sigma_ik[n,t]);
		}
	}
}
