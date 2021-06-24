// copied from 4-11
// replacing inv_logit with categorical_logit
// removing the normalizing I had to do with inv_logit 
// doing this for p and for later multinomial probabilities 
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
	// first_class_means is just the average wage in each class, as a place to start instead of randomly assigning wages in the first time point 
	real<lower=0> first_class_means[K];
}
transformed data {
	int job[N,T];
	real ly[N,T];
	real adj_x[Xnum, N*T];
	int<lower=1,upper=3> m[T-1,N];
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

	// creating categorical variable for staying vs. leaving 
	for (n in 1:N) { 
		for (t in 1:(T-1)) {
			m[t,n] = r[t,n] + 2*s[t,n] + 3*(1-r[t,n])*(1-s[t,n]);
		}
	}
}
parameters{
	real<lower=0> sigmasq_w0;
	real<lower=0> sigmasq_w;
	real w[K,T,2];
	real z[N,2];
	real<lower=0> sigmasq_bk;
	real beta_k[K];
	real<lower=0> sigmasq_dbz;
	vector[Xnum] delta;
	vector[Xnum] beta;
	vector[Xnum] zeta;
	real<lower=0> sigmasq_uvq;
	real u[K];
	real v[K];
	real q[K];
	real<lower=0> sigmasq_k[K];
}
transformed parameters {
	real<lower=0> sigma_w0;
	real<lower=0> sigma_w;
	real<lower=0> sigma_bk;
	real<lower=0> sigma_dbz;
	real<lower=0> sigma_uvq;
	real<lower=0> sigma_k[K];
	real p[N,T,K];
	real<lower=0> sigmasq_i[N,T];
	real<lower=0> sigma_i[N,T];
	real eta[N,T-1];
	real lambda[N,T-1];
	real kappa[N,T-1];

	sigma_w0 = sqrt(sigmasq_w0);
	sigma_w = sqrt(sigmasq_w);
	sigma_bk = sqrt(sigmasq_bk);
	sigma_dbz = sqrt(sigmasq_dbz);
	sigma_uvq = sqrt(sigmasq_uvq);
	sigma_k = sqrt(sigmasq_k);

// computing probabilities of connecting between each worker and firm 
	for (n in 1:N) {
		for (t in 1:T) {
			for (k in 1:K) {
				p[n,t,k] = beta_k[k] - distance(to_vector(z[n]),to_vector(w[k,t]));
			}
		}
	}
	for (n in 1:N) {
		sigmasq_i[n,1] = sigmasq_k[job[n,1]];
		for (t in 2:T) {
		//eta is probability to stay, lambda is probability to stay in class but switch firm, kappa is probability to switch classes
			eta[n,t-1] = dot_product(beta',to_vector(x_norm[n,t])) + v[job[n,t-1]];
			lambda[n,t-1] = dot_product(delta',to_vector(x[n,t])) + u[job[n,t-1]];
			kappa[n,t-1] = dot_product(zeta',to_vector(x_norm[n,t])) + q[job[n,t-1]];
			sigmasq_i[n,t] = sigmasq_k[job[n,t]];
		}
	}
	sigma_i = sqrt(sigmasq_i);

}
model{
	vector[K] probs;

	sigmasq_w0 ~ scaled_inv_chi_square(1,5);
	sigmasq_w ~ scaled_inv_chi_square(1,25);
	sigmasq_bk ~ scaled_inv_chi_square(1,5);
	// shared variance for intercept for eta, lambda, kappa probabilities 
	sigmasq_dbz ~ scaled_inv_chi_square(1,13);
	// shared variance for firm coefficient for eta, lambda, kappa probabilities 
	sigmasq_uvq ~ scaled_inv_chi_square(1,3);
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
	beta_k ~ normal(0.0,sigma_bk); 
	delta ~ normal(0.00005,sigma_dbz);
	beta ~ normal(0.0001,sigma_dbz);
	zeta ~ normal(0.0, sigma_dbz);
	u ~ normal(0.0, sigma_uvq);
	v ~ normal(0.0, sigma_uvq);
	q ~ normal(0.0, sigma_uvq);
	// drawing firm at first time-point 
	for (n in 1:N) {
		job[n,1] ~ categorical_logit(to_vector(p[n,1]));
	}
	for (t in 2:T) {
		for (n in 1:N) {
			vector[3] m_prob;
			m_prob[1] = eta[n,t-1];
			m_prob[2] = lambda[n,t-1];
			m_prob[3] = kappa[n,t-1];
			m[t-1,n] ~ categorical_logit(m_prob);
			// if person leaves the class, draw new class based on latent space probabilities 
			if (m[t-1,n]==3) {
				job[n,t] ~ categorical_logit(to_vector(p[n,t]));
			}
			// if person probably stays in class, draw new class with almost all probability on the same class as the previous time point 
			else {
				probs = rep_vector(0.0000001,K);
				probs[job[n,t-1]] = 1.0 - (K-1)*0.0000001;
				job[n,t] ~ categorical_logit(probs);
			}
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
