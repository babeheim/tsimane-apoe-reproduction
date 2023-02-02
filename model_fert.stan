functions {

  matrix GP(int K, real kappa, real tau, real delta) {

    matrix[K, K] Rho;
    real KR;

    KR = K;

    for (i in 1:(K-1)) {
    for (j in (i+1):K) {
    Rho[i, j] = kappa * exp(-tau * ((j-i)^2 / KR^2));
    Rho[j, i] = Rho[i, j];
    }
    }

    for (i in 1:K) {
    Rho[i, i] = 1;
    }

    return delta * cholesky_decompose(Rho);

  }

}

data {
	
	int n;
	int y;
	int allele [n];
  int c [n, y];
  int n_allele;
  real bmi [n];

}

parameters {
	
	vector [y] mu_raw;

  real m_mu;
 	real <lower = 0, upper = 1> kappa;
 	real <lower = 0> tau;
  real <lower = 0> delta;

	vector [n_allele] b_allele;
  real b_bmi;

}

transformed parameters{

 vector [y] mu;

 mu = m_mu +  GP(y, kappa, tau, delta) * mu_raw;

}

model {

  mu_raw ~ normal(0, 1);

  m_mu ~ normal(0, 10);
  kappa ~ beta(12, 2);
  tau ~ exponential(1);
  delta ~ exponential(1);

	b_allele ~ normal(0, 1);
  b_bmi ~ normal(0, 1);

	for (i in 1:n) {

  if (bmi[i] > -99) {

  for (j in 1:y) {

    if (c[i, j] > -99) {

      c[i, j] ~ bernoulli_logit(mu[j] + b_allele[allele[i]] + b_bmi * bmi[i]);

    }
    }
    }
    }

}
