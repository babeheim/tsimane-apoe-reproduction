data {
	
	int n;
	vector [n] outcome;
	int allele [n];
  real bmi [n];
  int n_allele;

}

parameters {

  real <lower = 0> sigma;

  real a;
	vector [n_allele] b_allele;
  real b_bmi;

}

model {

  sigma ~ exponential(1);

	a ~ normal(0, 1);
  b_allele ~ normal(0, 1);
  b_bmi ~ normal(0, 1);

	for (i in 1:n) {
  if (outcome[i] > -99) {
  if (bmi[i] > -99) {

    outcome[i] ~ normal(a + b_allele[allele[i]] + b_bmi * bmi[i], sigma);

  }
  }
  }

}
