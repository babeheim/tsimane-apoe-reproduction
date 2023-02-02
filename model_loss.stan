data {
	
	int n;
	int loss [n];
	int allele [n];
  real bmi [n];

}

parameters {

  real a;
	vector [2] b_allele;
  real b_bmi;

}

model {

	a ~ normal(0, 1);
  b_allele ~ normal(0, 1);
  b_bmi ~ normal(0, 1);

	for (i in 1:n) {
  if (loss[i] > -99) {
  if (bmi[i] > -99) {

    loss[i] ~ bernoulli_logit(a + b_allele[allele[i]] + b_bmi * bmi[i]);

  }
  }
  }

}
