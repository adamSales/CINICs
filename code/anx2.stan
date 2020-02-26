data{
 int<lower=1> Ndat;
 int<lower=1> Nsubj; // # subjects
 int<lower=1,upper=Nsubj> subj[Ndat];
 real<lower=0> time[Ndat];
 int<lower=1,upper=3> anx[Ndat];
 int<lower=1> p;
 row_vector[p] X[Nsubj];

}

parameters{
 real alpha[Nsubj]; // random intercept
 real gamma[Nsubj]; // random time slope
 vector[p] beta; //fixed effects

 ordered[2] C;

 real<lower=0> sig_alpha;
 real<lower=0> sig_gamma;
 real mu_gamma;

}

model{
 for(i in 1:Ndat){
  anx[i] ~ ordered_logistic(alpha[subj[i]]+gamma[subj[i]]*time[i]+X[subj[i]]*beta,C);
 }
 gamma ~ normal(mu_gamma,sig_gamma);
 alpha ~ normal(0,sig_alpha);
 beta ~ normal(0,1);
 mu_gamma ~ normal(0,1);
}
