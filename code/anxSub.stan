data{
 int<lower=1> Nanx;
 int<lower=1> Nsubj; // # subjects
 int<lower=1,upper=Nsubj> subjAnx[Nanx];
 real<lower=0> timeAnx[Nanx];
 int<lower=1,upper=3> anx[Nanx];
 int<lower=1> p;
 row_vector[p] X[Nsubj];

 int<lower=1> Nsub;
 int<lower=1,upper=Nsubj> subjSub[Nsub];
 real<lower=0> timeSub[Nsub];
 int<lower=0,upper=1> sub[Nsub];
 //int<lower=1,upper=3> anxLag[Nsub];
 //real<lower=0> lagTime[Nsub];

}

parameters{
 real b0;
 real b1;
 //real alphaAnx[Nsubj]; // random intercept
 //real gammaAnx[Nsubj]; // random time slope
 vector[p] betaAnx; //fixed effects

 ordered[2] C;

 //real<lower=0> sig_alphaAnx;
 //real<lower=0> sig_gammaAnx;
 //real mu_gammaAnx;

 vector[4] rand[Nsubj]; // random terms: 1= int anx, 2=trend anx, 3=int sub, 4= trend sub

 //real alphaSub[Nsubj]; // random intercept
 //real gammaSub[Nsubj]; // random time slope
 vector[p] betaSub; //fixed effects

 // real<lower=0> sig_alphaSub;
 // real<lower=0> sig_gammaSub;
 // real mu_gammaSub;

 corr_matrix[4] Omega;        // prior correlation
 vector<lower=0>[4] tau;      // prior scale
 vector[3] mu_prime; //

}

transformed parameters{
 vector[4] mu_rand;

 mu_rand[1]=0;
 for(i in 2:4) mu_rand[i]=mu_prime[i-1];

}

model{

 tau ~ cauchy(0, 2.5);
 Omega ~ lkj_corr(2);

 rand ~ multi_normal(mu_rand, quad_form_diag(Omega, tau));


 for(i in 1:Nanx){
  anx[i] ~ ordered_logistic(rand[subjAnx[i]][1]+rand[subjAnx[i]][2]*timeAnx[i]+X[subjAnx[i]]*betaAnx,C);
 }
 // gammaAnx ~ normal(mu_gammaAnx,sig_gammaAnx);
 // alphaAnx ~ normal(0,sig_alphaAnx);
 betaAnx ~ normal(0,1);
 //mu_gammaAnx ~ normal(0,1);

 for(i in 1:Nsub){
  sub[i] ~ bernoulli_logit(rand[subjSub[i]][3]+rand[subjSub[i]][4]*timeSub[i]+X[subjSub[i]]*betaSub);//+b0*anxLag[i]+b1*anxLag[i]*timeLag[i]);
 }
 // gammaSub ~ normal(mu_gammaSub,sig_gammaSub);
 // alphaSub ~ normal(0,sig_alphaSub);
 betaSub ~ normal(0,1);
 //mu_gammaSub ~ normal(0,1);

}
