data{
 int<lower=1> Npat;
 int<lower=1> Nobs[Npat];
 int<lower=1> N;
 int<lower=1> begin[Npat];
 real<lower=0> time[N];
 real cd4[N];
}
transformed data {
  real delta = 1e-9;
}
parameters{
 real<lower=0> rho[Npat];
 real<lower=0> alpha[Npat];
 real<lower=0> sigma[Npat];
 vector[N] eta;
 real<lower=0> sigInt;
}

model {
 for(i in 1:Npat){
  vector[Nobs[i]] f;
  {
    real timei[Nobs[i]]= segment(time,begin[i],begin[i+1]);
    matrix[Nobs[i], Nobs[i]] L_K;
    matrix[Nobs[i], Nobs[i]] K = cov_exp_quad(timei, alpha[i], rho[i]);

    // diagonal elements
    for (n in 1:Nobs[i])
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
  for(k in begin[i]:(begin[i+1]-1))
   cd4[k] ~ normal(intercept[i]+f[k-begin[i]+1], sigma[i]);
 }
 intercept~normal(0,sigInt);
 sigInt~cauchy(0, 2.5);
 rho ~ inv_gamma(5, 5);
 alpha ~ std_normal();
 sigma ~ std_normal();
 eta ~ std_normal();


}