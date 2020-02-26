// preliminary--can I get this to work?
// regress hard substance abuse ("anyHard") on anxiety level

data{
 int<lower=1> N; // # subjects
 int<lower=3> T; // # time-points
 int<lower=1> nobsAnx; // # observed anxiety scores
 int<lower=1,upper=3> anxObs[nobsAnx]; // observed anxiety scores, vectorized
 int<lower=1,upper=N> subjAnx[nobsAnx]; // subject ids for anxiety
 int<lower=1,upper=T> timeAnx[nobsAnx]; // time ids for anxiety
 int<lower=1> nobsSub; // # observed substance use
 int<lower=0,upper=1> subObs[nobsSub]; // observed substance use, vectorized
 int<lower=1,upper=N> subjSub[nobsSub]; // subject ids for substance use
 int<lower=1,upper=T> timeSub[nobsSub]; // time ids for substance use
 int<lower=1> p; // # stable covariates
 matrix[N,p] X; // covariate matrix

}

parameters{


 ordered[2] cAnx; // cutoffs for anxiety


 real b; // the "effect" of anxiety on substance abuse

 real gamAnx; // the autoregressive parameter for anxiety
 real gamSub; // " " for substance use
 vector[p] betaAnx; // coefficients for X
 vector[p] betaSub;
 real alphaAnx;
 real alphaSub;

 matrix[N,T] anx; // latent anxiety
 matrix[N,T] sub; // latent substance use


}

transformed parameters{
 vector[N] xbAnx;
 vector[N] xbSub;

 xbAnx=X*betaAnx;
 xbSub=X*betaSub;
}

//transformed parameters{
model{

 matrix[N,T] etaAnx;
 matrix[N,T] etaSub;

 real anxStar[nobsAnx];
 real subStar[nobsSub];

 vector[3] anxTheta[nobsAnx];

 real eta[nobsAnx];


// first time point. what assumptions am i encoding here??




 for(i in 1:N){
  etaAnx[i,1]=alphaAnx+xbAnx[i];
  etaSub[i,1]=alphaSub+xbSub[i];
  anx[i,1]~normal(etaAnx[i,1],1);
  sub[i,1]~normal(etaSub[i,1],1);
 }

 for(t in 2:T){
  etaAnx[,t]=alphaAnx+xbAnx+gamAnx*anx[,t-1];
  etaSub[,t]=alphaSub+xbSub+gamSub*sub[,t-1]+b*anx[,t-1];
  for(i in 1:N){
   anx[i,t]~normal(etaAnx[i,t],1);
   sub[,t]~normal(etaSub[i,t],1);
  }
 }

 for(i in 1:nobsAnx){
  eta[i]=etaAnx[subjAnx[i],timeAnx[i]];
  anxTheta[i][1]=1-Phi(eta[i]-cAnx[1]);
  anxTheta[i][2]=Phi(eta[i]-cAnx[1])-Phi(eta[i]-cAnx[2]);
  anxTheta[i][3]=Phi(eta[i]-cAnx[2]);
 }

 for(i in 1:nobsSub)
  subStar[i]=etaSub[subjSub[i],timeSub[i]];


 alphaAnx~normal(0,1);
 alphaSub~normal(0,1);
 b~normal(0,.5);
 gamAnx~normal(.25,.5);
 gamSub~normal(.25,.5);
 betaAnx~normal(0,1);
 betaSub~normal(0,1);

// now actual observed data
 subObs~bernoulli(Phi(subStar));

 for(i in 1:nobsAnx)
  anxObs[i]~categorical(anxTheta[i]);
}

