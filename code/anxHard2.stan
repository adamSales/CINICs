// preliminary--can I get this to work?
// regress hard substance abuse ("anyHard") on anxiety level

data{
 int<lower=1> N; // # subjects
 int<lower=3> nt; // # time-points
 int<lower=0,upper=1> isObsAnx[N,nt];
 int<lower=0,upper=1> isObsSub[N,nt];
 int<lower=1,upper=3> anxObs[N,nt]; // observed anxiety scores, vectorized
 int<lower=0,upper=1> subObs[N,nt]; // observed substance use, vectorized
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
// real alphaAnx;
 real alphaSub;

 matrix[N,nt] anx; // latent anxiety
 matrix[N,nt] sub; // latent substance use


}

transformed parameters{
 vector[N] xbAnx;
 vector[N] xbSub;

 xbAnx=X*betaAnx;
 xbSub=X*betaSub;
}

//transformed parameters{
model{

 matrix[N,nt] etaAnx;
 matrix[N,nt] etaSub;



// priors
// alphaAnx~normal(0,1);
 alphaSub~normal(0,1);
 b~normal(0,.5);
 gamAnx~normal(.25,.5);
 gamSub~normal(.25,.5);
 betaAnx~normal(0,1);
 betaSub~normal(0,1);


// first time point. what assumptions am i encoding here??
  etaAnx[,1]=xbAnx; //+alphaAnx;
  etaSub[,1]=alphaSub+xbSub;


 for(t in 2:nt){
  etaAnx[,t]=xbAnx+gamAnx*anx[,t-1]; //+alphaAnx;
  etaSub[,t]=alphaSub+xbSub+gamSub*sub[,t-1]+b*anx[,t-1];
 }

 for(i in 1:N){
  for(t in 1:nt){
   if(isObsAnx[i,t]==1){
    if(anxObs[i,t]==1){
     anx[i,t]~normal(etaAnx[i,t],1)T[,cAnx[1]];
    } else if(anxObs[i,t]==2){
     anx[i,t]~normal(etaAnx[i,t],1)T[cAnx[1],cAnx[2]];
    } else{
     anx[i,t]~normal(etaAnx[i,t],1)T[cAnx[2],];
    }
    // if(anxObs[i,t]>1){
    //  anx[i,t]~normal(etaAnx[i,t],1)T[0,];
    // } else{
    //  anx[i,t]~normal(etaAnx[i,t],1)T[,0];
    // }
   } else{
    anx[i,t]~normal(etaAnx[i,t],1);
   }
   if(isObsSub[i,t]==1){
    if(subObs[i,t]==1){
     sub[i,t]~normal(etaSub[i,t],1)T[0,];
    } else{
     sub[i,t]~normal(etaSub[i,t],1)T[,0];
    }
   } else{
    sub[i,t]~normal(etaSub[i,t],1);
   }
  }
 }
}

