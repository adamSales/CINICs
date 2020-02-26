attach(stanDat)
cAnx <- c(-.5,.5)
b <- 0
gamAnx <- 0
gamSub <- 0
betaAnx <- rep(0,p)
betaSub <- rep(0,p)
alphaAnx <- 0
alphaSub <- 0
anx <- matrix(0,N,nt)
sub <- matrix(0,N,nt)


anxInit <- function(isObsAnx,anxObs,cAnx){
  nt <- ncol(anxObs)
  N <- nrow(anxObs)
  anx <- matrix(rnorm(nt*N,0,1),N,nt)
  anx[isObsAnx==1&anxObs==1] <- cAnx[1]-.2
  anx[isObsAnx==1&anxObs==2] <- mean(cAnx)
  anx[isObsAnx==1&anxObs==3] <- cAnx[2]+.2
  anx
}

subInit <- function(isObsSub,subObs){
  nt <- ncol(subObs)
  N <- nrow(subObs)
  sub <- matrix(rnorm(nt*N,0,1),N,nt)
  sub[isObsSub==1&subObs==0] <- -.5
  sub[isObsSub==1&subObs==1] <- .5
  sub
}

initFun <- function(){
  attach(stanDat)
  on.exit(detach(stanDat))
  anxMod <- polr(ordered(anxObs[,1])~X,method='probit')
  subMod <- glm(subObs[,1]~X,family=binomial('probit'))
  list(
    b=rnorm(1,0,.5),
    gamAnx=rnorm(1,.25,.5),
    gamSub=rnorm(1,.25,.5),
    betaAnx=coef(anxMod),
    betaSub=coef(subMod)[-1],
    #alphaAnx=rnorm(1),
    alphaSub=coef(subMod)[1],
    cAnx=anxMod$zeta,
    anx=anxInit(isObsAnx,anxObs,anxMod$zeta),
    sub=subInit(isObsSub,subObs))
}

iii <- initFun()

attach(iii)

xbAnx <- X%*%betaAnx
xbSub <- X%*%betaSub

etaAnx <- matrix(NA,N,nt)
etaSub <- matrix(NA,N,nt)

  etaAnx[,1]=alphaAnx+xbAnx
etaSub[,1]=alphaSub+xbSub

for(t in 2:nt){
  etaAnx[,t]=alphaAnx+xbAnx+gamAnx*anx[,t-1];
  etaSub[,t]=alphaSub+xbSub+gamSub*sub[,t-1]+b*anx[,t-1];
}
prob <- matrix(1,N,nt)
for(i in 1:N) for(t in 1:12) if(isObsAnx[i,t]==1) prob[i,t] <- min(pnorm(anx[i,t],etaAnx[i,t],1),pnorm(anx[i,t],etaAnx[i,t],1,lower.tail=FALSE))
min(prob)
