library(tidyverse)
library(reshape2)
library(MASS)
library(R2jags)

load('processedDat.RData')




## transform data to fit anxHard.stan

jagsDat <- list()
dat$id <- as.numeric(as.factor(dat$StudyId))

stopifnot(max(dat$id)==n_distinct(dat$StudyId))
jagsDat$N <- max(dat$id)

dat$time <- dat$time+1

jagsDat$nt <- 12

jagsDat$anx <- as.matrix(dat%>%dplyr::select(id,time,anxlevel)%>%dcast(id~time,value.var='anxlevel')%>%dplyr::select(-id))


## fill in missing values by assumption
dat$race <- as.factor(dat$Race_R)
dat$race[is.na(dat$race)] <- 1
levels(dat$race) <- names(attr(datO$Race_R,'labels'))

dat$Hispanic_R[is.na(dat$Hispanic_R)] <- 2
dat$hisp <- ifelse(dat$Hispanic_R==1,1,0)

dat$gender[is.na(dat$gender)] <- 1

X <- model.matrix(~age+gender+race+hisp+Site,
  data=dat%>%group_by(id)%>%summarize_all(function(x) x[1]))

X <- scale(X[,-1])

jagsDat$X <- X

jagsDat$ncov <- ncol(X)

mod <- jags(jagsDat,parameters=c('betaAnx','C','gamAnx'),model='anxHard.bug',n.chains=1,n.iter=100)


#mod <- stan('anxHard2.stan',data=stanDat,chains=1,iter=10,pars=c('anx','sub'),include=FALSE)


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


mod <- stan('anxHard2.stan',data=stanDat,chains=4,iter=10000,pars=c('b','cAnx','gamAnx','gamSub','betaAnx','betaSub','alphaSub'),init=initFun,cores=4)
save(mod,file='anxHardMod.RData')
