library(tidyverse)
library(rstan)
library(reshape2)
library(MASS)


load('processedDat.RData')

options(mc.cores = min(parallel::detectCores(),10))
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')


## transform data to fit anxHard.stan
stanDat <- list()

dat$id <- as.numeric(as.factor(dat$StudyId))

stopifnot(max(dat$id)==n_distinct(dat$StudyId))
stanDat$N <- max(dat$id)

dat$time <- dat$time+1

stanDat$nt <- 12

stanDat$anxObs <- as.matrix(dat%>%select(id,time,anxlevel)%>%dcast(id~time,value.var='anxlevel')%>%select(-id))
stanDat$isObsAnx <- ifelse(is.na(stanDat$anxObs),0,1)
stanDat$anxObs[stanDat$isObsAnx==0] <- 1

stanDat$subObs <- as.matrix(dat%>%select(id,time,anyHard)%>%dcast(id~time,value.var='anyHard')%>%select(-id))
stanDat$isObsSub <- ifelse(is.na(stanDat$subObs),0,1)
stanDat$subObs[stanDat$isObsSub==0] <- 1


## fill in missing values by assumption
dat$race <- as.factor(dat$Race_R)
dat$race[is.na(dat$race)] <- 1
levels(dat$race) <- names(attr(datO$Race_R,'labels'))

dat$Hispanic_R[is.na(dat$Hispanic_R)] <- 2
dat$hisp <- ifelse(dat$Hispanic_R==1,1,0)

dat$gender[is.na(dat$gender)] <- 1

X <- model.matrix(~age+gender+race+hisp+Site,
  data=dat%>%group_by(id)%>%summarize_all(function(x) x[1]))

stanDat$X <- scale(X[,-1])

stanDat$p <- ncol(X)-1

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
