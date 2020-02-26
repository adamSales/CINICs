library(rstan)
library(tidyverse)
library(lubridate)

load('cleanedData.RData')


dat <- dat%>%filter(!is.na(anxlevel)|!is.na(anyHard)|!is.na(alcCat)|!is.na(Dep_total_G),
  !is.na(Date))%>%
  group_by(StudyId)%>%mutate(ntimes=n())%>%filter(ntimes>1)%>%
  arrange(Date)%>%mutate(time=as.numeric(Date-min(Date))/100)%>%
  ungroup()%>%
  mutate(id=as.numeric(as.factor(StudyId)),
    dep=ifelse(Dep_total_G<3,Dep_total_G,3))


datAnx <- filter(dat,!is.na(anxlevel))
datSub <- filter(dat,!is.na(anyHard))
datAlc <- filter(dat,!is.na(alcCat))
datDep <- filter(dat,!is.na(dep))

stanDat <- list(
  Nanx=nrow(datAnx),
  Nsubj=max(dat$id),
  subjAnx=datAnx$id,
  timeAnx=datAnx$time,
  anx=datAnx$anxlevel,
  Nsub=nrow(datSub),
  subjSub=datSub$id,
  timeSub=datSub$time,
  sub=datSub$anyHard,
  alc=datAlc$alcCat,
  Nalc=nrow(datAlc),
  subjAlc=datAlc$id,
  timeAlc=datAlc$time,
  dep=datDep$dep,
  Ndep=nrow(datDep),
  subjDep=datDep$id,
  timeDep=datDep$time

)

## subject-level covariates

dat$race <- as.factor(dat$Race_R)

covDat <- dat%>%
  dplyr::select(id,age,gender,race,Hispanic_R,Site,ntimes)%>%
  group_by(id)%>%
  summarize_all(function(x) if(any(!is.na(x))) x[!is.na(x)][1] else NA)

## fill in missing values by assumption
covDat$race[is.na(covDat$race)] <- 1

covDat$Hispanic_R[is.na(covDat$Hispanic_R)] <- 2
covDat$hisp <- ifelse(covDat$Hispanic_R==1,1,0)

covDat$gender[is.na(covDat$gender)] <- 1

X <- model.matrix(~age+gender+race+hisp+Site+ntimes,covDat)

stanDat$X <- scale(X[,-1])

stanDat$p <- ncol(X)-1

mod2 <- stan('anxDepSubAlc.stan',data=stanDat,chains=5,iter=5000,cores=5,pars='rand',include=FALSE); save(mod2,file='anxDepSubAlc.RData')

### correlations
omega <- rstan::extract(mod2,'Omega')[[1]]

omegaM <- apply(omega,c(2,3),mean)
omegaP <- apply(omega,c(2,3),function(x) min(mean(x<0),mean(x>0))*2)

omegaM <- omegaM[c(1,2,7,8,5,6,3,4),c(1,2,7,8,5,6,3,4)]
omegaP <- omegaP[c(1,2,7,8,5,6,3,4),c(1,2,7,8,5,6,3,4)]

colnames(omegaM) <- rownames(omegaM) <- colnames(omegaP) <- rownames(omegaP) <-
  paste(rep(c('anx','dep','alc','sub'),each=2),rep(c('int','slope'),4))

slopeR <- omegaM[c(2,4,6,8),c(2,4,6,8)]
intR <- omegaM[c(1,3,5,7),c(1,3,5,7)]
colnames(slopeR) <- rownames(slopeR) <- colnames(intR) <- rownames(intR) <- c('anx','dep','alc','sub')

intSlopeR <- omegaM[c(1,3,5,7),c(2,4,6,8)]
rownames(intSlopeR) <- paste(c('anx','dep','alc','sub'),'int')
colnames(intSlopeR) <- paste(c('anx','dep','alc','sub'),'slope')
