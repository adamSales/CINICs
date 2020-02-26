library(rstan)
library(tidyverse)
library(lubridate)

load('cleanedData.RData')


dat <- dat%>%filter(!is.na(anxlevel),!is.na(Date))%>%
  group_by(StudyId)%>%mutate(ntimes=n())%>%filter(ntimes>1)%>%
  arrange(Date)%>%mutate(time=as.numeric(Date-min(Date))/100)%>%
  ungroup()%>%
  mutate(id=as.numeric(as.factor(StudyId)))


stanDat <- list(
  Ndat=nrow(dat),
  Nsubj=max(dat$id),
  subj=dat$id,
  time=dat$time,
  anx=dat$anxlevel
  )

dat <- dat

dat$gender <- with(dat,ifelse(is.na(PresentSex_R),birthSex_R,PresentSex_R))
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

X <- model.matrix(~age+gender+race+hisp+Site,covDat)

stanDat$X <- scale(X[,-1])

stanDat$p <- ncol(X)-1

mod <- stan('anx2.stan',data=stanDat,chains=4,iter=1000,cores=10)
