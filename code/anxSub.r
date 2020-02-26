library(rstan)
library(tidyverse)
library(lubridate)

load('cleanedData.RData')


dat <- dat%>%filter(!is.na(anxlevel)|!is.na(anyHard),
  !is.na(Date))%>%
  group_by(StudyId)%>%mutate(ntimes=n())%>%filter(ntimes>1)%>%
  arrange(Date)%>%mutate(time=as.numeric(Date-min(Date))/100)%>%
  ungroup()%>%
  mutate(id=as.numeric(as.factor(StudyId)))


datAnx <- filter(dat,!is.na(anxlevel))
datSub <- filter(dat,!is.na(anyHard))


stanDat <- list(
  Nanx=nrow(datAnx),
  Nsubj=max(dat$id),
  subjAnx=datAnx$id,
  timeAnx=datAnx$time,
  anx=datAnx$anxlevel,
  Nsub=nrow(datSub),
  subjSub=datSub$id,
  timeSub=datSub$time,
  sub=datSub$anyHard
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

X <- model.matrix(~age+gender+race+hisp+Site,covDat)

stanDat$X <- scale(X[,-1])

stanDat$p <- ncol(X)-1

mod <- stan('anxSub.stan',data=stanDat,chains=4,iter=5000,cores=4); save(mod,file='anxSub.RData')
