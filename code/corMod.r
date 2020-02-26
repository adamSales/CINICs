library(rstan)
library(tidyverse)
library(lubridate)

load('cleanedData.RData')
dat$Site[dat$Site=='NA'] <- NA

dat <- dat%>%group_by(StudyId)%>%mutate(Site=ifelse(any(!is.na(Site)), Site[!is.na(Site)][1],NA))%>%ungroup()

dat <- dat%>%filter(!is.na(anxlevel)|!is.na(anyHard)|!is.na(cig3c)|!is.na(alcCat)|!is.na(cd4)|!is.na(Dep_total_G)|
                      !is.na(qol1_5_total),
  !is.na(Site),
  !is.na(Date))%>%
  group_by(StudyId)%>%mutate(ntimes=n())%>%filter(ntimes>1)%>%
arrange(Date)%>%
mutate(time=as.numeric(Date-min(Date))/100)%>%
  ungroup()%>%
  mutate(id=as.numeric(as.factor(StudyId)),
    dep=ifelse(Dep_total_G<3,Dep_total_G,3),
    qol=ifelse(qol1_5_total>10,qol1_5_total-10,1),
    anx=anxlevel,
    hard=anyHard,
    alc=alcCat,cig=cig3c)

## subject-level covariates

dat$race <- as.factor(dat$Race_R)



covDat <- dat%>%
  dplyr::select(id,age,gender,race,Hispanic_R,Site,ntimes)%>%
  filter(!is.na(Site))%>%
  group_by(id)%>%
  summarize_all(function(x) if(any(!is.na(x))) x[!is.na(x)][1] else NA)



## fill in missing values by assumption
covDat$race[is.na(covDat$race)] <- 1

covDat$Hispanic_R[is.na(covDat$Hispanic_R)] <- 2
covDat$hisp <- ifelse(covDat$Hispanic_R==1,1,0)

covDat$gender[is.na(covDat$gender)] <- 1

covDat$age[is.na(covDat$age)] <- mean(covDat$age,na.rm=TRUE)

X <- model.matrix(~age+gender+race+hisp+Site+ntimes,covDat)


stanDat <- list(
  p=ncol(X)-1,
  Npat=max(dat$id),
  X=scale(X[,-1])
)

capitalize <- function(x) paste0(toupper(substr(x,1,1)),tolower(substr(x,2,nchar(x))))

for(yy in c('anx','hard','alc','dep','cig','qol','cd4')){
  print(yy)
  ddd <- dat[!is.na(dat[[yy]]),]
  stanDat[[paste0('N',yy)]] <- nrow(ddd)
  stanDat[[paste0('pat',capitalize(yy))]] <- ddd$id
  stanDat[[paste0('time',capitalize(yy))]] <- ddd$time
  stanDat[[yy]] <- ddd[[yy]]
}


#mod <- stan('corMod.stan',data=stanDat,iter=100)

mod2 <- stan('corMod.stan',data=stanDat,iter=5000,cores=4,pars=c('psych','sub','alphaCd4','gammaCd4','alphaQol','gammaQol','trendCig','trendAlc','trendHard','subTot'), include=FALSE,thin=5)

save(mod2,file='corMod.RData')
## cd4Dat <- filter(dat,!is.na(cd4))%>%arrange(id)

## cd4Dat$cd4 <- scale(cd4Dat$cd4)

## begin <- vapply(unique(cd4Dat$id),function(id) which(cd4Dat$id==id)[1],1)

## nobs <- cd4Dat%>%group_by(id)%>%summarize(nobs=n())











## gpDat <- with(as.data.frame(cd4Dat),
##   list(
##     Npat=length(unique(id)),
##     N=nrow(cd4Dat),
##     Nobs=nobs$nobs,
##     begin=begin,
##     time=time,
##     cd4=cd4))
## gpDat$cd4 <- cd4Dat$cd4[,1]

## gp <- stan('cd4gp.stan')



## datAnx <- filter(dat,!is.na(anxlevel))
## datSub <- filter(dat,!is.na(anyHard))


## stanDat <- list(
##   Nanx=nrow(datAnx),
##   Nsubj=max(dat$id),
##   subjAnx=datAnx$id,
##   timeAnx=datAnx$time,
##   anx=datAnx$anxlevel,
##   Nsub=nrow(datSub),
##   subjSub=datSub$id,
##   timeSub=datSub$time,
##   sub=datSub$anyHard
##   )


## stanDat$X <- scale(X[,-1])

## stanDat$p <- ncol(X)-1

## mod <- stan('anxSub.stan',data=stanDat,chains=4,iter=1000,cores=10); save(mod,file='anxSub.RData')
