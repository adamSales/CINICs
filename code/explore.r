library(lme4)
library(tidyverse)
library(ordinal)
## load('data/cleanedData.RData')

## naDate <- as.Date(NA)

## dat <- dat%>%
##   mutate(rownum=1:n())%>%
##   filter(Date>as.Date("1900-01-01"))%>%
##   arrange(Date)%>%
##   group_by(StudyId)%>%
##   mutate(lagDate=lag(Date))%>%
##   mutate_at(
##     vars(-Date,-StudyId,-lagDate,-BirthYear,-ends_with('_R'),-Site,-age,-rownum),
##     list(
##       lastObs=~na.locf(if_else(!is.na(lag(.)),lagDate,naDate),na.rm=FALSE),
##       lag=~na.locf(lag(.),na.rm=FALSE)
##     )
##   )%>%
##   select(-lagDate)%>%
##   arrange(rownum)%>%
##   select(-rownum)


## save(dat,file='cleanedDataWithLags.RData')

load('cleanedDataWithLags.RData')

dat <- dat%>%
  group_by(StudyId)%>%
  mutate_at(vars(BirthYear,ends_with('_R'),Site,age),~if(any(!is.na(.))) na.omit(.)[1] else NA)%>%
  mutate(cd4_100=cd4/100,cd4_100_lag=cd4_lag/100)%>%
  ungroup()

dat <- droplevels(dat)

## adhere as a function of lagged depression

ggplot(dat,aes(log(dep_total_lag+1),adhere))+geom_jitter()+geom_smooth(method='lm')

dat%>%group_by(adhere,dep_total_lag)%>%mutate(n=n())%>%
  ggplot(aes(log(dep_total_lag+1),adhere))+geom_point(aes(size=n),show.legend=FALSE)+
  geom_smooth(method='lm',method.args=list(weights=n))+scale_size(range=c(.5,100))

 dat%>%
   mutate(depLagCat=paste('Lagged Depression=',ifelse(dep_total_lag==0,0,ifelse(dep_total_lag<4,'1-3','4+'))))%>%
   ggplot(aes(adhere))+geom_bar(aes(y = (..count..)/sum(..count..)))+coord_flip()+facet_wrap(~depLagCat,nrow=1)+
   scale_y_continuous('Percentage',labels=scales::percent)

dat%>%group_by(dep_total_lag)%>%summarize(n=n(),avg.adhere=mean(adhere,na.rm=TRUE))%>%
  ggplot(aes(log(dep_total_lag+1),avg.adhere,size=n))+geom_point(show.legend=FALSE)+
  geom_smooth(method='lm',show.legend=FALSE)+scale_size(range=c(.5,20))
ggsave('laggedDepressionVsMeanAdhere.jpg')



dat%>%group_by(dep_total_lag)%>%summarize(n=n(),avg.adhere=mean(adhere,na.rm=TRUE))%>%
  ggplot(aes(dep_total_lag,avg.adhere,size=n))+geom_point(show.legend=FALSE)+
  geom_smooth(method='lm',show.legend=FALSE)+scale_size(range=c(.5,20))
ggsave('laggedDepressionVsMeanAdhereNoLog.jpg')


dat$adDepDist <- (dat$Date-dat$dep_total_lastObs)/1000

#aaa <- lmer(adhere~log(dep_total_lag+1)+adDepDist+I(as.numeric(adDepDist)^2)+as.factor(adhere_lag)+(adDepDist+I(as.numeric(adDepDist)^2)|StudyId),data=dat)
AdDep2 <- update(AdDep,.~.+age+Site+PresentSex_R+Hispanic_R+as.factor(Race_R))
#AdDep0 <- lmer(adhere~log(dep_total_lag+1)+adDepDist+as.factor(adhere_lag)+(1|StudyId),data=dat)

#AdDepInt <- lmer(adhere~dep_total_lag*as.factor(adhere_lag)*adDepDist+(adDepDist|StudyId),data=dat)

AdDepNoLog <- lmer(adhere~dep_total_lag+adDepDist+as.factor(adhere_lag)+age+Site+birthSex_R+Hispanic_R+Race_R+(adDepDist|StudyId),data=dat)

AdDepOrd <- clmm(factor(adhere,ordered=TRUE)~log(dep_total_lag+1)+adDepDist+as.factor(adhere_lag)+(adDepDist|StudyId),data=dat)
save(AdDepOrd,AdDep,AdDep2,file='fittedModels/AdDepOrd.RData')
AdDepOrd0 <- clmm(factor(adhere,ordered=TRUE)~log(dep_total_lag+1)+adDepDist+as.factor(adhere_lag)+(1|StudyId),data=dat)
anova(AdDepOrd,AdDepOrd0)
AdDepOrd2 <- clmm(factor(adhere,ordered=TRUE)~log(dep_total_lag+1)+adDepDist+as.factor(adhere_lag)+(adDepDist|StudyId)+age+Site+birthSex_R+Hispanic_R+Transgendered_R+Race_R,data=dat)
save(AdDepOrd2,file='fittedModels/AdDepOrd2.RData')

nd <- within(dat,{adDepDist=0;adhere=3})[rownames(model.frame(AdDep)),]
ppp <- predict(AdDep,nd)
nd%>%
  mutate(adhere=ppp+resid(AdDep))%>%
  group_by(dep_total_lag)%>%summarize(n=n(),avg.adhere=mean(adhere,na.rm=TRUE))%>%
  ggplot(aes(log(dep_total_lag+1),avg.adhere,size=n))+geom_point(show.legend=FALSE)+
  geom_smooth(method='lm',show.legend=FALSE)+scale_size(range=c(.5,5))+ylim(1,5)
ggsave('laggedDepressionVsMeanAdhereAdj.jpg')


## cd4 count as function of depression and adhere

p <- dat%>%
  mutate(depLagCat=paste('Lagged Depression=',ifelse(dep_total_lag==0,0,ifelse(dep_total_lag<4,'1-3','4+'))))%>%
  ggplot(aes(factor(adhere_lag,exclude=NA),cd4_100))+
  geom_boxplot()+
  geom_smooth(aes(adhere_lag,cd4_100),method='lm')+
  labs(x='Lagged Adherence',y='CD4 Count/100')+
  scale_x_discrete(na.translate=FALSE)
ggsave('adhereLagVsCD4.jpg',p)

ggsave('adhereLagVsCD4byDepLag.jpg',p+facet_wrap(~depLagCat,nrow=1))

dat$depDist <- (dat$Date-dat$dep_total_lastObs)/1000
dat$adhereDist <- (dat$Date-dat$adhere_lastObs)/1000
dat$cd4Dist <- (dat$Date-dat$cd4_lastObs)/1000

cd4AdDep <- lmer(cd4_100~log(dep_total_lag+1)+depDist+adhere_lag+adhereDist+cd4_100_lag+(depDist+adhereDist+cd4Dist|StudyId),data=dat)
summary(update(cd4AdDep,subset=cd4Dist>=adhereDist))
cd4AdDep2 <- update(cd4AdDep,.~.+age+Site+PresentSex_R+Hispanic_R+as.factor(Race_R)



save(cd4AdDep,cd4AdDep2,file='fittedModels/cd4AdDep.RData')


mmm <- colMeans(model.matrix(cd4AdDep))

nd <- within(dat,{dep_total_lag=exp(1.06)-1;depDist=0;adhereDist=0;cd4_100_lag=5.3})[rownames(model.frame(cd4AdDep)),]

ppp <- predict(cd4AdDep,nd)

p <- nd%>%
  mutate(cd4_100=ppp+resid(cd4AdDep))%>%
  ggplot(aes(factor(adhere_lag,exclude=NA),cd4_100))+
  geom_boxplot()+
  geom_smooth(aes(adhere_lag,cd4_100),method='lm')+
  labs(x='Lagged Adherence',y='CD4 Count/100')+
  scale_x_discrete(na.translate=FALSE)+ylim(1,5)
#ggsave('adhereLagVsCD4.jpg',p)


## qol as function of depression, adhere, cd4

dat%>%
  mutate(
    #depLagCat=paste('Lagged Depression=',ifelse(dep_total_lag==0,0,ifelse(dep_total_lag<4,'1-3','4+'))),
    cd4LagCat=round(as.numeric(as.character(cut(dat$cd4_lag,quantile(dat$cd4_lag,seq(0,1,.05),na.rm=TRUE),labels=quantile(dat$cd4_lag,seq(0.05,1,.05)-.025,na.rm=TRUE)))),-1)
    #adhere_lag=paste('Lagged Adhere=\n',adhere_lag)
  )%>%
  filter(!is.na(dep_total_lag),!is.na(cd4_lag),!is.na(adhere_lag))%>%
  gather("predictor","value",dep_total_lag,cd4LagCat,adhere_lag)%>%
  group_by(predictor,value)%>%
  summarize(n=n(),avg.qol=mean(qol1_5_total,na.rm=TRUE))%>%
  ungroup()%>%
  ggplot(aes(value,avg.qol))+
  geom_point(aes(size=n),show.legend=FALSE)+
  geom_smooth(method='lm',formula=y~poly(x,2))+
  facet_wrap(~predictor,scales="free_x")

ggsave('qolByDepCd4Adhere.jpg')

dat$cd4Dist <- (dat$Date-dat$cd4_lastObs)/1000

qolAdDepCd4 <- clmm(factor(qol1_5_total,ordered=TRUE)~log(dep_total_lag+1)+depDist+adhere_lag+adhereDist+cd4_100_lag+cd4Dist+qol1_5_total_lag+(depDist+adhereDist+cd4Dist|StudyId),data=dat)


qolAdDepCd4int <- clmm(factor(qol1_5_total,ordered=TRUE)~log(dep_total_lag+1)*cd4_100_lag*adhere_lag+
                      depDist+adhereDist+cd4Dist+qol1_5_total_lag+(depDist+adhereDist+cd4Dist|StudyId),data=dat)

qolAdDepCd4.1 <- update(qolAdDepCd4,.~.+age+Site+birthSex_R+Hispanic_R+Transgendered_R+Race_R)

qolAdDepCd4int.1 <- update(qolAdDepCd4int,.~.+age+Site+birthSex_R+Hispanic_R+Transgendered_R+Race_R)


save(list=grep('qolAdDepCd4',ls(),value=TRUE),file='fittedModels/qolAdDepCd4.RData')
  ## group_by(dep_total_lag,adhere_lag,cd4LagCat)%>%
  ## summarize(n=n(),avg.qol=mean(qol1_5_total,na.rm=TRUE),cd4Lag=mean(cd4,na.rm=TRUE))%>%


cat('Increasing log depression by 1 sd decreases adherence by 0.06 sd')
cat('Increasing adherence by 1 (=1SD) increases cd4 count by 8')
# increasing cd4 by 100 increases odds of moving to higher qol level by 1%
# increasing adherences by 1 increases odds "" by 10%

