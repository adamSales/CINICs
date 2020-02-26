library(lme4)
library(tidyverse)
library(ordinal)
library(texreg)
library(gridExtra)


load('data/cleanedDataWithLags.RData')

dat <- dat%>%
  group_by(StudyId)%>%
  mutate_at(vars(BirthYear,ends_with('_R'),Site,age),~if(any(!is.na(.))) na.omit(.)[1] else NA)%>%
  mutate(cd4_100=cd4/100,cd4_100_lag=cd4_lag/100,
    depDist=(Date-dep_total_lastObs)/1000,
    adDist=(Date-adhere_lastObs)/1000,
    cd4Dist=(Date-cd4_lastObs)/1000,
    qolDist=(Date-qol1_5_total_lastObs)/1000,
    Race_R=as.factor(Race_R),
    logDepLag=log(dep_total_lag+1)
  )%>%
  ungroup()

dat <- droplevels(dat)



##################################################
### adherence vs depression
##################################################

###### models
AdDep00 <- lmer(adhere~logDepLag+(1|StudyId),data=dat)
AdDep0 <- lmer(adhere~logDepLag+as.factor(adhere_lag)+(1|StudyId),data=dat)
AdDep <- lmer(adhere~logDepLag+depDist+adDist+as.factor(adhere_lag)+(adDist+depDist|StudyId),data=dat)
AdDep2 <- update(AdDep,.~.+age+Site+PresentSex_R+Hispanic_R+Race_R)
AdDep3 <- update(AdDep0,subset=adDist>=depDist)

save(AdDep00,AdDep0,AdDep,AdDep2,AdDep3,file='fittedModels/AdDepModels.RData')
htmlreg(list(AdDep00,AdDep0,AdDep,AdDep2,AdDep3),file='coefTables/adhereVsDepression.html')

###### plots

dat%>%group_by(dep_total_lag)%>%summarize(n=n(),avg.adhere=mean(adhere,na.rm=TRUE))%>%
  ggplot(aes(log(dep_total_lag+1),avg.adhere,size=n))+geom_point(show.legend=FALSE)+
  geom_smooth(method='lm',show.legend=FALSE)+scale_size(range=c(.5,20))+ylim(1,5)+
  xlab('log(Depression+1) (lagged)')+ylab('Avg. Adherence')
ggsave('plots/adhereVsDepression.jpg')

nd <- within(dat,{depDist=0;adhere_lag=3;adDist=0})[rownames(model.frame(AdDep)),]
ppp <- predict(AdDep,nd)
nd%>%
  mutate(adhere=ppp+resid(AdDep))%>%
  group_by(dep_total_lag)%>%summarize(n=n(),avg.adhere=mean(adhere,na.rm=TRUE))%>%
  ggplot(aes(log(dep_total_lag+1),avg.adhere,size=n))+geom_point(show.legend=FALSE)+
  geom_smooth(method='lm',show.legend=FALSE)+scale_size(range=c(.5,5))+ylim(1,5)+
  xlab('log(Depression+1) (lagged)')+ylab('Avg. Adherence (Adjusted)')
ggsave('plots/adhereVsDepressionAdj.jpg')


##################################################
### cd4 vs adherence & depression
##################################################

##### Models
cd4AdDep00 <- lmer(cd4~logDepLag+adhere_lag+(1|StudyId),data=dat)
cd4AdDep0 <- lmer(cd4~logDepLag+adhere_lag+cd4_lag+(1|StudyId),data=dat)
cd4AdDep <- lmer(cd4~logDepLag+depDist+adhere_lag+adDist+cd4_lag+cd4Dist+(depDist+adDist+cd4Dist|StudyId),data=dat)
cd4AdDep2 <- update(cd4AdDep,.~.+age+Site+PresentSex_R+Hispanic_R+Race_R)
cd4AdDep3 <- update(cd4AdDep0,subset=cd4Dist>=adDist&cd4Dist>=depDist)

save(list=grep('cd4AdDep',ls(),value=TRUE),file='fittedModels/cd4AdDepModels.RData')
htmlreg(list(cd4AdDep00,cd4AdDep0,cd4AdDep,cd4AdDep2,cd4AdDep3),file='coefTables/cd4VsadhereAndDepression.html')

##### Plots
dat%>%
  select(cd4,'vs. Adherence (lagged)'=adhere_lag,'vs. Depression (lagged)'=dep_total_lag)%>%
  gather("iv","value",-cd4)%>%
  mutate(bpVal=as.factor(value))%>%
  ggplot(aes(bpVal,cd4))+
  geom_boxplot()+
  geom_smooth(aes(value,cd4),method='lm')+
  facet_wrap(~iv,scale="free_x")+
  labs(x=NULL,y='CD4 Count')+
  scale_x_discrete(na.translate=FALSE)
ggsave('plots/cd4VsadhereAndDepression.jpg')


nd1 <- within(dat,{logDepLag=1.06;depDist=0;adhereDist=0;cd4_lag=531;cd4Dist=0})[rownames(model.frame(cd4AdDep)),]
ppp1 <- predict(cd4AdDep,nd1)
nd2 <- within(dat,{adhere_lag=mean(model.frame(cd4AdDep)$adhere_lag);depDist=0;adhereDist=0;cd4_lag=531;cd4Dist=0})[rownames(model.frame(cd4AdDep)),]
ppp2 <- predict(cd4AdDep,nd2)

dat[rownames(model.frame(cd4AdDep)),]%>%
  select(cd4,'vs. Adherence (lagged)'=adhere_lag,'vs. Depression (lagged)'=dep_total_lag)%>%
  gather("iv","value",-cd4)%>%
  mutate(bpVal=as.factor(value),cd4=ifelse(iv=='vs. Adherence (lagged)',ppp1,ppp2)+resid(cd4AdDep))%>%
  ggplot(aes(bpVal,cd4))+
  geom_boxplot()+
  geom_smooth(aes(value,cd4),method='lm')+
  facet_wrap(~iv,scale="free_x")+
  labs(x=NULL,y='CD4 Count (Adjusted)')+
  scale_x_discrete(na.translate=FALSE)
ggsave('plots/cd4VsadhereAndDepressionAdj.jpg')


##################################################
### qol vs cd4 & adherence & depression
##################################################

##### Models

qolAdDepCd400 <- clmm(factor(qol1_5_total,ordered=TRUE)~log(dep_total_lag+1)+adhere_lag+cd4_100_lag+(1|StudyId),data=dat)
qolAdDepCd40 <- clmm(factor(qol1_5_total,ordered=TRUE)~log(dep_total_lag+1)+adhere_lag+cd4_100_lag+as.factor(qol1_5_total_lag)+(1|StudyId),data=dat)
qolAdDepCd4 <- clmm(factor(qol1_5_total,ordered=TRUE)~log(dep_total_lag+1)+depDist+adhere_lag+adDist+cd4_100_lag+cd4Dist+as.factor(qol1_5_total_lag)+qolDist+(depDist+adDist+cd4Dist+qolDist|StudyId),data=dat)
qolAdDepCd4quad <- clmm(factor(qol1_5_total,ordered=TRUE)~poly(log(dep_total_lag+1),2)+depDist+poly(adhere_lag,2)+adDist+poly(cd4_100_lag,2)+cd4Dist+as.factor(qol1_5_total_lag)+qolDist+(depDist+adDist+cd4Dist+qolDist|StudyId),data=dat[rownames(model.frame(qolAdDepCd4)),])

qolAdDepCd42 <- update(qolAdDepCd4,.~.+age+Site+PresentSex_R+Hispanic_R+Race_R)
qolAdDepCd43 <- update(qolAdDepCd40,subset=qolDist>=adDist&qolDist>=depDist&qolDist>=cd4Dist)

save(list=grep('qolAdDepCd4',ls(),value=TRUE),file='fittedModels/qolAdDepCd4Models.RData')
htmlreg(list(qolAdDepCd400,qolAdDepCd40,qolAdDepCd4,qolAdDepCd4quad,qolAdDepCd42,qolAdDepCd43),file='coefTables/qolVscd4AndadhereAndDepression.html')

dat%>%
  mutate(
    cd4LagCat=round(as.numeric(as.character(cut(dat$cd4_lag,quantile(dat$cd4_lag,seq(0,1,.05),na.rm=TRUE),labels=quantile(dat$cd4_lag,seq(0.05,1,.05)-.025,na.rm=TRUE)))),-1)
  )%>%
  filter(!is.na(dep_total_lag),!is.na(cd4_lag),!is.na(adhere_lag))%>%
  select(qol1_5_total,'vs. Adherence (lagged)'=adhere_lag,'vs. Depression (lagged)'=dep_total_lag,'vs CD4 Count (lagged)'=cd4LagCat)%>%
  gather("predictor","value",-qol1_5_total)%>%
  group_by(predictor,value)%>%
  summarize(n=n(),avg.qol=mean(qol1_5_total,na.rm=TRUE))%>%
  ungroup()%>%
  ggplot(aes(value,avg.qol))+
  geom_point(aes(size=n),show.legend=FALSE)+
  geom_smooth(method='lm',formula=y~poly(x,2))+
  facet_wrap(~predictor,scales="free_x")+
  labs(x=NULL,y='Avg. QOL')+
ggsave('plots/qolByDepCd4Adhere.jpg')

simpQolMod <- lmer(qol1_5_total~poly(logDepLag,2)+poly(adhere_lag,2)+poly(cd4_100_lag,2)+qol1_5_total_lag+(1|StudyId),data=dat[rownames(model.frame(qolAdDepCd40)),])
ndAd <- within(dat[rownames(model.frame(qolAdDepCd40)),],{logDepLag=1.06;cd4_100_lag=5.31;qol1_5_total_lag=13})#[rownames(model.frame(simpQolMod)),]
pppAd <- predict(simpQolMod,ndAd)
ndDep <- within(dat[rownames(model.frame(qolAdDepCd40)),],{adhere_lag=mean(model.frame(cd4AdDep)$adhere_lag);cd4_100_lag=5.31;qol1_5_total_lag=13})#[rownames(model.frame(simpQolMod)),]
pppDep <- predict(simpQolMod,ndDep)
ndCd4 <- within(dat[rownames(model.frame(qolAdDepCd40)),],{adhere_lag=mean(model.frame(cd4AdDep)$adhere_lag);logDepLag=1.06;qol1_5_total_lag=13})#[rownames(model.frame(simpQolMod)),]
pppCd4 <- predict(simpQolMod,ndCd4)

dat[rownames(model.frame(qolAdDepCd40)),]%>%
  mutate(
    cd4LagCat=round(as.numeric(as.character(cut(cd4_lag,quantile(cd4_lag,seq(0,1,.05),na.rm=TRUE),labels=quantile(cd4_lag,seq(0.05,1,.05)-.025,na.rm=TRUE)))),-1)
  )%>%
  select(qol1_5_total,'vs. Adherence (lagged)'=adhere_lag,'vs. Depression (lagged)'=dep_total_lag,'vs CD4 Count (lagged)'=cd4LagCat)%>%
  gather("predictor","value",-qol1_5_total)%>%
  mutate(qol=ifelse(predictor=='vs. Adherence (lagged)',pppAd,ifelse(predictor=='vs. Depression (lagged)',pppDep,pppCd4))+resid(simpQolMod))%>%
  group_by(predictor,value)%>%
  summarize(n=n(),avg.qol=mean(qol,na.rm=TRUE))%>%
  ungroup()%>%
  ggplot(aes(value,avg.qol))+
  geom_point(aes(size=n),show.legend=FALSE)+
  geom_smooth(method='lm',formula=y~poly(x,2))+
  facet_wrap(~predictor,scales="free_x")+
  labs(x=NULL,y='Avg. QOL (adjusted)')+
ggsave('plots/qolByDepCd4AdhereAdj.jpg')

dat2 <- dat%>%mutate(logDepLag=logDepLag-mean(logDepLag,na.rm=TRUE),
  adhere_lag=adhere_lag-mean(adhere_lag,na.rm=TRUE),
  cd4_100_lag=cd4_100_lag-mean(cd4_100_lag,na.rm=TRUE))

qolAdDepCd400l <- lmer(qol1_5_total~logDepLag+adhere_lag+cd4_100_lag+(1|StudyId),data=dat2)
qolAdDepCd40l <- lmer(qol1_5_total~logDepLag+adhere_lag+cd4_100_lag+as.factor(qol1_5_total_lag)+(1|StudyId),data=dat2)
qolAdDepCd4l <- lmer(qol1_5_total~logDepLag+depDist+adhere_lag+adDist+cd4_100_lag+cd4Dist+as.factor(qol1_5_total_lag)+qolDist+(cd4Dist+adDist|StudyId),data=dat2)
qolAdDepCd4quadl <- update(qolAdDepCd4l,.~.+I(logDepLag^2)+I(adhere_lag^2)+I(cd4_100_lag^2))

qolAdDepCd42l <- update(qolAdDepCd4l,.~.+age+Site+PresentSex_R+Hispanic_R+Race_R)
qolAdDepCd43l <- update(qolAdDepCd40l,subset=qolDist>=adDist&qolDist>=depDist&qolDist>=cd4Dist)

htmlreg(list(qolAdDepCd400l,qolAdDepCd40l,qolAdDepCd4l,qolAdDepCd4quadl,qolAdDepCd42l,qolAdDepCd43l),file='coefTables/qolVscd4AndadhereAndDepressionLmer.html')
