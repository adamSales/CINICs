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

save(list=grep('cd4AdDep',ls()),file='fittedModels/cd4AdDepModels.RData')
htmlreg(list(cd4AdDep00,cd4AdDep0,cd4AdDep,cd4AdDep2,cd4AdDep3),file='coefTables/cd4VsadhereAndDepression.html')

##### Plots
dat%>%
  mutate(rldep=rank(logDepLag))%>%
  select(cd4,'vs. Adherence (lagged)'=adhere_lag,'vs. rank(Depression+1) (lagged)'=rldep)%>%
  gather("iv","value",-cd4)%>%
  ggplot(aes(factor(value,exclude=NA),cd4))+
  geom_boxplot()+
  geom_smooth(aes(value,cd4),method='lm')+
  facet_wrap(~iv,scale="free_x")+
  labs(x=NULL,y='CD4 Count')+
  scale_x_discrete(na.translate=FALSE)
ggsave('plots/cd4VsadhereAndDepression.jpg')


nd1 <- within(dat,{logDepLag=1.06;depDist=0;adhereDist=0;cd4_lag=531;cd4Dist=0})[rownames(model.frame(cd4AdDep)),]
ppp1 <- predict(cd4AdDep,nd)
nd2 <- within(dat,{adhere_lag=mean(model.frame(cd4AdDep)$adhere_lag);depDist=0;adhereDist=0;cd4_lag=531;cd4Dist=0})[rownames(model.frame(cd4AdDep)),]
ppp2 <- predict(cd4AdDep,nd)

bind_rows(
  tibble(
  mutate(cd4_100=ppp+resid(cd4AdDep))%>%
  ggplot(aes(factor(adhere_lag,exclude=NA),cd4_100))+
  geom_boxplot()+
  geom_smooth(aes(adhere_lag,cd4_100),method='lm')+
  labs(x='Lagged Adherence',y='CD4 Count/100')+
  scale_x_discrete(na.translate=FALSE)+ylim(1,5)

p2 <- dat%>%
  ggplot(

grid.arrange(p1,p2,nrow=1)
