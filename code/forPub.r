library(lme4)
library(lmerTest)
library(r2glmm)
library(tidyverse)
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
  ungroup()%>%
  filter(is.na(cd4_lag)|cd4_lag<3000,is.na(cd4)|cd4<3000)%>%  ## these are extreme outliers
  droplevels()



## dat2 <- dat%>%mutate(logDepLag=logDepLag-mean(logDepLag,na.rm=TRUE),
##   adhere_lag=adhere_lag-mean(adhere_lag,na.rm=TRUE),
##   cd4_100_lag=cd4_100_lag-mean(cd4_100_lag,na.rm=TRUE))

## load('fittedModels/AdDepModels.RData')
## load('fittedModels/cd4AdDepModels.RData')
## load('fittedModels/qolAdDepCd4LinearModels.RData')

### refit models for report
## first standardize main study variables: logDepLag, adhere,cd4,qol

#dat <- dat[-c(which.max(datS$cd4),which.min(datS$cd4)),]

datS <- dat%>%
  mutate(logDepLag=logDepLag/sd(logDepLag,na.rm=TRUE))
for(vv in c('adhere','cd4','qol1_5_total')){
  SD <- sd(dat[[vv]],na.rm=TRUE)
  datS[[vv]] <- datS[[vv]]/SD
  datS[[paste0(vv,'_lag')]] <-   datS[[paste0(vv,'_lag')]]/SD
}

mods <- list(
  Adhere = lmer(adhere~logDepLag+depDist+adDist+adhere_lag+(adDist+depDist|StudyId),data=dat),
  CD4 = lmer(cd4~logDepLag+depDist+adhere_lag+adDist+cd4_lag+cd4Dist+(depDist+adDist+cd4Dist|StudyId),data=dat),
  QOL = lmer(qol1_5_total~logDepLag+depDist+adhere_lag+adDist+cd4_lag+cd4Dist+qol1_5_total_lag+qolDist+(cd4Dist+adDist|StudyId),data=dat)
    )

modsS <- map(mods,update,data=datS)

### adjust mods$cd4 for heteroskedasticity in adherence?
mf <- model.frame(modsS$CD4)
adVal <- sort(unique(mf$adhere_lag))
for(i in 2:5){
    val <- adVal[i]
    mf[[paste0('ad',i)]] <- ifelse(mf$adhere_lag==val,seq(nrow(mf)),0)
}

cd4Mod <- update(modsS$CD4,.~.+(1|ad2)+(1|ad3)+(1|ad4)+(1|ad5),data=mf)

save(cd4Mod,file='fittedModels/cd4hetero.RData')


save(mods,modsS,file='fittedModels/modsForPub.RData')

### table 1
pn <- function(x) formatC(x,digits=1,format='f')
cap <- function(x) paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
tab1dat <-
  dat%>%
    filter(StudyId%in%model.frame(mods$Adhere)$StudyId)%>%
      group_by(StudyId)%>%
        summarize(Age=mean(age,na.rm=TRUE),
                  Gender=na.omit(PresentSex_R)[1],
                  Race=ifelse(all(is.na(Hispanic_R))|na.omit(Hispanic_R)[1]=='no',na.omit(Race_R)[1],'Latino'),
                  Adherence=mean(adhere,na.rm=TRUE),
                  CD4=mean(cd4,na.rm=TRUE),
                  Depression=mean(dep_total,na.rm=TRUE),
                  `Log(Depression)`=mean(log(dep_total+1),na.rm=TRUE),
                  QOL=mean(qol1_5_total,na.rm=TRUE)
                  )%>%
          select(-StudyId)

tab1 <- map_dfr(
  names(tab1dat),
  function(nn) {
    x <- tab1dat[[nn]]
    if(is.numeric(x))
      return(
        tibble(
          Variable=nn,
          Report='Mean (SD)',
          Value=paste0(pn(mean(x,na.rm=TRUE)),' (',pn(sd(x,na.rm=TRUE)),')')
          )
        )
    map_dfr(
      unique(na.omit(x)),
      ~tibble(
        Variable=ifelse(.==unique(na.omit(x))[1],nn,''),
        Report='N (%)',
        Value=paste0(cap(.),': ',prettyNum(sum(x==.,na.rm=TRUE),','),' (',
                            pn(mean(x==.,na.rm=TRUE)*100),'%)'))
      )
  }
  )

write.csv(tab1,'table1.csv',row.names=FALSE)

coefNames <- c(
    logDepLag="log(Dep.) (Lagged)",
    depDist="Dep. Meas. Gap",
    adDist="Adhere Meas. Gap",
    adhere_lag="Adherence (Lagged)",
    cd4_lag="CD4 (Lagged)",
    cd4Dist="CD4 Meas. Gap",
    qol1_5_total_lag="QOL (Lagged)",
    qolDist="QOL Meas. Gap"
    )

makeTrObj <- function(mod){
  tr <- extract(mod)
  tr@gof <- tr@gof[tr@gof.names%in%c("AIC","BIC","Log Likelihood","Num. obs.","Num. groups: StudyId")]
  tr@gof.decimal <- tr@gof.decimal[tr@gof.names%in%c("AIC","BIC","Log Likelihood","Num. obs.","Num. groups: StudyId")]
  tr@gof.names <- tr@gof.names[tr@gof.names%in%c("AIC","BIC","Log Likelihood","Num. obs.","Num. groups: StudyId")]
  tr@gof <- c(tr@gof,r2beta(mod,partial=FALSE)[1,'Rsq'])
  tr@gof.decimal <- c(tr@gof.decimal,TRUE)
  tr@gof.names <- c(tr@gof.names,'R2')

  tr@coef.names <- ifelse(tr@coef.names%in% names(coefNames),
                          coefNames[tr@coef.names],
                          tr@coef.names
                          )
  tr
}

modCoefs <- map(modsS,makeTrObj)

screenreg(modCoefs)

modTable <- htmlreg(modCoefs,file='coefTables/forPub.doc',doctype=TRUE,html.tag=TRUE,head.tag=TRUE,body.tag=TRUE,
                    reorder.coef=c(1,2,5,6,8,3,4,7,9))


### plots
predRange <- list()

pr <- function(mod,nd,iv){
    iv <- sym(iv)
     nd%>%arrange(!! iv)%>%select(names(fixef(mod))[-1])%>%slice(1,n())%>%as.matrix()%>%cbind(c(1,1),.)%*%fixef(mod)
}
### adherence vs depression
nd <- within(model.frame(mods$Adhere),{depDist=0;adhere_lag=mean(dat$adhere_lag,na.rm=TRUE);adDist=0})
ppp <- predict(mods$Adhere,nd)

datPoint <- model.frame(mods$Adhere)%>%
    mutate(adhere=ppp+resid(mods$Adhere),
           dep_total_lag=exp(logDepLag))

datAv <- datPoint%>%
    group_by(dep_total_lag)%>%summarize(n=n(),avg.adhere=mean(adhere,na.rm=TRUE))

predRange$ad <- pr(mods$Adhere,nd,"logDepLag")
## ggplot(datPoint,aes(dep_total_lag+1,adhere))+
##     geom_jitter(alpha=0.1)+
ggplot(datAv,aes(dep_total_lag+1,avg.adhere,size=n))+
    geom_point(show.legend=FALSE)+
            geom_smooth(method='lm',show.legend=FALSE,color='grey',se=FALSE,fullrange=TRUE)+
                scale_size(range=c(.5,5))+
                    ylim(1,5)+
                        scale_x_continuous(trans='log',breaks=c(1:8,10,15,20,28),minor_breaks=NULL)+
                            xlab('Depression+1 (lagged)')+ylab('Avg. Adherence (Adjusted)')
ggsave('plots/forPub/adhereVsDepressionAdj.png',width=3,height=3)

### cd4 vs adherence and depression
nd1 <- within(model.frame(mods$CD4),{logDepLag=mean(logDepLag,na.rm=TRUE);depDist=0;adDist=0;cd4_lag=mean(cd4_lag,na.rm=TRUE);cd4Dist=0})
ppp1 <- predict(mods$CD4,nd1)
predRange$cd4$ad <- pr(mods$CD4,nd1,"adhere_lag")
nd2 <- within(model.frame(mods$CD4),{adhere_lag=mean(dat$adhere_lag,na.rm=TRUE);depDist=0;adDist=0;cd4_lag=mean(cd4_lag,na.rm=TRUE);cd4Dist=0})
ppp2 <- predict(mods$CD4,nd2)
predRange$cd4$dep <- pr(mods$CD4,nd2,"logDepLag")



pdat <-
    model.frame(mods$CD4)%>%
  select(cd4,'vs. Adherence (lagged)'=adhere_lag,'vs. Depression (lagged)'=logDepLag)%>%
  gather("iv","value",-cd4)%>%
  mutate(
      width=.4*(log(exp(value)+1)-value),
      jitter=value+map_dbl(1:n(), ~if(iv[.]=='vs. Adherence (lagged)') runif(1,-.4,.4) else runif(1,-width[.],width[.])),
      bpVal=value,#as.factor(value),
      cd4=(ifelse(iv=='vs. Adherence (lagged)',ppp1+resid(mods$CD4),ppp2)+resid(mods$CD4)))#%>%

p <- ggplot(pdat,aes(jitter,cd4))+
  geom_point(alpha=0.2)+
  geom_smooth(aes(value,cd4),method='lm',color='grey',se=FALSE,fullrange=TRUE)+
  facet_wrap(~iv,scale="free_x")+
  labs(x=NULL,y='CD4 Count (Adjusted)')+
  scale_x_continuous(
      breaks=function(x) if(max(x)>4) 1:5 else log(c(1:8,10,15,20,28)),
      labels=function(x) if(max(x)>4) 1:5 else c(1:8,10,15,20,28)
      )+
ylim(0,max(c(ppp1,ppp2)+resid(mods$CD4)))

ggsave('plots/forPub/cd4VsadhereAndDepressionAdj.png',plot=p,width=6,height=3)


### QOL vs cd4, adherence and depression
preds <-     c('cd4_lag','adhere_lag','logDepLag')%>%setNames(.,.)
ndQOL <- map(preds,
             function(x)
                 model.frame(mods$QOL)%>%
                     mutate_at(vars(ends_with('lag'),-!! sym(x)),mean,na.rm=TRUE)%>%
                         mutate_at(vars(ends_with('Dist')),~0)
             )
predRange$qol <- map(preds,~pr(mods$QOL,ndQOL[[.]],.))
ppp <- map_dfc(preds,~predict(mods$QOL,newdata=nd$QOL[[.]]))

   ## pivot_longer(cols=everything(),names_to="iv",values_to="qol")%>%
   ## arrange(iv)



p <- model.frame(mods$QOL)%>%
    ## mutate(
    ##     cd4_lag=as.numeric(as.character(cut(cd4_lag,seq(0,3000,50),labels=seq(0,2950,50)+25))),
    ## )%>%
    select(qol=qol1_5_total,cd4_lag,adhere_lag,logDepLag)%>%
    gather("iv","value",-qol)%>%
    mutate(
      #width=.4*(log(exp(value)+1)-value),
      #jitter=value+map_dbl(1:n(), ~if(iv[.]=='adhere_lag') runif(1,-.4,.4) else if(iv[.]=='logDepLag') runif(1,-width[.],width[.]) else 0),
      bpVal=value,#as.factor(value),
      qol=ifelse(iv=='cd4_lag',ppp[['cd4_lag']]+resid(mods$QOL),
          ifelse(iv=='adhere_lag',ppp[['adhere_lag']]+resid(mods$QOL),
                 ppp[['logDepLag']]+resid(mods$QOL))
                 )-10
    )%>%
    group_by(iv,value)%>%
    summarize(n=n(),avg.qol=mean(qol,na.rm=TRUE))%>%
    ungroup()%>%
        mutate(iv=c(cd4_lag='vs. CD4 Count (lagged)',adhere_lag='vs. Adherence (lagged)',logDepLag='vs. Depression (lagged)')[iv])%>%
  ggplot(aes(value,avg.qol))+
  geom_point(aes(size=n),show.legend=FALSE)+
  geom_smooth(aes(value,avg.qol,weight=n),method=lm,color='grey',se=FALSE,fullrange=TRUE)+#,formula=y~poly(x,2))+
  #geom_segment(mapping=aes(x=xstart,y=ystart,xend=xend,yend=yend),color='grey',inherit.aes=FALSE)+
      facet_wrap(~iv,scales="free_x",ncol=2)+
  labs(x=NULL,y='Avg. QOL (adjusted)')+
        scale_x_continuous(
      breaks=function(x) if(max(x,na.rm=TRUE)>100)  seq(0,2500,500) else if(max(x,na.rm=TRUE)>4) 1:5 else log(c(1:8,10,15,20,28)),
      labels=function(x) if(max(x,na.rm=TRUE)>100)  seq(0,2500,500) else if(max(x,na.rm=TRUE)>4) 1:5 else c(1:8,10,15,20,28),
      minor_breaks=function(x) if(max(x,na.rm=TRUE)<5) log(1:28) else NULL)+
  scale_size(range=c(.5,5))+
  ylim(0,5)

ggsave('plots/forPub/qolVscd4adhereAndDepressionAdj.png',width=6,height=6)


save(predRange,file='fittedModels/predRanges.RData')




modsS$CD4 <- cd4Mod


