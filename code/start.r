library(ggplot2)
library(foreign)
library(dplyr)
library(readr)


#cd4 <- read.csv('cd4.csv')
#dat <- read.spss('CINICS altogether+HbA1C.sav',to.data.frame=TRUE)
dat <- read.spss('CINICS altogether no HbA1C.sav',to.data.fram=TRUE)

dat$anxlevel <- factor(as.numeric(gsub(' ','',as.character(dat$anxlevel))),ordered=TRUE)
dat$Dep_total_G <- ordered(dat$Dep_total_G)
dat$cig3c <- factor(as.numeric(gsub(' ','',as.character(dat$cig3c))),ordered=TRUE)

spss2date <- function(x) as.POSIXct(x, origin = "1582-10-14")
dat$Date <- spss2date(dat$Date)


### note (maybe for later): take all dates seriously? w/ diagnosis unclear

#### plot depression levels
dep <- dat%>%filter(!is.na(dep_total),!is.na(Date),!is.na(StudyId))%>%group_by(StudyId)%>%
    mutate(days=Date-min(Date),dep_total,nmeas=n(),timeSpan=max(Date)-min(Date),depLag=lag(Dep_total_G))



fit <- lowess(dep$days,dep$dep_total)

fit <- as.data.frame(fit)%>%group_by(x)%>%summarize(y=mean(y))


dep$resid <- dep$dep_total-fit$y[match(dep$days,fit$x)]

ggplot(filter(dep,nmeas>2),aes(days,resid,group=StudyId))+geom_jitter()

## peopple w maximum residual at day 0
dep%>%
  ungroup()%>%
  filter(nmeas>2)%>%
  mutate(mr=max(resid[days==0]))%>%
  group_by(StudyId)%>%
  filter(resid[days==0]==mr)%>%
  ungroup()%>%
  ggplot(aes(days,resid,group=factor(StudyId),color=factor(StudyId)))+geom_line()

## choosing lines
sid <- dep%>%filter(nmeas>2)%>%arrange(resid)
sid <- sid$StudyId[c(1,floor(nrow(sid)*seq(0,1,length=10))[-1])]

plot(jitter(as.numeric(dep$days)),jitter(dep$dep_total),cex=0.3,col=alpha('grey',.5))

for(id in sid){
    ddd <- subset(dep,StudyId==id)
    lines(ddd$days,ddd$dep_total,type='b')
}


ggplot()+geom_jitter(data=dep,mapping=aes(days,dep_total,size=sid,alpha=sid),size=0.3,alpha=0.5)+
    geom_line(data=subset(dep,StudyId%in%sid),mapping=aes(days,dep_total,group=StudyId),size=1)+
        geom_point(data=subset(dep,StudyId%in%sid),mapping=aes(days,dep_total),inherit.aes=FALSE,size=2)

### depression & alcohol
depAlc <- dat%>%filter(!is.na(dep_total)|!is.na(alc_Sum),!is.na(Date),!is.na(StudyId))%>%
    group_by(StudyId)%>%
        mutate(days=Date-min(Date),ndep=sum(!is.na(dep_total)),nalc=sum(!is.na(alc_Sum)))%>%
            filter(ndep>1,nalc>1)%>%arrange(days)%>%mutate(depLag=lag(dep_total),alcLag=lag(alc_Sum))

### for each measurement of alc, find most recent prior measurement of dep
prevDep <- function(x){
    if((nrow(x)==1)){
        x$prevDate <- NA
        x$prevDep <- NA
        return(x)
    }
    if(all(is.na(x$dep_total))){
        x$prevDate <- NA
        x$prevDep <- NA
        return(x)
    }
    if(min(x$Date[!is.na(x$dep_total)])>=max(x$Date[!is.na(x$alc_Sum)])){
        x$prevDate <- NA
        x$prevDep <- NA
        return(x)
    }

    dep <- filter(x,!is.na(dep_total))
    alc <- filter(x,!is.na(alc_Sum))

    x$prevDate[!is.na(x$alc_Sum)] <- vapply(alc$Date, function(d) max(dep$Date[dep$Date<d]),dat$Date[1])
    pd <- dep$dep_total[match( na.omit(x$prevDate),dep$Date)]
    x$prevDep[!is.na(x$alc_Sum)]  <- pd



    x
}

depAlc <- dat%>%filter(!is.na(dep_total)|!is.na(alc_Sum),!is.na(Date),!is.na(StudyId))%>%
    select(StudyId,Date,alc_Sum,dep_total)%>%arrange(Date)%>%group_by(StudyId)%>%prevDep()

depAlc$Date <- spss2date(depAlc$Date)
depAlc$prevDate <- spss2date(depAlc$prevDate)


#depAlc <- subset(depAlc,is.finite(alc_Sum)&is.finite(prevDep))

ggplot(depAlc,aes(prevDep,alc_Sum))+geom_jitter(width=0.25,height=0.25)+geom_smooth()

library(lme4)

mod <- lmer(alc_Sum~prevDep+(1|StudyId),data=depAlc)
depAlc$diff <- depAlc$Date-depAlc$prevDate
depAlc <- depAlc%>%group_by(StudyId)%>%mutate(alcLag=lag(alc_Sum))

depAlc$diff <- depAlc$diff-mean(depAlc$diff[is.finite(depAlc$diff)])

mod2 <- lmer(alc_Sum~prevDep*diff+(1|StudyId),data=depAlc)
mod3 <- update(mod2,.~.+alcLag)


### add covariates
mna <- function(x) mean(x,na.rm=TRUE)
nonmis1 <- function(x) na.omit(x)[1]
covs <- dat%>%group_by(StudyId)%>%summarize(BirthYear=mna(BirthYear),birthSex=na.omit(birthSex_R)[1],
  trans=any(birthSex_R!=PresentSex_R,na.rm=TRUE),latinx=nonmis1(Hispanic_R),race=nonmis1(Race_R))

Mode <- function(x) levels(na.omit(as.factor(x)))[which.max(table(as.factor(x)))]
covs <- as.data.frame(covs)
for(cc in 2:ncol(covs)) covs[is.na(covs[,cc]),cc] <- Mode(covs[,cc])

depAlc <- left_join(depAlc,covs)
depAlc$age <- lubridate::year(depAlc$Date)-as.numeric(depAlc$BirthYear)
mod4 <- update(mod3,.~.+age+birthSex+trans+latinx+race)

depAlc <- depAlc%>%mutate(nalc=sum(!is.na(alc_Sum)))

mod5 <- update(mod4,.~.+nalc)


### diabetes?
diagnosis <- read.csv('Diagnosis.csv')
diagId <- unique(diagnosis$StudyId[diagnosis$DxCategory%in%c('Diabetes: Type 2','Diabetes: Type unspecified')])

depAlc$diab <- depAlc$StudyId%in%diagId

summary(mod6 <- update(mod5,.~.+prevDep*diab))


### depression & alcohol
depAlc <- dat%>%filter(!is.na(dep_total)|!is.na(alc_Sum),!is.na(Date),!is.na(StudyId))%>%
    group_by(StudyId)%>%
        mutate(days=Date-min(Date),ndep=sum(!is.na(dep_total)),nalc=sum(!is.na(alc_Sum)))%>%
            filter(ndep>1,nalc>1)%>%arrange(days)%>%mutate(depLag=lag(dep_total),alcLag=lag(alc_Sum))

### for each measurement of alc, find most recent prior measurement of dep
prevDep <- function(x){
     x$prevDate <- x$Date
     x$prevDep <- NA
     x$prevAlc0 <- NA
     x$prevAlc1 <- NA
    if((nrow(x)==1)){
        return(x)
    }
    if(all(is.na(x$dep_total))){
        return(x)
    }
    if(min(x$Date[!is.na(x$dep_total)])>=max(x$Date[!is.na(x$alc_Sum)])){
        return(x)
    }


     for(i in 2:nrow(x)){  ### works cuz arranged by date
         j <- max(which(!is.na(x$dep_total[1:(i-1)])))
         if(!is.finite(j)) next
         x$prevDate[i] <- x$Date[j]
         x$prevDep[i] <- x$dep_total[j]
         x$prevAlc0[i] <- x$alc_Sum[j]
         if(j>1) x$prevAlc1[i] <- x$alc_Sum[max(which(!is.na(x$alc_Sum[1:(j-1)])))]
     }

     ## aaa <- sapply(x$Date, function(dd){
     ##                   j <- which.max(x$Date[x$Date<dd & !is.na(x$dep_total)])
     ##                   if(length(j)) if(is.finite(j)) return(j)
     ##                   NA
     ##               })
     ## x$PrevDate <- x$Date[aaa]
     ## x$prevDep <- x$dep_total[aaa]
     x

 }


depAlc <- dat%>%filter(!is.na(dep_total)|!is.na(alc_Sum),!is.na(Date),!is.na(StudyId))%>%
    select(StudyId,Date,alc_Sum,dep_total)%>%arrange(Date)%>%group_by(StudyId,Date)%>%
        summarize(alc_Sum=ifelse(any(!is.na(alc_Sum)),mean(alc_Sum,na.rm=TRUE),NA),
                  dep_total=ifelse(any(!is.na(dep_total)),mean(dep_total,na.rm=TRUE),NA))%>%ungroup()

### checking
as.data.frame(depAlc%>%filter(StudyId==1007264654)%>%prevDep())

depAlc <- do.call('rbind',lapply(split(depAlc,depAlc$StudyId),prevDep))
depAlc$prevDate[depAlc$prevDate==depAlc$Date] <- NA
depAlc$diff <- depAlc$Date-depAlc$prevDate
depAlc$prevAlc <- ifelse(is.na(depAlc$prevAlc0),depAlc$prevAlc1,depAlc$prevAlc0)

ggplot(depAlc,aes(prevDep,alc_Sum))+geom_jitter(width=0.25,height=0.25)+geom_smooth()

depAlc$diff <- as.numeric(depAlc$diff)
with(depAlc,plot(table(diff[diff<200])))

depAlc$diffC <- cut(depAlc$diff,c(0,91,182,294,371,735,Inf))
depAlc <- depAlc%>%group_by(StudyId)%>%mutate(nobs=sum(!is.na(alc_Sum)&!is.na(prevDep)))

with(depAlc,plot(jitter(nobs),jitter(dep_total)))
with(depAlc,plot(jitter(nobs),jitter(alc_Sum)))

depAlc$nobsC <- depAlc$nobs
depAlc$nobsC[depAlc$nobsC>3] <- 3

ggplot(depAlc%>%filter(!is.na(diffC),nobsC>0),aes(dep_total-prevDep,alc_Sum-prevAlc))+geom_jitter(width=0.25,height=0.25)+geom_smooth(method='lm')+facet_grid(nobsC~diffC)

library(lme4)
depAlc$diff2 <- depAlc$diff/100
mod1 <- lmer(alc_Sum~prevDep*diff2+prevAlc+(1|StudyId),data=depAlc)

ggplot(dat,aes(dep_total,alc_Sum))+geom_jitter(width=0.25,height=0.25)+geom_smooth(method='lm')



mod <- lmer(alc_Sum~prevDep+(1|StudyId),data=depAlc)
depAlc$diff <- depAlc$Date-depAlc$prevDate
depAlc <- depAlc%>%group_by(StudyId)%>%mutate(alcLag=lag(alc_Sum))

depAlc$diff <- depAlc$diff-mean(depAlc$diff[is.finite(depAlc$diff)])

mod2 <- lmer(alc_Sum~prevDep*diff+(1|StudyId),data=depAlc)
mod3 <- update(mod2,.~.+alcLag)


### add covariates
mna <- function(x) mean(x,na.rm=TRUE)
nonmis1 <- function(x) na.omit(x)[1]
covs <- dat%>%group_by(StudyId)%>%summarize(BirthYear=mna(BirthYear),birthSex=na.omit(birthSex_R)[1],
  trans=any(birthSex_R!=PresentSex_R,na.rm=TRUE),latinx=nonmis1(Hispanic_R),race=nonmis1(Race_R))

Mode <- function(x) levels(na.omit(as.factor(x)))[which.max(table(as.factor(x)))]
covs <- as.data.frame(covs)
for(cc in 2:ncol(covs)) covs[is.na(covs[,cc]),cc] <- Mode(covs[,cc])

depAlc <- left_join(depAlc,covs)
depAlc$age <- lubridate::year(depAlc$Date)-as.numeric(depAlc$BirthYear)
mod4 <- update(mod3,.~.+age+birthSex+trans+latinx+race)

depAlc <- depAlc%>%mutate(nalc=sum(!is.na(alc_Sum)))

mod5 <- update(mod4,.~.+nalc)


### diabetes?
diagnosis <- read.csv('Diagnosis.csv')
diagId <- unique(diagnosis$StudyId[diagnosis$DxCategory%in%c('Diabetes: Type 2','Diabetes: Type unspecified')])

depAlc$diab <- depAlc$StudyId%in%diagId

summary(mod6 <- update(mod5,.~.+prevDep*diab))




##depression
depAnx <- read_csv('Pro_DepressionAnxiety.csv')
dep <- depAnx[,paste0('dep',1:9)]
dep[dep==''] <- NA

for(i in 1:ncol(dep))
    levels(dep[,i]) <- c('not at all','several days','more than half of days','nearly every day')

depNum <- matrix(nrow=nrow(dep),ncol=ncol(dep))
depNum[dep=='not at all'] <- 0
depNum[dep=='several days'] <- 1
depNum[dep=='more than half of days'] <- 2
depNum[dep=='nearly every day'] <- 3

dep <- rowMeans(depNum,na.rm=TRUE)*ncol(depNum)



orna <- function(x,y) mean(is.na(x)*is.na(y))

files <- list.files()

csvs <- files[grep('.csv',files,fixed=TRUE)]

dataFiles <- list()

for(ff in csvs){
    print(ff)
    dataFiles[[ff]] <- read_csv(ff)
}

for(ff in csvs){
    print(ff)
    print(head(dataFiles[[ff]]))
}



## abstract due 2/18

#### coarsen time variable (month), turn into regular time points w/ missing data


################ similar thing, but with adhere
prevDepAd <- function(x){
     x$prevDate <- x$Date
     x$prevDep <- NA
     x$prevAd0 <- NA
     x$prevAd1 <- NA
    if((nrow(x)==1)){
        return(x)
    }
    if(all(is.na(x$dep_total))){
        return(x)
    }
    if(min(x$Date[!is.na(x$dep_total)])>=max(x$Date[!is.na(x$adhere)])){
        return(x)
    }


     for(i in 2:nrow(x)){  ### works cuz arranged by date
         j <- max(which(!is.na(x$dep_total[1:(i-1)])))
         if(!is.finite(j)) next
         x$prevDate[i] <- x$Date[j]
         x$prevDep[i] <- x$dep_total[j]
         x$prevAd0[i] <- x$adhere[j]
         if(j>1) x$prevAd1[i] <- x$adhere[max(which(!is.na(x$adhere[1:(j-1)])))]
     }

     ## aaa <- sapply(x$Date, function(dd){
     ##                   j <- which.max(x$Date[x$Date<dd & !is.na(x$dep_total)])
     ##                   if(length(j)) if(is.finite(j)) return(j)
     ##                   NA
     ##               })
     ## x$PrevDate <- x$Date[aaa]
     ## x$prevDep <- x$dep_total[aaa]
     x

 }



depAdhere <- dat%>%filter(!is.na(dep_total)|!is.na(adhere),!is.na(Date),!is.na(StudyId))%>%
    select(StudyId,Date,adhere,dep_total)%>%arrange(Date)%>%group_by(StudyId)%>%prevDepAd()

depAdhere$Date <- spss2date(depAdhere$Date)
depAdhere$prevDate <- spss2date(depAdhere$prevDate)


#depAlc <- subset(depAlc,is.finite(alc_Sum)&is.finite(prevDep))

ggplot(depAdhere,aes(prevDep,adhere))+geom_jitter(width=0.25,height=0.25)+geom_smooth()

library(lme4)

modAd <- lmer(adhere~prevDep+(1|StudyId),data=depAdhere)
depAdhere$diff <- depAdhere$Date-depAdhere$prevDate
depAdhere <- depAdhere%>%group_by(StudyId)%>%mutate(adLag=lag(adhere))

depAdhere$diff <- depAdhere$diff-mean(depAdhere$diff[is.finite(depAdhere$diff)])

mod2 <- lmer(adhere~prevDep*diff+(1|StudyId),data=depAdhere)
mod3 <- update(mod2,.~.+adLag)

