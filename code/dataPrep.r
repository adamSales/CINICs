#source('cleanData.r')
library(zoo)
library(tidyverse)
library(haven)
library(lubridate)

load('data/cleanedData.RData')

naDate <- as.Date(NA)

dat <- dat%>%
  mutate(rownum=1:n())%>%
  filter(Date>as.Date("1900-01-01"))%>%
  arrange(Date)%>%
  group_by(StudyId)%>%
  mutate(lagDate=lag(Date))%>%
  mutate_at(
    vars(-Date,-StudyId,-lagDate,-BirthYear,-ends_with('_R'),-Site,-age,-rownum),
    list(
      lastObs=~na.locf(if_else(!is.na(lag(.)),lagDate,naDate),na.rm=FALSE),
      lag=~na.locf(lag(.),na.rm=FALSE)
    )
  )%>%
  select(-lagDate)%>%
  arrange(rownum)%>%
  select(-rownum)




save(dat,file='data/cleanedDataWithLags.RData')

## load('data/cleanedData.RData')

## ### coarsen date to 1/2 year within subject (~9% of between-visit intervals are less than this)
## dat <- dat%>%
##   group_by(StudyId)%>%
##   mutate(time=floor(as.numeric((Date-min(Date))/182.5)))%>%
##   ungroup()%>%group_by(StudyId,time)%>%
##   summarize(
##    BirthYear=BirthYear[1],
##    PresentSex_R=PresentSex_R[1],
##    birthSex_R=birthSex_R[1],
##    Hispanic_R=Hispanic_R[1],
##    Transgendered_R=Transgendered_R[1],
##    Race_R=Race_R[1],
##    Age_2018_birthyear=Age_2018_birthyear[1],
##    Race_T=Race_T[1],
##    Site=Site[1],
##    InsuranceType_R.1=InsuranceType_R.1[1],
##    visrt=mean(visrt,na.rm=TRUE),
##    hrqol6=mean(hrqol6,na.rm=TRUE),
##    qol1_5_total=mean(qol1_5_total,na.rm=TRUE),
##    pact4_R=mean(pact4_R,na.rm=TRUE),
##    age=mean(age,na.rm=TRUE),
##    ndate=ndate[1],
##    dep_total=max(dep_total,na.rm=TRUE),
##    Dep_total_G=max(Dep_total_G,na.rm=TRUE),
##    cig3c=max(cig3c,na.rm=TRUE),
##    alc_Sum=max(alc_Sum,na.rm=TRUE),
##    AMPH_SUM=max(AMPH_SUM,na.rm=TRUE),
##    COC_SUM=max(COC_SUM,na.rm=TRUE),
##    OPI_SUM=max(OPI_SUM,na.rm=TRUE),
##    POT_SUM=max(POT_SUM,na.rm=TRUE),
##    anxlevel=max(anxlevel,na.rm=TRUE))


## dat <- dat%>%group_by(StudyId)%>%mutate(ntimes=n())%>%filter(ntimes>1)%>%ungroup()


## ### cut out visits after 11th (.95 quantile)

## dat <- filter(dat,time<12)




## ### -Inf -> NA
## dat[dat==-Inf] <- NA

## dat$qol1_5_total <- round(dat$qol1_5_total)-10
## dat$pact4_R <- floor(dat$pact4_R) ## i.e. ever-yes in time period (2=no)

## ### distributions
## stabCat <- dat%>%group_by(StudyId)%>%select(gender,Hispanic_R:InsuranceType_R.1,-Transgendered_R)%>%
##   summarize_if(function(x) n_distinct(x,na.rm=TRUE)<10,function(x) x[1])%>%gather(demo,category,-StudyId)

## ggplot(stabCat,aes(demo))+geom_bar(aes(fill=category))+coord_flip()

## varyCat <- dat%>%select(visrt:anxlevel)%>%select_if(function(x) n_distinct(x,na.rm=TRUE)<10)%>%gather()%>%filter(is.finite(value))

## ggplot(varyCat,aes(key))+geom_bar(aes(fill=as.factor(value)))+coord_flip()

## varCont <- dat%>%select(visrt:anxlevel)%>%select_if(function(x) n_distinct(x,na.rm=TRUE)>9)%>%gather()%>%filter(is.finite(value))

## ggplot(varCont,aes(value))+geom_histogram()+facet_wrap(~key,scales='free')

## varCont0 <- dat%>%select(visrt:anxlevel)%>%select_if(function(x) n_distinct(x,na.rm=TRUE)>9)%>%
##   summarize_all(function(x) mean(x==0,na.rm=TRUE))




## varCont <- dat%>%select(visrt:anxlevel)%>%select_if(function(x) n_distinct(x,na.rm=TRUE)>9)%>%gather()%>%filter(is.finite(value))
## ggplot(varCont,aes(value))+geom_histogram()+facet_wrap(~key,scales='free')
## ggplot(varCont,aes(log(value+1)))+geom_histogram()+facet_wrap(~key,scales='free')


## varyCat <- dat%>%select(visrt:anxlevel,anyHard)%>%select_if(function(x) n_distinct(x,na.rm=TRUE)<10)%>%gather()%>%filter(is.finite(value))
## ggplot(varyCat,aes(key))+geom_bar(aes(fill=as.factor(value)))+coord_flip()

## save(dat,file='processedDat.RData')

## ## now, insert missingness for each possible time point for each subject
