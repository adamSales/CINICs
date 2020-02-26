library(zoo)
library(tidyverse)
library(haven)
library(lubridate)

#datO <- read_sav('data/CINICS altogether no HbA1C.sav')
#dat <- datO
dat <- read_sav('data/CINICS altogether no HbA1C.sav')
for(i in which(sapply(dat,is.character))) dat[[i]][dat[[i]]==''] <- NA

naclass <- function(cls){
  x <- NA
  class(x) <- cls
  x
}
naForCls <- sapply(unique(sapply(dat,class)),naclass,simplify=FALSE)

dat <-
  dat%>%filter(!is.na(Date))%>%
  mutate_at(c('cig3c','anxlevel'),as.numeric)

recodeHavenLabelled <- function(x){
  if(!inherits(x,"haven_labelled" )) return(x)
  names(attr(x,'labels'))[match(x,attr(x,'labels'))]
}


### what subject characteristics are stable?
## stable <- dat%>%group_by(StudyId)%>%mutate(nobs=n())%>%filter(nobs>1)%>%
##  summarize_at(vars(BirthYear,ends_with('_R'),Site),n_distinct,na.rm=TRUE)

## lapply(stable[,-1],table)

## create stable covariate dataset, one line per subject
covs <-
  dat%>%
  group_by(StudyId)%>%
  mutate_if(~inherits(.,'haven_labelled'),recodeHavenLabelled)%>%
  mutate(BirthYear=as.numeric(BirthYear))%>%
  summarize_at(
    vars(BirthYear:Race_R,Site),
    ~ if(all(is.na(.))) .[1] else unique(na.omit(.))[1]
  )%>%
  mutate(
    Site=ifelse(is.na(Site),"NA",Site),
    gender=ifelse(is.na(PresentSex_R),birthSex_R,PresentSex_R)
    )

## for variables that change over time, we want one row per subject/date:
dat <- dat%>%
  select(!!c('StudyId',setdiff(names(dat),names(covs)),'BirthYear'))%>%
  group_by(StudyId,Date)%>%
  summarize_all(~ ifelse(all(is.na(.)),.[1],unique(na.omit(.))[1]))%>%
  mutate(age=lubridate::year(Date)-BirthYear)%>%
  ungroup()%>%
  select(-BirthYear)

## ndate <- dat%>%group_by(StudyId)%>%summarize(ndate=n_distinct(Date))

## dateSep <- dat%>%group_by(StudyId)%>%arrange(Date)%>%mutate(ndate=n(),dateGap=Date-lag(Date))%>%
##   summarize(ndate=ndate[1],
##     minGap=ifelse(ndate==1,-1,min(dateGap,na.rm=TRUE)),
##     maxGap=ifelse(ndate==1,-1,max(dateGap,na.rm=TRUE)))

## plotSep <- dateSep%>%filter(ndate>1)%>%gather(minMax,gap,minGap,maxGap)

## ggplot(plotSep%>%filter(gap<50),aes(gap))+geom_histogram(bins=100)+facet_wrap(~minMax)

## quantile(dateSep$minGap[dateSep$minGap>0],seq(0,1,.05))


### numeric variables that are almost always (>90%) equal to 0 -> dichotomous
dat <- dat%>%
    mutate_if(
        function(x) is.numeric(x)&mean(x==0,na.rm=TRUE)>0.9,
        function(x) ifelse(x>0,1,0)
    )

### combine hard drugs
dat <- mutate(dat,anyHard=ifelse(OPI_SUM>0|COC_SUM>0|AMPH_SUM>0,1,0))

alc <- read_csv('data/PRO_AlcoholAudit.csv')

alc <- alc%>%
  mutate(
    alcCat=ifelse(alcriskhi=='',NA,
      ifelse(alcrisklo=='not at risk',1,
      ifelse(alcriskhi=='not at risk',2,3))),
    StudyId=as.numeric(StudyId)
  )%>%
  select(StudyId,Date,alcCat)%>%
  filter(StudyId%in%dat$StudyId)

dat <- full_join(dat,alc)


adhere <- read_csv('data/PRO_Adherence.csv')%>%
  mutate(
    adhere=ifelse(selfrt=='',NA,
      ifelse(selfrt=='Very poor'|selfrt=='Poor',1,
      ifelse(selfrt=='Fair',2,
      ifelse(selfrt=='Good',3,
      ifelse(selfrt=='Very good',4,5))))),
    StudyId=as.numeric(StudyId))%>%
  select(StudyId,Date,adhere)%>%
  filter(StudyId%in%dat$StudyId)

dat <- full_join(dat,adhere)

diagnosis <- read_csv('data/Diagnosis.csv')%>%
  group_by(StudyId)%>%
  summarize(diab=ifelse(any(DxCategory=='Diabetes: Type 2'|DxCategory=='Diabetes: Type unspecified'),1,0))%>%
  mutate(StudyId=as.numeric(StudyId))

dat <- left_join(dat,diagnosis)

cd4 <- read_csv('data/Lab.csv')%>%
  filter(TestName=='CD4 cell absolute')%>%
  transmute(
    Date=date(ResultDate),
    cd4=as.numeric(Result),
    StudyId=as.numeric(StudyId))%>%
  filter(!is.na(cd4),cd4<100000,StudyId%in%dat$StudyId)
## one cd4 count is almost 1M which is not possible and skews everything

### some duplicate dates for some reason:
# ns <- cd4%>%group_by(StudyId)%>%summarize(nrec=n(),ndate=n_distinct(Date))
# sum(ns$nrec-ns$ndate)
# [1] 369
cd4 <- cd4%>%group_by(StudyId,Date)%>%summarize(cd4=round(mean(cd4,na.rm=TRUE)))

dat <- full_join(dat,cd4)

ns <- dat%>%group_by(StudyId)%>%summarize(nrec=n(),ndate=n_distinct(Date))
mean(ns$nrec==ns$ndate)

## merge covs back in
dat <- full_join(dat,covs,by='StudyId')%>%
    mutate(age=lubridate::year(Date)-BirthYear)


save(dat,file='data/cleanedData.RData')
