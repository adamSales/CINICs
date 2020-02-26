library(nlme)
library(ipw)
library(tidyverse)


#Multilevel panel regression with cross-lagged effects will examine the associations among psychological and behavioral factors and health outcomes over time for HIV positive patients with and without diabetes.

load('data/cleanedData.RData')

table(duplicated(dat[,c('StudyId','time')]))

## do people switch sites?
dat%>%group_by(StudyId)%>%summarize(nsite=n_distinct(Site,na.rm=TRUE))%>%xtabs(~nsite,.)

site <- dat%>%group_by(StudyId)%>%summarize(Site=na.omit(Site)[1])
## function to combine records from the same day
comb <- function(x){
  if(all(is.na(x))) return(NA)
  x <- na.omit(x)
  if(length(x)==1) return(x)
  if(n_distinct(x)==1) return(x[1])
  if(!is.numeric(x)) return(x[1])
  mean(x)
}

dat <- dat%>%filter(!is.na(Date))%>%select(-Date,-BirthYear,-Site)%>%group_by(StudyId,time)%>%summarize_all(comb)

#### effect of dep on adhere controlling for previous cd4, substance abuse, dep
dat$cd4Sqrt <- sqrt(dat$cd4)

