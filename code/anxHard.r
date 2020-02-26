library(rstan)

options(mc.cores = min(parallel::detectCores(),10))
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')


## transform data to fit anxHard.stan
stanDat <- list()

dat$id <- as.numeric(as.factor(dat$StudyId))

stopifnot(max(dat$id)==n_distinct(dat$StudyId))
stanDat$N <- max(dat$id)

dat$time <- dat$time+1

stanDat$T <- 12

anxDat <- dat%>%select(id,time,anxlevel)%>%filter(is.finite(anxlevel))

stanDat$nobsAnx <- nrow(anxDat)
stanDat$anxObs <- anxDat$anxlevel
stanDat$subjAnx <- anxDat$id
stanDat$timeAnx <- anxDat$time

subDat <- dat%>%select(id,time,anyHard)%>%filter(is.finite(anyHard))

stanDat$nobsSub <- nrow(subDat)
stanDat$subObs <- subDat$anyHard
stanDat$subjSub <- subDat$id
stanDat$timeSub <- subDat$time

## fill in missing values by assumption
dat$race <- as.factor(dat$Race_R)
dat$race[is.na(dat$race)] <- 1
levels(dat$race) <- names(attr(datO$Race_R,'labels'))

dat$Hispanic_R[is.na(dat$Hispanic_R)] <- 2
dat$hisp <- ifelse(dat$Hispanic_R==1,1,0)

dat$gender[is.na(dat$gender)] <- 1

X <- model.matrix(~age+gender+race+hisp+Site,
  data=dat%>%group_by(id)%>%summarize_all(function(x) x[1]))

stanDat$X <- scale(X[,-1])

stanDat$p <- ncol(X)-1

mod <- stan('anxHard.stan',data=stanDat,chains=1,iter=100)
