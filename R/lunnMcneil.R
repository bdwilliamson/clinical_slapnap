# Lunn and McNeil (1995, Biometrics)
# this implementation assumes 2 failure types in addition to censoring
# the 2-sided test of {H0: VE against failure type 0 = VE against failure type 1} is evaluated using the
# the 2-sided interaction p-value for the term 'dbltreatment:dbldelta'

# INPUT:
# flrtime- right-censored failure time
# flrstatus- indicator of observed infection (1=yes, 0=no)
# flrtype- indicator of failure type (positive integers; cencode=0; no NAs)
# for censorings, flrtype can be either NA or the value of one of the two failure types (in the duplicated
# data-set, there will be two entries for each censoring, one for each failure type)
# Vx- Vaccination status (0=placebo, 1=vaccine)

lunnMcneilTest <- function(flrtime,flrstatus,flrtype,Vx) {
  flrtype <- ifelse(is.na(flrtype), flrtype, ifelse(flrtype>0, flrtype-1, NA))
  # flrtype = 0 or 1
  # Evaluate  H0: VE(type 0) = VE(type 1) via the Lunn and McNeil (1995, Biometrics) trick:
  currdelta <- flrstatus
  dblfutime <- c(flrtime,flrtime)
  dblfustat <- c(currdelta,rep(0,length(currdelta)))
  dbltreatment <- c(Vx,Vx)
  #firstdelta <- ifelse(currdelta==0 | (!is.na(flrtype) & flrtype==0),1,0)
  firstdelta <- ifelse(is.na(flrtype), 0, flrtype)
  seconddelta <- 1 - firstdelta
  dbldelta <- c(firstdelta,seconddelta)
    
  fit <- coxph(Surv(dblfutime,dblfustat) ~ dbltreatment*dbldelta)
  return(summary(fit))
}

# stratified Lunn and McNeil test
lunnMcneilTestS <- function(flrtime,flrstatus,flrtype,Vx,stratVar) {
  flrtype <- ifelse(is.na(flrtype), flrtype, ifelse(flrtype>0, flrtype-1, NA))
  # flrtype = 0 or 1
  # Evaluate  H0: VE(type 0) = VE(type 1) via the Lunn and McNeil (1995, Biometrics) trick:
  currdelta <- flrstatus
  dblfutime <- c(flrtime,flrtime)
  dblfustat <- c(currdelta,rep(0,length(currdelta)))
  dbltreatment <- c(Vx,Vx)
  dblStratVar <- rep(stratVar, 2)
  #firstdelta <- ifelse(currdelta==0 | (!is.na(flrtype) & flrtype==0),1,0)
  firstdelta <- ifelse(is.na(flrtype), 0, flrtype)
  seconddelta <- 1 - firstdelta
  dbldelta <- c(firstdelta,seconddelta)
  
  fit <- coxph(Surv(dblfutime,dblfustat) ~ dbltreatment*dbldelta + strata(dblStratVar))
  return(summary(fit))
}

lunnMcneilTest3gt <- function(ftime, fstatus, ftype, vaccine){
  # x is a scalar (one of 0,1,2); r is a scalar (one of 1,2)
  f <- function(x, r){
    s <- 0:2
    idx <- which(s==x)
    return(s[-idx][r])
  }
  
  # the '...M' variables denote the multiplicated data-set
  ftimeM <- rep(ftime, 3)
  fstatusM <- c(fstatus, rep(0,2*length(fstatus)))
  ftype <- ifelse(is.na(ftype), 0, ftype)
  ftypeM <- c(ftype, sapply(ftype, f, r=1), sapply(ftype, f, r=2))
  vaccineM <- rep(vaccine, 3)
  
  fit <- coxph(Surv(ftimeM, fstatusM) ~ vaccineM*ftypeM)
  return(summary(fit))
}
