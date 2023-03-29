# Use simple power formula under a proportional hazards
# with a single treatment indicator covariate
# and event-driven design (Fleming and Harrington, 1991, page 395)
# Input parameters:
# n = total number of events during the 2-year follow-up period
# VE0, VE1: values of VE under the null and alternative
# alpha = 2-sided type I error rate of the test
# power = power required for rejecting the null if VE = VE1
# AR = Infection attack rate over the 2-year follow-up period
# CR = Censoring rate over the 2-year follow-up period
# Output:
# n, nv, np: Numbers of events for having the specified power
# N: Total sample size
# Assumptions:
# 1:1 randomization
# random censoring
# ~ rare event
# the analysis counts all infections between month 0 and 24
computeSS <- function(VE0,VE1,AR,CR,alpha,power) {
  num <- qnorm(1-alpha/2) + qnorm(power)
  den <- 0.5*(-log(1-VE1) + log(1-VE0))
  n <- ceiling((num/den)^2)
  np <- round(n/(2-VE1))
  nv <- n - np
  N <- ceiling((2*np)/(AR*(1-CR/2)))
  cat("\n")
  cat(paste("Required numbers of events n, nv, np = ",n,nv,np),"\n")
  cat("\n")
  cat(paste("Total sample size = ",N),"\n")
  vect <- c(n,nv,np,N,VE0,VE1,AR,CR,alpha,power)
  return(vect)
}
alpha <- 0.05
power <- 0.90
# censoring rate
CR <- .1
# Background designs null of zero
VE0 <- 0
VE1 <- 0.5
AR <- 0.03*2
ans4 <- computeSS(VE0,VE1,AR,CR,alpha,power)
AR <- 0.015*2
ans4 <- computeSS(VE0,VE1,AR,CR,alpha,power)
AR <- 0.005*2
ans4 <- computeSS(VE0,VE1,AR,CR,alpha,power)
