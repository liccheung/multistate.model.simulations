#This program generates data under current US recommended screening intervals: 5 years following HPV-negative and 1 year following HPV-positive results.
#Output are true simulated precancer risk for initially HPV-negative, initially HPV-positive, and new HPV-positive

#increases number of infections to 3
#only 1 type
library(tidyverse)
library(survival)
library(Icens)

#General simulation settings
nsim <- 1 #number of simulation datasets 
n <- 10000000 #number of samples in each dataset
num <- nsim*n
set.seed(12726)

#Parameter values - set up so time represents years
p1 <- .1
p2 <- .5
p3 <- .055 #Didem's paper - 5.5% of HPV+ have prevalent CIN3+
#parameter for acquisition for HPV
lambda1 <- 0.055
#clearance parameters
shape21 <- 1 #0.702
scale21 <- 1.5
#progression parameters
shape22 <- 1
scale22 <- 60

#Part 1
#Initial state t0
h1 <- rbinom(num,1,(p1)) #binary where 1 signifies HPV+ individuals
h1[h1==1] <- 2*rbinom(length(h1[h1==1]),1,p3)+1 #for low risk types, 1=low risk 3=pre-cancer

#Part 2
#Times to acquisition, clearance, and progression allowing for reinfection
t1 <- rexp(num,lambda1) #time from HPV- to HPV+
#Time to t2 - Progression or Clearance
t21 <- rweibull(num,shape21,scale21)  #time to clearance
t22 <- rweibull(num,shape22,scale22) #time to progression
tlow <- pmin(t21,t22) #actual time is event that occurred first 

t12 <- rexp(num,lambda1) #time from HPV- to HPV+
#Time to t2 - Progression or Clearance
t212 <- rweibull(num,shape21,scale21)  #time to clearance
t222 <- rweibull(num,shape22,scale22) #time to progression
tlow2 <- pmin(t212,t222) #actual time is event that occurred first 

t13 <- rexp(num,lambda1) #time from HPV- to HPV+
#Time to t2 - Progression or Clearance
t213 <- rweibull(num,shape21,scale21)  #time to clearance
t223 <- rweibull(num,shape22,scale22) #time to progression
#LCC: corrected this - it should be min(t21,t22) if low-risk and min(t21,t23) if high risk
tlow3 <- pmin(t213,t223) #actual time is event that occurred first 

t14 <- rexp(num,lambda1) #time from HPV- to HPV+
#Time to t2 - Progression or Clearance
t214 <- rweibull(num,shape21,scale21)  #time to clearance
t224 <- rweibull(num,shape22,scale22) #time to progression
#LCC: corrected this - it should be min(t21,t22) if low-risk and min(t21,t23) if high risk
tlow4 <- pmin(t214,t224) #actual time is event that occurred first 

#Observation Process
#for flexibility, I assume visit_time is a vector
time_point <- function(visit_time){
  visit_fx <- ifelse(h1==3, 3, #pre-cancers cannot regress
               ifelse(h1==1 & tlow>visit_time, 1, #initially HPV+, no change       
                ifelse(h1==1 & tlow<=visit_time & tlow==t22, 3, #low-risk types progress
                  ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)>visit_time, 0, #low/high-risk types clear and is not reacquired
                   ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)>visit_time, 1, #reacquisition unresolved
                    ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t222, 3, #reacqusition progressed
                     ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)>visit_time, 0, #reacqisition cleared
                      ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)<=visit_time & (tlow+t12+tlow2+t13+tlow3)>visit_time, 1, #reacquisition 2x unresolved
                       ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)<=visit_time & (tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t223, 3, #reacqisition 2x progressed       
                         ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)<=visit_time & (tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (tlow+t12+tlow2+t13+tlow3+t14)>visit_time, 0, #reacqisition 2x cleared
                          ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)<=visit_time & (tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (tlow+t12+tlow2+t13+tlow3+t14)<=visit_time & (tlow+t12+tlow2+t13+tlow3+t14+tlow4)>visit_time, 1, #reacqisition 3x unresolved
                           ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)<=visit_time & (tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (tlow+t12+tlow2+t13+tlow3+t14)<=visit_time & (tlow+t12+tlow2+t13+tlow3+t14+tlow4)<=visit_time & tlow4==t224, 3, #reacqisition 3x unresolved
                            ifelse(h1==1 & tlow<=visit_time & tlow==t21 & (tlow+t12)<=visit_time & (tlow+t12+tlow2)<=visit_time & tlow2==t212 & (tlow+t12+tlow2+t13)<=visit_time & (tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (tlow+t12+tlow2+t13+tlow3+t14)<=visit_time & (tlow+t12+tlow2+t13+tlow3+t14+tlow4)<=visit_time & tlow4==t214, 0, #reacqisition 3x cleared, no more reacquisition allowed
                             ifelse(h1==0 & t1>visit_time, 0, #initially HPV-, stays HPV-
                              ifelse(h1==0 & t1<=visit_time & (t1+tlow)>visit_time, 1, #first HPV+, no change
                               ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t22, 3, #first HPV+ progress
                                ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)>visit_time, 0, #low/high-risk types clear and is not reacquired
                                 ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)>visit_time, 1, #reacquisition unresolved
                                  ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t222, 3, #reacqusition progressed
                                   ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)>visit_time, 0, #reacqisition cleared
                                    ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3)>visit_time, 1, #reacquisition 2x unresolved
                                     ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t223, 3, #reacqisition 2x progressed       
                                      ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (t1+tlow+t12+tlow2+t13+tlow3+t14)>visit_time, 0, #reacqisition 2x cleared
                                       ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (t1+tlow+t12+tlow2+t13+tlow3+t14)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3+t14+tlow4)>visit_time, 1, #reacqisition 3x unresolved
                                        ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (t1+tlow+t12+tlow2+t13+tlow3+t14)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3+t14+tlow4)<=visit_time & tlow4==t224, 3, #reacqisition 3x unresolved
                                         ifelse(h1==0 & t1<=visit_time & (t1+tlow)<=visit_time & tlow==t21 & (t1+tlow+t12)<=visit_time & (t1+tlow+t12+tlow2)<=visit_time & tlow2==t212 & (t1+tlow+t12+tlow2+t13)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3)<=visit_time & tlow3==t213 & (t1+tlow+t12+tlow2+t13+tlow3+t14)<=visit_time & (t1+tlow+t12+tlow2+t13+tlow3+t14+tlow4)<=visit_time & tlow4==t214, 0, NA)))))))))))))))))))))))))) #reacqisition 3x cleared, no more reacquisition allowed
  return(visit_fx)
}

#dependent visits - US recommendations: 5 years for HPV- result, 1 year for HPV+ results
r0v3 <- h1
v1_intervals <- ifelse(h1==0,rnorm(sum(h1==0),5,0.5),rnorm(sum(h1>0),1,0.25))
v1_intervals[v1_intervals<0] <- 0
v1v3 <- v1_intervals
r1v3 <- time_point(v1v3)

v2_intervals <- ifelse(r1v3==0,rnorm(sum(r1v3==0),5,0.5),rnorm(sum(r1v3>0),1,0.25))
v2_intervals[v2_intervals<0] <- 0
v2v3 <- v1v3 + v2_intervals
r2v3 <- time_point(v2v3)

v3_intervals <- ifelse(r2v3==0,rnorm(sum(r2v3==0),5,0.5),rnorm(sum(r2v3>0),1,0.25))
v3_intervals[v3_intervals<0] <- 0
v3v3 <- v2v3 + v3_intervals
r3v3 <- time_point(v3v3)

v4_intervals <- ifelse(r3v3==0,rnorm(sum(r3v3==0),5,0.5),rnorm(sum(r3v3>0),1,0.25))
v4_intervals[v4_intervals<0] <- 0
v4v3 <- v3v3 + v4_intervals
r4v3 <- time_point(v4v3)

v5_intervals <- ifelse(r4v3==0,rnorm(sum(r4v3==0),5,0.5),rnorm(sum(r4v3>0),1,0.25))
v5_intervals[v5_intervals<0] <- 0
v5v3 <- v4v3 + v5_intervals
r5v3 <- time_point(v5v3)

v6_intervals <- ifelse(r5v3==0,rnorm(sum(r5v3==0),5,0.5),rnorm(sum(r5v3>0),1,0.25))
v6_intervals[v6_intervals<0] <- 0
v6v3 <- v5v3 + v6_intervals
r6v3 <- time_point(v6v3)

v7_intervals <- ifelse(r6v3==0,rnorm(sum(r6v3==0),5,0.5),rnorm(sum(r6v3>0),1,0.25))
v7_intervals[v7_intervals<0] <- 0
v7v3 <- v6v3 + v7_intervals
r7v3 <- time_point(v7v3)

v8_intervals <- ifelse(r7v3==0,rnorm(sum(r7v3==0),5,0.5),rnorm(sum(r7v3>0),1,0.25))
v8_intervals[v8_intervals<0] <- 0
v8v3 <- v7v3 + v8_intervals
r8v3 <- time_point(v8v3)

v9_intervals <- ifelse(r8v3==0,rnorm(sum(r8v3==0),5,0.5),rnorm(sum(r8v3>0),1,0.25))
v9_intervals[v9_intervals<0] <- 0
v9v3 <- v8v3 + v9_intervals
r9v3 <- time_point(v9v3)

v10_intervals <- ifelse(r9v3==0,rnorm(sum(r9v3==0),5,0.5),rnorm(sum(r9v3>0),1,0.25))
v10_intervals[v10_intervals<0] <- 0
v10v3 <- v9v3 + v10_intervals
r10v3 <- time_point(v10v3)

#time of detection of new incident HPV
fpos <- ifelse(r0v3>0,NA,
               ifelse(r1v3>=1,v1v3,
                      ifelse(r2v3>=1,v2v3,
                             ifelse(r3v3>=1,v3v3,
                                    ifelse(r4v3>=1,v4v3,
                                           ifelse(r5v3>=1,v5v3,
                                                  ifelse(r6v3>=1,v6v3,
                                                         ifelse(r7v3>=1,v7v3,
                                                                ifelse(r8v3>=1,v8v3,
                                                                       ifelse(r9v3>=1,v9v3,
                                                                              ifelse(r10v3>=1,v10v3,NA)))))))))))

#true time to precancer
t_precancer <- ifelse(h1==3,0.001,
                ifelse(h1==1 & tlow==t22,t22,
                 ifelse(h1==1 & tlow==t21 & tlow2==t222,t21+t12+t222,
                  ifelse(h1==1 & tlow==t21 & tlow2==t212 & tlow3==t223,t21+t12+t212+t13+t223,
                   ifelse(h1==1 & tlow==t21 & tlow2==t212 & tlow3==t213 & tlow4==t224,t21+t12+t212+t13+t213+t14+t224,
                    ifelse(h1==0 & tlow==t22,t1+t22,
                     ifelse(h1==0 & tlow==t21 & tlow2==t222,t1+t21+t12+t222,
                      ifelse(h1==0 & tlow==t21 & tlow2==t212 & tlow3==t223,t1+t21+t12+t212+t13+t223,
                       ifelse(h1==0 & tlow==t21 & tlow2==t212 & tlow3==t213 & tlow4==t224,t1+t21+t12+t212+t13+t213+t14+t224,t1+t21+t12+t212+t13+t213+t14+t214)))))))))
precanc <- ifelse(t_precancer==(t1+t21+t12+t212+t13+t213+t14+t214),0,1) #indicator variable for precancer

#fit non-parametric survival models to true times to obtain cumulative precancer risk (accounting for prevalent precancers and right-censoring)
hpvpos <- coxph(Surv(t_precancer[h1>=1],precanc[h1>=1])~1)
hpvneg <- coxph(Surv(t_precancer[h1==0],precanc[h1==0])~1)
nhpvpos <- coxph(Surv(pmax((t_precancer-fpos),0.001)[is.na(fpos)==0],precanc[is.na(fpos)==0])~1)

hpvposfit <- survfit(hpvpos)
hpvnegfit <- survfit(hpvneg)
nhpvposfit <- survfit(nhpvpos)

save(hpvposfit,hpvnegfit,nhpvposfit,file="~/Desktop/BJC methods paper/targetedrisk.RData")