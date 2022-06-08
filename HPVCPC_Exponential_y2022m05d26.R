#This program generates data under fixed intervals, random intervals (mixed-cases interval censoring), and
#test-dependent (doctor's care) observation schemes.
#Precancer risk estimates are from healthy vs. precancer survival approach (EM-ICM to account for left-, interval-, and right-censoring)
#Precancer risk estimates are for initially HPV-negative, initially HPV-positive, and new HPV-positive

#increases number of infections to 3
#only 1 type
library(tidyverse)
library(survival)
library(Icens)

#General simulation settings
nsim <- 1 #number of simulation datasets 
n <- 10000 #number of samples in each dataset
num <- nsim*n
set.seed(12726)

#Parameter values - set up so time represents years
p1 <- .1
p2 <- .5
p3 <- .055 #Didem Egemen's paper - 5.5% of HPV+ have prevalent CIN3+
#parameter for acquisition of HPV
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
#Time to t1
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

#fixed intervals
r0 <- h1
r0.5 <- time_point(rep(0.5,num))
r1 <- time_point(rep(1,num))
r1.5 <- time_point(rep(1.5,num))
r2<- time_point(rep(2,num))
r2.5 <- time_point(rep(2.5,num))
r3 <- time_point(rep(3,num))
r3.5 <- time_point(rep(3.5,num))
r4 <- time_point(rep(4,num))
r4.5 <- time_point(rep(4.5,num))
r5 <- time_point(rep(5,num))
r5.5 <- time_point(rep(5.5,num))
r6 <- time_point(rep(6,num))
r6.5 <- time_point(rep(6.5,num))
r7 <- time_point(rep(7,num))
r7.5 <- time_point(rep(7.5,num))
r8 <- time_point(rep(8,num))
r8.5 <- time_point(rep(8.5,num))
r9 <- time_point(rep(9,num))
r9.5 <- time_point(rep(9.5,num))
r10 <- time_point(rep(10,num))

#numbers are more indicative of a moderately-high incidence setting
sum(r3>0 & r0==0)/sum(r0==0)  #7.3% incident infections after 3 years
sum(r0==0 & r3==3)/sum(r0==0 & r3>0)  #3.2% of incident HPV+ are CIN3+ vs. 2% reported in Didem's paper for well-screened KPNC cohort
sum(r0==0 & r3>0 & r8==3)/sum(r0==0 & r3>0)  #6.5% of incident HPV+ are CIN3+ after 5 years vs. 5.7% after 5-years reported in Didem's paper for KPNC Cohort
sum(r0==0 & r5==3)/sum(r0==0) #0.52% of HPV-negatives are CIN3+ after 5 years - this is similar to Li's paper on BD Onclarity trial (Didem's paper reports 0.15% for KPNC)


#mixed cases interval censoring
n_visits <- ifelse(rbinom(num,1,.3),11,sample.int(11, num, replace= TRUE))-1
intervals <- rnorm(sum(n_visits),1,0.25)
intervals[intervals<0] <- 0

a1 <- sum(n_visits>=1) #people who have 1 or more visits 
a2 <- sum(n_visits>=2) #people who have 2 or more visits
a3 <- sum(n_visits>=3)
a4 <- sum(n_visits>=4)
a5 <- sum(n_visits>=5)
a6 <- sum(n_visits>=6)
a7 <- sum(n_visits>=7)
a8 <- sum(n_visits>=8)
a9 <- sum(n_visits>=9)
a10 <- sum(n_visits>=10)

v1<- ifelse(n_visits>=1, intervals[1:a1],0) #for people with 1+ visits, interval between initial visit and 1st follow up 
v2 <- v1 + ifelse(n_visits>=2, intervals[(a1+1):(a1+a2)],0) #interval between visits for people with 1-2 follow ups
v3 <- v2 + ifelse(n_visits>=3, intervals[(a1+a2+1):(a1+a2+a3)],0)
v4 <- v3 + ifelse(n_visits>=4, intervals[(a1+a2+a3+1):(a1+a2+a3+a4)],0)
v5 <- v4 + ifelse(n_visits>=5, intervals[(a1+a2+a3+a4+1):(a1+a2+a3+a4+a5)],0)
v6 <- v5 + ifelse(n_visits>=6, intervals[(a1+a2+a3+a4+a5+1):(a1+a2+a3+a4+a5+a6)],0)
v7 <- v6 + ifelse(n_visits>=7, intervals[(a1+a2+a3+a4+a5+a6+1):(a1+a2+a3+a4+a5+a6+a7)],0)
v8 <- v7 + ifelse(n_visits>=8, intervals[(a1+a2+a3+a4+a5+a6+a7+1):(a1+a2+a3+a4+a5+a6+a7+a8)],0)
v9 <- v8 + ifelse(n_visits>=9, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9)],0)
v10 <- v9 + ifelse(n_visits>=10, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)],0)

r0v2 <- h1 #initial visit
r1v2 <- time_point(v1) #status of patients with 1+ visits at times listed in v1
r2v2 <- time_point(v2) #status of patients with 1-2 visits at times listed in v2
r3v2 <- time_point(v3)
r4v2 <- time_point(v4)
r5v2 <- time_point(v5)
r6v2 <- time_point(v6)
r7v2 <- time_point(v7)
r8v2 <- time_point(v8)
r9v2 <- time_point(v9)
r10v2 <- time_point(v10)

#dependent visits - doctor's care: 3 years after HPV-/1 year after HPV+
maxtime <- ifelse(rbinom(num,1,.3),12,runif(num,0,12))  #increased maxtime to 12
r0v3 <- h1
v1_intervals <- ifelse(h1==0,rnorm(sum(h1==0),3,0.5),rnorm(sum(h1>0),1,0.25))
v1_intervals[v1_intervals<0] <- 0
v1v3 <- ifelse(v1_intervals<=maxtime,v1_intervals,0)
nomore <- rep(0,num)
nomore <- ifelse(v1_intervals>maxtime | nomore==1,1,0)
r1v3 <- time_point(v1v3)

v2_intervals <- ifelse(r1v3==0,rnorm(sum(r1v3==0),3,0.5),rnorm(sum(r1v3>0),1,0.25))
v2_intervals[v2_intervals<0] <- 0
v2v3 <- ifelse((v1v3 + v2_intervals)<=maxtime,(v1v3 + v2_intervals),v1v3)
nomore <- ifelse((v1v3 + v2_intervals)>maxtime | nomore==1,1,0)
r2v3 <- time_point(v2v3)

v3_intervals <- ifelse(r2v3==0,rnorm(sum(r2v3==0),3,0.5),rnorm(sum(r2v3>0),1,0.25))
v3_intervals[v3_intervals<0] <- 0
v3v3 <- ifelse((v2v3 + v3_intervals)<=maxtime,(v2v3 + v3_intervals),v2v3)
nomore <- ifelse((v2v3 + v3_intervals)>maxtime | nomore==1,1,0)
r3v3 <- time_point(v3v3)

v4_intervals <- ifelse(r3v3==0,rnorm(sum(r3v3==0),3,0.5),rnorm(sum(r3v3>0),1,0.25))
v4_intervals[v4_intervals<0] <- 0
v4v3 <- ifelse((v3v3 + v4_intervals)<=maxtime,(v3v3 + v4_intervals),v3v3)
nomore <- ifelse((v3v3 + v4_intervals)>maxtime | nomore==1,1,0)
r4v3 <- time_point(v4v3)

v5_intervals <- ifelse(r4v3==0,rnorm(sum(r4v3==0),3,0.5),rnorm(sum(r4v3>0),1,0.25))
v5_intervals[v5_intervals<0] <- 0
v5v3 <- ifelse((v4v3 + v5_intervals)<=maxtime,(v4v3 + v5_intervals),v4v3)
nomore <- ifelse((v4v3 + v5_intervals)>maxtime | nomore==1,1,0)
r5v3 <- time_point(v5v3)

v6_intervals <- ifelse(r5v3==0,rnorm(sum(r5v3==0),3,0.5),rnorm(sum(r5v3>0),1,0.25))
v6_intervals[v6_intervals<0] <- 0
v6v3 <- ifelse((v5v3 + v6_intervals)<=maxtime,(v5v3 + v6_intervals),v5v3)
nomore <- ifelse((v5v3 + v6_intervals)>maxtime | nomore==1,1,0)
r6v3 <- time_point(v6v3)

v7_intervals <- ifelse(r6v3==0,rnorm(sum(r6v3==0),3,0.5),rnorm(sum(r6v3>0),1,0.25))
v7_intervals[v7_intervals<0] <- 0
v7v3 <- ifelse((v6v3 + v7_intervals)<=maxtime,(v6v3 + v7_intervals),v6v3)
nomore <- ifelse((v6v3 + v7_intervals)>maxtime | nomore==1,1,0)
r7v3 <- time_point(v7v3)

v8_intervals <- ifelse(r7v3==0,rnorm(sum(r7v3==0),3,0.5),rnorm(sum(r7v3>0),1,0.25))
v8_intervals[v8_intervals<0] <- 0
v8v3 <- ifelse((v7v3 + v8_intervals)<=maxtime,(v7v3 + v8_intervals),v7v3)
nomore <- ifelse((v7v3 + v8_intervals)>maxtime | nomore==1,1,0)
r8v3 <- time_point(v8v3)

v9_intervals <- ifelse(r8v3==0,rnorm(sum(r8v3==0),3,0.5),rnorm(sum(r8v3>0),1,0.25))
v9_intervals[v9_intervals<0] <- 0
v9v3 <- ifelse((v8v3 + v9_intervals)<=maxtime,(v8v3 + v9_intervals),v8v3)
nomore <- ifelse((v8v3 + v9_intervals)>maxtime | nomore==1,1,0)
r9v3 <- time_point(v9v3)

v10_intervals <- ifelse(r9v3==0,rnorm(sum(r9v3==0),3,0.5),rnorm(sum(r9v3>0),1,0.25))
v10_intervals[v10_intervals<0] <- 0
v10v3 <- ifelse((v9v3 + v10_intervals)<=maxtime,(v9v3 + v10_intervals),v9v3)
nomore <- ifelse((v9v3 + v10_intervals)>maxtime | nomore==1,1,0)
r10v3 <- time_point(v10v3)

#create data set of observed
obsdat1 <- data.frame(r0,r0.5,r1,r1.5,r2,r2.5,r3,r3.5,r4,r4.5,r5, r5.5, r6, r6.5, r7, r7.5, r8, r8.5, r9, r9.5, r10)
obsdat2 <- data.frame(r0v2, r1v2, r2v2, r3v2, r4v2, r5v2, r6v2, r7v2, r8v2, r9v2, r10v2, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
obsdat3 <- data.frame(r0v3, r1v3, r2v3, r3v3, r4v3, r5v3, r6v3, r7v3, r8v3, r9v3, r10v3, v1v3, v2v3, v3v3, v4v3, v5v3, v6v3, v7v3, v8v3, v9v3, v10v3)

#estimate for baseline HPV+, baseline HPV-, and new HPV+ (where the previous recorded visit was HPV-)
simsum <- data.frame()
for (j in 1:nsim){
  print(j)
  lower <- (j-1)*n+1
  upper <- j*n
  xsam <- data.frame(obsdat1[lower:upper,])
  
  #Fixed Observation Intervals Risk (Obs1)
  
  #risk for hpv- and hpv+ at enrollment in fixed intervals
  frisk_neg0 <- 0 #by definition
  frisk_pos0 <- sum(xsam$r0[xsam$r0>0]==3)/sum(xsam$r0>0)  #baseline HPV+, risk at 0; everyone who has precancer vs. everyone who's hpv+ at r0
  
  frisk_neg0.5 <- sum(xsam$r0.5[xsam$r0==0]==3)/sum(xsam$r0==0)  #baseline HPV-, risk at 0.5; everyone who developed precancer at r0.5 after being hpv- at r0 compared to everyone who is hpv- at r0
  frisk_pos0.5 <- sum(xsam$r0.5[xsam$r0>0]==3)/sum(xsam$r0>0)  #baseline HPV+, risk at 0.5; everyone who develops pre-cancer at r0.5 from hpv+ at r0
  
  frisk_neg1 <- sum(xsam$r1[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 1
  frisk_pos1 <- sum(xsam$r1[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 1
  
  frisk_neg1.5 <- sum(xsam$r1.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 1.5
  frisk_pos1.5 <- sum(xsam$r1.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 1.5
  
  frisk_neg2 <- sum(xsam$r2[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 2
  frisk_pos2 <- sum(xsam$r2[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 2
  
  frisk_neg2.5 <- sum(xsam$r2.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 2.5
  frisk_pos2.5 <- sum(xsam$r2.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 2.5
  
  frisk_neg3 <- sum(xsam$r3[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 3
  frisk_pos3 <- sum(xsam$r3[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 3
  
  frisk_neg3.5 <- sum(xsam$r3.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 3.5
  frisk_pos3.5 <- sum(xsam$r3.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 3.5
  
  frisk_neg4 <- sum(xsam$r4[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 4
  frisk_pos4 <- sum(xsam$r4[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 4
  
  frisk_neg4.5 <- sum(xsam$r4.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 4.5
  frisk_pos4.5 <- sum(xsam$r4.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 4.5
  
  frisk_neg5 <- sum(xsam$r5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos5 <- sum(xsam$r5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg5.5 <- sum(xsam$r5.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos5.5 <- sum(xsam$r5.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg6 <- sum(xsam$r6[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos6 <- sum(xsam$r6[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg6.5 <- sum(xsam$r6.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos6.5 <- sum(xsam$r6.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg7 <- sum(xsam$r7[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos7 <- sum(xsam$r7[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg7.5 <- sum(xsam$r7.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos7.5 <- sum(xsam$r7.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg8 <- sum(xsam$r8[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos8 <- sum(xsam$r8[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg8.5 <- sum(xsam$r8.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos8.5 <- sum(xsam$r8.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg9 <- sum(xsam$r9[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos9 <- sum(xsam$r9[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg9.5 <- sum(xsam$r9.5[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos9.5 <- sum(xsam$r9.5[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  frisk_neg10 <- sum(xsam$r10[xsam$r0==0]==3)/sum(xsam$r0==0) #baseline HPV-, risk at 5
  frisk_pos10 <- sum(xsam$r10[xsam$r0>0]==3)/sum(xsam$r0>0) #baseline HPV+, risk at 5
  
  
  
  #Fixed New HPV
  newhpv_t <- ifelse(xsam$r0>0,0, #these were HPV+ at start
                     ifelse(xsam$r0.5>0,0.5,
                            ifelse(xsam$r1>0,1,
                                   ifelse(xsam$r1.5>0,1.5,
                                          ifelse(xsam$r2>0,2,
                                                 ifelse(xsam$r2.5>0,2.5,
                                                        ifelse(xsam$r3>0,3,
                                                               ifelse(xsam$r3.5>0,3.5,
                                                                      ifelse(xsam$r4>0,4,
                                                                             ifelse(xsam$r5>0,5,
                                                                                    ifelse(xsam$r5.5>0,5.5,
                                                                                           ifelse(xsam$r6>0,6,
                                                                                                  ifelse(xsam$r6.5>0,6.5,
                                                                                                         ifelse(xsam$r7>0,7,
                                                                                                                ifelse(xsam$r7.5>0,7.5,
                                                                                                                       ifelse(xsam$r8>0,8,
                                                                                                                              ifelse(xsam$r8.5>0,8.5,
                                                                                                                                     ifelse(xsam$r9>0,9,
                                                                                                                                            ifelse(xsam$r9.5>0,9.5,
                                                                                                                                                   ifelse(xsam$r10>0,10,NA)))))))))))))))))))) #these were never HPV+
  
  
  #last time before precancer detection
  fnlower <-  ifelse(is.na(newhpv_t)==1 | newhpv_t==0,NA,  #these were never HPV+ or were HPV+ at time 0; make all NA's (HPV-)and 0's (HPV+ at r0) NA
                     ifelse(xsam$r10 < 3 & newhpv_t <= 10, 10-newhpv_t, #if patient tested HPV+ at year 10 and they first tested HPV+ at or before year 10, time between first HPV+ test and their latest visit
                            ifelse(xsam$r9.5 < 3 & newhpv_t <= 9.5, 9.5-newhpv_t,
                                   ifelse(xsam$r9 < 3 & newhpv_t <= 9, 9-newhpv_t,
                                          ifelse(xsam$r8.5 < 3 & newhpv_t <= 8.5, 8.5-newhpv_t,
                                                 ifelse(xsam$r8 < 3 & newhpv_t <= 8, 8-newhpv_t,
                                                        ifelse(xsam$r7.5 < 3 & newhpv_t <= 7.5, 7.5-newhpv_t,
                                                               ifelse(xsam$r7 < 3 & newhpv_t <= 7, 7-newhpv_t,
                                                                      ifelse(xsam$r6.5 < 3 & newhpv_t <= 6.5, 6.5-newhpv_t,
                                                                             ifelse(xsam$r6 < 3 & newhpv_t <= 6, 6-newhpv_t,
                                                                                    ifelse(xsam$r5.5 < 3 & newhpv_t <= 5.5, 5.5-newhpv_t,
                                                                                           ifelse(xsam$r5 < 3 & newhpv_t <= 5, 5-newhpv_t,
                                                                                                  ifelse(xsam$r4.5 < 3 & newhpv_t <= 4.5, 4.5-newhpv_t,
                                                                                                         ifelse(xsam$r4 < 3 & newhpv_t <= 4, 4-newhpv_t,
                                                                                                                ifelse(xsam$r3.5 < 3 & newhpv_t <= 3.5, 3.5-newhpv_t,
                                                                                                                       ifelse(xsam$r3 < 3 & newhpv_t <= 3, 3-newhpv_t,
                                                                                                                              ifelse(xsam$r2.5 < 3 & newhpv_t <= 2.5, 2.5-newhpv_t,
                                                                                                                                     ifelse(xsam$r2 < 3 & newhpv_t <= 2, 2-newhpv_t,
                                                                                                                                            ifelse(xsam$r1.5 < 3 & newhpv_t <= 1.5, 1.5-newhpv_t,
                                                                                                                                                   ifelse(xsam$r1 < 3 & newhpv_t <= 1, 1-newhpv_t,
                                                                                                                                                          ifelse(xsam$r0.5 < 3 & newhpv_t <= 0.5, 0.5-newhpv_t,
                                                                                                                                                                 ifelse(xsam$r0 < 3 & newhpv_t <= 0, 0, -0.01)))))))))))))))))))))) 

  #time of precancer detection
  fnupper <- ifelse(is.na(newhpv_t)==1 | newhpv_t==0,NA, #these were never HPV+ or were HPV+ at time 0
                    ifelse(xsam$r0.5==3 & newhpv_t<=0.5,0.5-newhpv_t,
                           ifelse(xsam$r1==3 & newhpv_t<=1,1-newhpv_t,
                                  ifelse(xsam$r1.5==3 & newhpv_t<=1.5,1.5-newhpv_t,
                                         ifelse(xsam$r2==3 & newhpv_t<=2,2-newhpv_t,
                                                ifelse(xsam$r2.5==3 & newhpv_t<=2.5,2.5-newhpv_t,
                                                       ifelse(xsam$r3==3 & newhpv_t<=3,3-newhpv_t,
                                                              ifelse(xsam$r3.5==3 & newhpv_t<=3.5,3.5-newhpv_t,
                                                                     ifelse(xsam$r4==3 & newhpv_t<=4,4-newhpv_t,
                                                                            ifelse(xsam$r4.5==3 & newhpv_t<=4.5,4.5-newhpv_t,
                                                                                   ifelse(xsam$r5==3 & newhpv_t<=5,5-newhpv_t,
                                                                                          ifelse(xsam$r5.5==3 & newhpv_t<=5.5,5.5-newhpv_t,
                                                                                                 ifelse(xsam$r6==3 & newhpv_t<=6,6-newhpv_t,
                                                                                                        ifelse(xsam$r6.5==3 & newhpv_t<=6.5,6.5-newhpv_t,
                                                                                                               ifelse(xsam$r7==3 & newhpv_t<=7,7-newhpv_t,
                                                                                                                      ifelse(xsam$r7.5==3 & newhpv_t<=7.5,7.5-newhpv_t,
                                                                                                                             ifelse(xsam$r8==3 & newhpv_t<=8,8-newhpv_t,
                                                                                                                                    ifelse(xsam$r8.5==3 & newhpv_t<=8.5,8.5-newhpv_t,
                                                                                                                                           ifelse(xsam$r9==3 & newhpv_t<=9,9-newhpv_t,
                                                                                                                                                  ifelse(xsam$r9.5==3 & newhpv_t<=9.5,9.5-newhpv_t,
                                                                                                                                                         ifelse(xsam$r10==3 & newhpv_t<=10,10-newhpv_t,Inf)))))))))))))))))))))
  
  #fit EM-ICM model to interval-censored data, with adjustments to allow for prevalent precancer estimation
  A <- cbind(fnlower[is.na(fnlower)==0], fnupper[is.na(fnupper)==0])
  fit <- EMICM(A)
  
  #risk of new hpv+ at fixed intervals
  fnrisk_0 <- fit$sigma[min(which(fit$intmap[2,]>0))] #estimated risk at t0; select columns where bottom row of intmap (time)>0, min takes earliest occurrence of that
  fnrisk_0.5 <- fit$sigma[min(which(fit$intmap[2,]>0.5))]
  fnrisk_1 <- fit$sigma[min(which(fit$intmap[2,]>1))]
  fnrisk_1.5 <- fit$sigma[min(which(fit$intmap[2,]>1.5))]
  fnrisk_2 <- fit$sigma[min(which(fit$intmap[2,]>2))]
  fnrisk_2.5 <- fit$sigma[min(which(fit$intmap[2,]>2.5))]
  fnrisk_3 <- fit$sigma[min(which(fit$intmap[2,]>3))]
  fnrisk_3.5 <- fit$sigma[min(which(fit$intmap[2,]>3.5))]
  fnrisk_4 <- fit$sigma[min(which(fit$intmap[2,]>4))]
  fnrisk_4.5 <- fit$sigma[min(which(fit$intmap[2,]>4.5))]
  fnrisk_5 <- fit$sigma[min(which(fit$intmap[2,]>5))]
  fnrisk_5.5 <- fit$sigma[min(which(fit$intmap[2,]>5.5))]
  fnrisk_6 <- fit$sigma[min(which(fit$intmap[2,]>6))]
  fnrisk_6.5 <- fit$sigma[min(which(fit$intmap[2,]>6))]
  fnrisk_7 <- fit$sigma[min(which(fit$intmap[2,]>7))]
  fnrisk_7.5 <- fit$sigma[min(which(fit$intmap[2,]>7.5))]
  fnrisk_8 <- fit$sigma[min(which(fit$intmap[2,]>8))]
  fnrisk_8.5 <- fit$sigma[min(which(fit$intmap[2,]>8.5))]
  fnrisk_9 <- fit$sigma[min(which(fit$intmap[2,]>9))]
  fnrisk_9.5 <- fit$sigma[min(which(fit$intmap[2,]>9.5))]
  fnrisk_10 <- fit$sigma[min(which(fit$intmap[2,]>10))]
    
  #Mixed Case Intervals Risk (Obs2)
  xsam2 <- data.frame(obsdat2[lower:upper,])
  
  #HPV- at enrollment
  lower2 <- ifelse(xsam2$r10v2<3,xsam2$v10,
                   ifelse(xsam2$r9v2<3,xsam2$v9, 
                          ifelse(xsam2$r8v2<3,xsam2$v8,
                                 ifelse(xsam2$r7v2<3,xsam2$v7,
                                        ifelse(xsam2$r6v2<3,xsam2$v6,
                                               ifelse(xsam2$r5v2<3,xsam2$v5,
                                                      ifelse(xsam2$r4v2<3,xsam2$v4,
                                                             ifelse(xsam2$r3v2<3,xsam2$v3,
                                                                    ifelse(xsam2$r2v2<3,xsam2$v2,
                                                                           ifelse(xsam2$r1v2<3,xsam2$v1,
                                                                                  ifelse(xsam2$r0v2<3,0,-0.01)))))))))))
  #upper is the first time = 3
  upper2 <- ifelse(xsam2$r0v2==3,0,
                   ifelse(xsam2$r1v2==3,xsam2$v1, 
                          ifelse(xsam2$r2v2==3,xsam2$v2,
                                 ifelse(xsam2$r3v2==3,xsam2$v3,
                                        ifelse(xsam2$r4v2==3,xsam2$v4,
                                               ifelse(xsam2$r5v2==3,xsam2$v5,
                                                      ifelse(xsam2$r6v2==3,xsam2$v6,
                                                             ifelse(xsam2$r7v2==3,xsam2$v7,
                                                                    ifelse(xsam2$r8v2==3,xsam2$v8,
                                                                           ifelse(xsam2$r9v2==3,xsam2$v9,
                                                                                  ifelse(xsam2$r10v2==3,xsam2$v10,Inf)))))))))))
  B <- cbind(lower2[xsam2$r0v2==0],upper2[xsam2$r0v2==0])
  mrisk_neg <- EMICM(B)
  
  C <- cbind(lower2[xsam2$r0v2>0], upper2[xsam2$r0v2>0])
  mrisk_pos <- EMICM(C)
  
  #initially HPV- and HPV+ cumulative risk at mixed intervals
  mrisk_neg0 <- ifelse(length(which(mrisk_neg$intmap[2,]<=0))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=0))])
  mrisk_pos0 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=0))] #estimated risk at t0; select columns where bottom row of intmap (time)>0, min takes earliest occurrence of that
  
  mrisk_neg0.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=0.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=0.5))])
  mrisk_pos0.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=0.5))]
  
  mrisk_neg1 <- ifelse(length(which(mrisk_neg$intmap[2,]<=1))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=1))])
  mrisk_pos1 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=1))]
  
  mrisk_neg1.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=1.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=1.5))])
  mrisk_pos1.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=1.5))]
  
  mrisk_neg2 <- ifelse(length(which(mrisk_neg$intmap[2,]<=2))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=2))])
  mrisk_pos2 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=2))]
  
  mrisk_neg2.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=2.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=2.5))])
  mrisk_pos2.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=2.5))]
  
  mrisk_neg3 <- ifelse(length(which(mrisk_neg$intmap[2,]<=3))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=3))])
  mrisk_pos3 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=3))]
  
  mrisk_neg3.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=3.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=3.5))])
  mrisk_pos3.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=3.5))]
  
  mrisk_neg4 <- ifelse(length(which(mrisk_neg$intmap[2,]<=4))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=4))])
  mrisk_pos4 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=4))]
  
  mrisk_neg4.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=4.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=4.5))])
  mrisk_pos4.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=4.5))]
  
  mrisk_neg5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=5))])
  mrisk_pos5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=5))]
  
  mrisk_neg5.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=5.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=5.5))])
  mrisk_pos5.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=5.5))]
  
  mrisk_neg6 <- ifelse(length(which(mrisk_neg$intmap[2,]<=6))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=6))])
  mrisk_pos6 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=6))]
  
  mrisk_neg6.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=6.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=6.5))])
  mrisk_pos6.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=6.5))]
  
  mrisk_neg7 <- ifelse(length(which(mrisk_neg$intmap[2,]<=7))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=7))])
  mrisk_pos7 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=7))]
  
  mrisk_neg7.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=7.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=7.5))])
  mrisk_pos7.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=7.5))]
  
  mrisk_neg8 <- ifelse(length(which(mrisk_neg$intmap[2,]<=8))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=8))])
  mrisk_pos8 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=8))]
  
  mrisk_neg8.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=8.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=8.5))])
  mrisk_pos8.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=8.5))]
  
  mrisk_neg9 <- ifelse(length(which(mrisk_neg$intmap[2,]<=9))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=9))])
  mrisk_pos9 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=9))]
  
  mrisk_neg9.5 <- ifelse(length(which(mrisk_neg$intmap[2,]<=9.5))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=9.5))])
  mrisk_pos9.5 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=9.5))]
  
  mrisk_neg10 <- ifelse(length(which(mrisk_neg$intmap[2,]<=10))==0,0,mrisk_neg$sigma[max(which(mrisk_neg$intmap[2,]<=10))])
  mrisk_pos10 <- mrisk_pos$sigma[max(which(mrisk_pos$intmap[2,]<=10))]
  
  #Mixed-cases Newly HPV+
  mnewhpv_t <- ifelse(xsam2$r0v2>0,0, #these were HPV+ at start
                      ifelse(xsam2$r1v2>0,xsam2$v1,
                             ifelse(xsam2$r2v2>0,xsam2$v2,
                                    ifelse(xsam2$r3v2>0,xsam2$v3,
                                           ifelse(xsam2$r4v2>0,xsam2$v4,
                                                  ifelse(xsam2$r5v2>0,xsam2$v5,
                                                         ifelse(xsam2$r6v2>0,xsam2$v6,
                                                                ifelse(xsam2$r7v2>0,xsam2$v7,
                                                                       ifelse(xsam2$r8v2>0,xsam2$v8,
                                                                              ifelse(xsam2$r9v2>0,xsam2$v9,
                                                                                     ifelse(xsam2$r10v2>0,xsam2$v10,NA))))))))))) #these were never HPV+
  
  mnlower <-  ifelse(is.na(mnewhpv_t)==1 | mnewhpv_t==0,NA,  #these were never HPV+ or were HPV+ at time 0; make all NA's (HPV-)and 0's (HPV+ at r0) NA
                     ifelse(xsam2$r10v2 < 3 & mnewhpv_t <= xsam2$v10, xsam2$v10-mnewhpv_t, #if patient tested HPV+ at year 10 and they first tested HPV+ at or before year 10, time between first HPV+ test and their latest visit
                            ifelse(xsam2$r9v2 < 3 & mnewhpv_t <= xsam2$v9, xsam2$v9-mnewhpv_t,
                                   ifelse(xsam2$r8v2 < 3 & mnewhpv_t <= xsam2$v8, xsam2$v8-mnewhpv_t,
                                          ifelse(xsam2$r7v2 < 3 & mnewhpv_t <= xsam2$v7, xsam2$v7-mnewhpv_t,
                                                 ifelse(xsam2$r6v2 < 3 & mnewhpv_t <= xsam2$v6, xsam2$v6-mnewhpv_t,
                                                        ifelse(xsam2$r5v2 < 3 & mnewhpv_t <= xsam2$v5, xsam2$v5-mnewhpv_t,
                                                               ifelse(xsam2$r4v2 < 3 & mnewhpv_t <= xsam2$v4, xsam2$v4-mnewhpv_t,
                                                                      ifelse(xsam2$r3v2 < 3 & mnewhpv_t <= xsam2$v3, xsam2$v3-mnewhpv_t,
                                                                             ifelse(xsam2$r2v2 < 3 & mnewhpv_t <= xsam2$v2, xsam2$v2-mnewhpv_t,
                                                                                    ifelse(xsam2$r1v2 < 3 & mnewhpv_t <= xsam2$v1, xsam2$v1-mnewhpv_t,
                                                                                           ifelse(xsam2$r0v2 < 3 & mnewhpv_t <= 0, 0, -0.01)))))))))))) 
  
  mnupper <- ifelse(is.na(mnewhpv_t)==1 | mnewhpv_t==0,NA, #these were never HPV+ or were HPV+ at time 0
                    ifelse(xsam2$r1v2==3 & mnewhpv_t<=xsam2$v1,xsam2$v1-mnewhpv_t,
                           ifelse(xsam2$r2v2==3 & mnewhpv_t<=xsam2$v2,xsam2$v2-mnewhpv_t,
                                  ifelse(xsam2$r3v2==3 & mnewhpv_t<=xsam2$v3,xsam2$v3-mnewhpv_t,
                                         ifelse(xsam2$r4v2==3 & mnewhpv_t<=xsam2$v4,xsam2$v4-mnewhpv_t,
                                                ifelse(xsam2$r5v2==3 & mnewhpv_t<=xsam2$v5,xsam2$v5-mnewhpv_t,
                                                       ifelse(xsam2$r6v2==3 & mnewhpv_t<=xsam2$v6,xsam2$v6-mnewhpv_t,
                                                              ifelse(xsam2$r7v2==3 & mnewhpv_t<=xsam2$v7,xsam2$v7-mnewhpv_t,
                                                                     ifelse(xsam2$r8v2==3 & mnewhpv_t<=xsam2$v8,xsam2$v8-mnewhpv_t,
                                                                            ifelse(xsam2$r9v2==3 & mnewhpv_t<=xsam2$v9,xsam2$v9-mnewhpv_t,
                                                                                   ifelse(xsam2$r10v2==3 & mnewhpv_t<=xsam2$v10,xsam2$v10-mnewhpv_t,Inf)))))))))))
  
  E <- cbind(mnlower[is.na(mnlower)==0], mnupper[is.na(mnupper)==0])
  mnfit <- EMICM(E)
  
  mnrisk_0 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=0))]
  mnrisk_0.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=0.5))]
  mnrisk_1 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=1))]
  mnrisk_1.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=1.5))]
  mnrisk_2 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=2))]
  mnrisk_2.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=2.5))]
  mnrisk_3 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=3))]
  mnrisk_3.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=3.5))]
  mnrisk_4 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=4))]
  mnrisk_4.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=4.5))]
  mnrisk_5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=5))]
  mnrisk_5.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=5.5))]
  mnrisk_6 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=6))]
  mnrisk_6.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=6))]
  mnrisk_7 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=7))]
  mnrisk_7.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=7.5))]
  mnrisk_8 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=8))]
  mnrisk_8.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=8.5))]
  mnrisk_9 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=9))]
  mnrisk_9.5 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=9.5))]
  mnrisk_10 <- mnfit$sigma[max(which(mnfit$intmap[2,]<=10))]
  
  #Dependent Observation Intervals - Doctor's Care - Risk (Obs3)
  
  xsam3 <- data.frame(obsdat3[lower:upper,])
  
  #HPV- at enrollment
  lower3 <- ifelse(xsam3$r10v3<3,xsam3$v10v3,
                   ifelse(xsam3$r9v3<3,xsam3$v9v3, 
                          ifelse(xsam3$r8v3<3,xsam3$v8v3,
                                 ifelse(xsam3$r7v3<3,xsam3$v7v3,
                                        ifelse(xsam3$r6v3<3,xsam3$v6v3,
                                               ifelse(xsam3$r5v3<3,xsam3$v5v3,
                                                      ifelse(xsam3$r4v3<3,xsam3$v4v3,
                                                             ifelse(xsam3$r3v3<3,xsam3$v3v3,
                                                                    ifelse(xsam3$r2v3<3,xsam3$v3,
                                                                           ifelse(xsam3$r1v3<3,xsam3$v1v3,
                                                                                  ifelse(xsam3$r0v3<3,0,-0.01)))))))))))
  #upper is the first time = 3
  upper3 <- ifelse(xsam3$r0v3==3,0,
                   ifelse(xsam3$r1v3==3,xsam3$v1v3, 
                          ifelse(xsam3$r2v3==3,xsam3$v2v3,
                                 ifelse(xsam3$r3v3==3,xsam3$v3v3,
                                        ifelse(xsam3$r4v3==3,xsam3$v4v3,
                                               ifelse(xsam3$r5v3==3,xsam3$v5v3,
                                                      ifelse(xsam3$r6v3==3,xsam3$v6v3,
                                                             ifelse(xsam3$r7v3==3,xsam3$v7v3,
                                                                    ifelse(xsam3$r8v3==3,xsam3$v8v3,
                                                                           ifelse(xsam3$r9v3==3,xsam3$v9v3,
                                                                                  ifelse(xsam3$r10v3==3,xsam3$v10v3,Inf)))))))))))
  
  F <- cbind(lower3[xsam3$r0v3==0], upper3[xsam3$r0v3==0])  
  drisk_neg <- EMICM(F)  
  
  G <- cbind(lower3[xsam3$r0v3>0], upper3[xsam3$r0v3>0])  
  drisk_pos <- EMICM(G)  
  
  drisk_neg0 <- ifelse(length(which(drisk_neg$intmap[2,]<=0))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=0))])
  drisk_pos0 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=0))] #estimated risk at t0; select columns where bottom row of intmap (time)>0, min takes earliest occurrence of that
  
  drisk_neg0.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=0.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=0.5))])
  drisk_pos0.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=0.5))]
  
  drisk_neg1 <- ifelse(length(which(drisk_neg$intmap[2,]<=1))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=1))])
  drisk_pos1 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=1))]
  
  drisk_neg1.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=1.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=1.5))])
  drisk_pos1.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=1.5))]
  
  drisk_neg2 <- ifelse(length(which(drisk_neg$intmap[2,]<=2))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=2))])
  drisk_pos2 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=2))]
  
  drisk_neg2.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=2.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=2.5))])
  drisk_pos2.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=2.5))]
  
  drisk_neg3 <- ifelse(length(which(drisk_neg$intmap[2,]<=3))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=3))])
  drisk_pos3 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=3))]
  
  drisk_neg3.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=3.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=3.5))])
  drisk_pos3.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=3.5))]
  
  drisk_neg4 <- ifelse(length(which(drisk_neg$intmap[2,]<=4))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=4))])
  drisk_pos4 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=4))]
  
  drisk_neg4.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=4.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=4.5))])
  drisk_pos4.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=4.5))]
  
  drisk_neg5 <- ifelse(length(which(drisk_neg$intmap[2,]<=5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=5))])
  drisk_pos5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=5))]
  
  drisk_neg5.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=5.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=5.5))])
  drisk_pos5.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=5.5))]
  
  drisk_neg6 <- ifelse(length(which(drisk_neg$intmap[2,]<=6))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=6))])
  drisk_pos6 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=6))]
  
  drisk_neg6.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=6.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=6.5))])
  drisk_pos6.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=6.5))]
  
  drisk_neg7 <- ifelse(length(which(drisk_neg$intmap[2,]<=7))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=7))])
  drisk_pos7 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=7))]
  
  drisk_neg7.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=7.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=7.5))])
  drisk_pos7.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=7.5))]
  
  drisk_neg8 <- ifelse(length(which(drisk_neg$intmap[2,]<=8))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=8))])
  drisk_pos8 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=8))]
  
  drisk_neg8.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=8.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=8.5))])
  drisk_pos8.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=8.5))]
  
  drisk_neg9 <- ifelse(length(which(drisk_neg$intmap[2,]<=9))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=9))])
  drisk_pos9 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=9))]
  
  drisk_neg9.5 <- ifelse(length(which(drisk_neg$intmap[2,]<=9.5))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=9.5))])
  drisk_pos9.5 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=9.5))]
  
  drisk_neg10 <- ifelse(length(which(drisk_neg$intmap[2,]<=10))==0,0,drisk_neg$sigma[max(which(drisk_neg$intmap[2,]<=10))])
  drisk_pos10 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=10))]
  
  
  #Dependent New HPV
  dnewhpv_t <- ifelse(xsam3$r0v3>0,0, #these were HPV+ at start
                      ifelse(xsam3$r1v3>0,xsam3$v1v3,
                             ifelse(xsam3$r2v3>0,xsam3$v2v3,
                                    ifelse(xsam3$r3v3>0,xsam3$v3v3,
                                           ifelse(xsam3$r4v3>0,xsam3$v4v3,
                                                  ifelse(xsam3$r5v3>0,xsam3$v5v3,
                                                         ifelse(xsam3$r6v3>0,xsam3$v6v3,
                                                                ifelse(xsam3$r7v3>0,xsam3$v7v3,
                                                                       ifelse(xsam3$r8v3>0,xsam3$v8v3,
                                                                              ifelse(xsam3$r9v3>0,xsam3$v9v3,
                                                                                     ifelse(xsam3$r10v3>0,xsam3$v10v3,NA))))))))))) #these were never HPV+
  
  
  #Then code last time before precancer detection
  dnlower <-  ifelse(is.na(dnewhpv_t)==1 | dnewhpv_t==0,NA,  #these were never HPV+ or were HPV+ at time 0; make all NA's (HPV-)and 0's (HPV+ at r0) NA
                     ifelse(xsam3$r10v3 < 3 & dnewhpv_t <= xsam3$v10v3, xsam3$v10v3-dnewhpv_t, #if patient tested HPV+ at year 10 and they first tested HPV+ at or before year 10, time between first HPV+ test and their latest visit
                            ifelse(xsam3$r9v3 < 3 & dnewhpv_t <= xsam3$v9v3, xsam3$v9v3-dnewhpv_t,
                                   ifelse(xsam3$r8v3 < 3 & dnewhpv_t <= xsam3$v8v3, xsam3$v8v3-dnewhpv_t,
                                          ifelse(xsam3$r7v3 < 3 & dnewhpv_t <= xsam3$v7v3, xsam3$v7v3-dnewhpv_t,
                                                 ifelse(xsam3$r6v3 < 3 & dnewhpv_t <= xsam3$v6v3, xsam3$v6v3-dnewhpv_t,
                                                        ifelse(xsam3$r5v3 < 3 & dnewhpv_t <= xsam3$v5v3, xsam3$v5v3-dnewhpv_t,
                                                               ifelse(xsam3$r4v3 < 3 & dnewhpv_t <= xsam3$v4v3, xsam3$v4v3-dnewhpv_t,
                                                                      ifelse(xsam3$r3v3 < 3 & dnewhpv_t <= xsam3$v3v3, xsam3$v3v3-dnewhpv_t,
                                                                             ifelse(xsam3$r2v3 < 3 & dnewhpv_t <= xsam3$v2v3, xsam3$v2v3-dnewhpv_t,
                                                                                    ifelse(xsam3$r1v3 < 3 & dnewhpv_t <= xsam3$v1v3, xsam3$v1v3-dnewhpv_t,
                                                                                           ifelse(xsam3$r0v3 < 3 & dnewhpv_t <= 0, 0, -0.01))))))))))))
  
  #Then code time of precancer detection
  dnupper <- ifelse(is.na(dnewhpv_t)==1 | dnewhpv_t==0,NA, #these were never HPV+ or were HPV+ at time 0
                    ifelse(xsam3$r1v3==3 & dnewhpv_t<=xsam3$v1v3,xsam3$v1v3-dnewhpv_t,
                           ifelse(xsam3$r2v3==3 & dnewhpv_t<=xsam3$v2v3,xsam3$v2v3-dnewhpv_t,
                                  ifelse(xsam3$r3v3==3 & dnewhpv_t<=xsam3$v3v3,xsam3$v3v3-dnewhpv_t,
                                         ifelse(xsam3$r4v3==3 & dnewhpv_t<=xsam3$v4v3,xsam3$v4v3-dnewhpv_t,
                                                ifelse(xsam3$r5v3==3 & dnewhpv_t<=xsam3$v5v3,xsam3$v5v3-dnewhpv_t,
                                                       ifelse(xsam3$r6v3==3 & dnewhpv_t<=xsam3$v6v3,xsam3$v6v3-dnewhpv_t,
                                                              ifelse(xsam3$r7v3==3 & dnewhpv_t<=xsam3$v7v3,xsam3$v7v3-dnewhpv_t,
                                                                     ifelse(xsam3$r8v3==3 & dnewhpv_t<=xsam3$v8v3,xsam3$v8v3-dnewhpv_t,
                                                                            ifelse(xsam3$r9v3==3 & dnewhpv_t<=xsam3$v9v3,xsam3$v9v3-dnewhpv_t,
                                                                                   ifelse(xsam3$r10v3==3 & dnewhpv_t<=xsam3$v10v3,xsam3$v10v3-dnewhpv_t,Inf)))))))))))
  
  H <- cbind(dnlower[is.na(dnlower)==0], dnupper[is.na(dnupper)==0])
  dnfit <- EMICM(H)
  
  #risk for new hpv+ at depedent intervals
  dnrisk_0 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=0))] #estimated risk at t0; select columns where bottom row of intmap (time)>0, min takes earliest occurrence of that
  dnrisk_0.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=0.5))]
  dnrisk_1 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=1))]
  dnrisk_1.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=1.5))]
  dnrisk_2 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=2))]
  dnrisk_2.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=2.5))]
  dnrisk_3 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=3))]
  dnrisk_3.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=3.5))]
  dnrisk_4 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=4))]
  dnrisk_4.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=4.5))]
  dnrisk_5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=5))]
  dnrisk_5.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=5.5))]
  dnrisk_6 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=6))]
  dnrisk_6.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=6))]
  dnrisk_7 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=7))]
  dnrisk_7.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=7.5))]
  dnrisk_8 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=8))]
  dnrisk_8.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=8.5))]
  dnrisk_9 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=9))]
  dnrisk_9.5 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=9.5))]
  dnrisk_10 <- dnfit$sigma[max(which(dnfit$intmap[2,]<=10))]
  
  saveres <- data.frame(frisk_neg0,frisk_neg0.5,frisk_neg1,frisk_neg1.5,frisk_neg2,frisk_neg2.5,frisk_neg3,frisk_neg3.5,frisk_neg4,frisk_neg4.5,frisk_neg5,frisk_neg5.5,frisk_neg6,frisk_neg6.5,frisk_neg7,frisk_neg7.5,frisk_neg8,frisk_neg8.5,frisk_neg9,frisk_neg9.5,frisk_neg10, 
                        frisk_pos0, frisk_pos0.5,frisk_pos1,frisk_pos1.5,frisk_pos2,frisk_pos2.5,frisk_pos3,frisk_pos3.5,frisk_pos4,frisk_pos4.5,frisk_pos5,frisk_pos5.5,frisk_pos6,frisk_pos6.5,frisk_pos7,frisk_pos7.5,frisk_pos8,frisk_pos8.5,frisk_pos9,frisk_pos9.5,frisk_pos10,
                        fnrisk_0,fnrisk_0.5,fnrisk_1,fnrisk_1.5,fnrisk_2,fnrisk_2.5,fnrisk_3,fnrisk_3.5,fnrisk_4,fnrisk_4.5,fnrisk_5,fnrisk_5.5,fnrisk_6,fnrisk_6.5,fnrisk_7,fnrisk_7.5,fnrisk_8,fnrisk_8.5,fnrisk_9,fnrisk_9.5,fnrisk_10,
                        mrisk_neg0, mrisk_neg0.5,mrisk_neg1,mrisk_neg1.5,mrisk_neg2,mrisk_neg2.5,mrisk_neg3,mrisk_neg3.5,mrisk_neg4,mrisk_neg4.5,mrisk_neg5,mrisk_neg5.5,mrisk_neg6,mrisk_neg6.5,mrisk_neg7,mrisk_neg7.5,mrisk_neg8,mrisk_neg8.5,mrisk_neg9,mrisk_neg9.5,mrisk_neg10,
                        mrisk_pos0, mrisk_pos0.5,mrisk_pos1,mrisk_pos1.5,mrisk_pos2,mrisk_pos2.5,mrisk_pos3,mrisk_pos3.5,mrisk_pos4,mrisk_pos4.5,mrisk_pos5,mrisk_pos5.5,mrisk_pos6,mrisk_pos6.5,mrisk_pos7,mrisk_pos7.5,mrisk_pos8,mrisk_pos8.5,mrisk_pos9,mrisk_pos9.5,mrisk_pos10,
                        mnrisk_0, mnrisk_0.5,mnrisk_1,mnrisk_1.5,mnrisk_2,mnrisk_2.5,mnrisk_3,mnrisk_3.5,mnrisk_4,mnrisk_4.5,mnrisk_5,mnrisk_5.5,mnrisk_6,mnrisk_6.5,mnrisk_7,mnrisk_7.5,mnrisk_8,mnrisk_8.5,mnrisk_9,mnrisk_9.5,mnrisk_10,
                        drisk_neg0, drisk_neg0.5,drisk_neg1,drisk_neg1.5,drisk_neg2,drisk_neg2.5,drisk_neg3,drisk_neg3.5,drisk_neg4,drisk_neg4.5,drisk_neg5,drisk_neg5.5,drisk_neg6,drisk_neg6.5,drisk_neg7,drisk_neg7.5,drisk_neg8,drisk_neg8.5,drisk_neg9,drisk_neg9.5,drisk_neg10,
                        drisk_pos0, drisk_pos0.5,drisk_pos1,drisk_pos1.5,drisk_pos2,drisk_pos2.5,drisk_pos3,drisk_pos3.5,drisk_pos4,drisk_pos4.5,drisk_pos5,drisk_pos5.5,drisk_pos6,drisk_pos6.5,drisk_pos7,drisk_pos7.5,drisk_pos8,drisk_pos8.5,drisk_pos9,drisk_pos9.5,drisk_pos10,
                        dnrisk_0, dnrisk_0.5,dnrisk_1,dnrisk_1.5,dnrisk_2,dnrisk_2.5,dnrisk_3,dnrisk_3.5,dnrisk_4,dnrisk_4.5,dnrisk_5,dnrisk_5.5,dnrisk_6,dnrisk_6.5,dnrisk_7,dnrisk_7.5,dnrisk_8,dnrisk_8.5,dnrisk_9,dnrisk_9.5,dnrisk_10)
  
  simsum <- rbind(simsum, saveres)
}

save(simsum,file="~/Desktop/HPVCPC_Exp.RData")