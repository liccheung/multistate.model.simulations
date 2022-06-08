#This program generates data under patient advoidant observation schemes.
#Precancer risk estimates are from healthy vs. precancer survival approach (EM-ICM to account for left-, interval-, and right-censoring)
#Precancer risk estimates are for initially HPV-negative, initially HPV-positive, and new HPV-positive

#increases number of infections to 3
#only 1 type
library(tidyverse)
library(survival)
library(Icens)

#General simulation settings
nsim <- 500 #number of simulation datasets 
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

#dependent visits - 3 year interval following HPV-negative and 5-year following HPV-positive results
maxtime <- ifelse(rbinom(num,1,.3),12,runif(num,0,12))  #increased maxtime to 12
r0v3 <- h1
v1_intervals <- ifelse(h1==0,rnorm(sum(h1==0),3,0.5),rnorm(sum(h1>0),5,0.5))
v1_intervals[v1_intervals<0] <- 0
v1v3 <- ifelse(v1_intervals<=maxtime,v1_intervals,0)
nomore <- rep(0,num)
nomore <- ifelse(v1_intervals>maxtime | nomore==1,1,0)
r1v3 <- time_point(v1v3)

v2_intervals <- ifelse(r1v3==0,rnorm(sum(r1v3==0),3,0.5),rnorm(sum(r1v3>0),5,0.5))
v2_intervals[v2_intervals<0] <- 0
v2v3 <- ifelse((v1v3 + v2_intervals)<=maxtime,(v1v3 + v2_intervals),v1v3)
nomore <- ifelse((v1v3 + v2_intervals)>maxtime | nomore==1,1,0)
r2v3 <- time_point(v2v3)

v3_intervals <- ifelse(r2v3==0,rnorm(sum(r2v3==0),3,0.5),rnorm(sum(r2v3>0),5,0.5))
v3_intervals[v3_intervals<0] <- 0
v3v3 <- ifelse((v2v3 + v3_intervals)<=maxtime,(v2v3 + v3_intervals),v2v3)
nomore <- ifelse((v2v3 + v3_intervals)>maxtime | nomore==1,1,0)
r3v3 <- time_point(v3v3)

v4_intervals <- ifelse(r3v3==0,rnorm(sum(r3v3==0),3,0.5),rnorm(sum(r3v3>0),5,0.5))
v4_intervals[v4_intervals<0] <- 0
v4v3 <- ifelse((v3v3 + v4_intervals)<=maxtime,(v3v3 + v4_intervals),v3v3)
nomore <- ifelse((v3v3 + v4_intervals)>maxtime | nomore==1,1,0)
r4v3 <- time_point(v4v3)

v5_intervals <- ifelse(r4v3==0,rnorm(sum(r4v3==0),3,0.5),rnorm(sum(r4v3>0),5,0.5))
v5_intervals[v5_intervals<0] <- 0
v5v3 <- ifelse((v4v3 + v5_intervals)<=maxtime,(v4v3 + v5_intervals),v4v3)
nomore <- ifelse((v4v3 + v5_intervals)>maxtime | nomore==1,1,0)
r5v3 <- time_point(v5v3)

v6_intervals <- ifelse(r5v3==0,rnorm(sum(r5v3==0),3,0.5),rnorm(sum(r5v3>0),5,0.5))
v6_intervals[v6_intervals<0] <- 0
v6v3 <- ifelse((v5v3 + v6_intervals)<=maxtime,(v5v3 + v6_intervals),v5v3)
nomore <- ifelse((v5v3 + v6_intervals)>maxtime | nomore==1,1,0)
r6v3 <- time_point(v6v3)

v7_intervals <- ifelse(r6v3==0,rnorm(sum(r6v3==0),3,0.5),rnorm(sum(r6v3>0),5,0.5))
v7_intervals[v7_intervals<0] <- 0
v7v3 <- ifelse((v6v3 + v7_intervals)<=maxtime,(v6v3 + v7_intervals),v6v3)
nomore <- ifelse((v6v3 + v7_intervals)>maxtime | nomore==1,1,0)
r7v3 <- time_point(v7v3)

v8_intervals <- ifelse(r7v3==0,rnorm(sum(r7v3==0),3,0.5),rnorm(sum(r7v3>0),5,0.5))
v8_intervals[v8_intervals<0] <- 0
v8v3 <- ifelse((v7v3 + v8_intervals)<=maxtime,(v7v3 + v8_intervals),v7v3)
nomore <- ifelse((v7v3 + v8_intervals)>maxtime | nomore==1,1,0)
r8v3 <- time_point(v8v3)

v9_intervals <- ifelse(r8v3==0,rnorm(sum(r8v3==0),3,0.5),rnorm(sum(r8v3>0),5,0.5))
v9_intervals[v9_intervals<0] <- 0
v9v3 <- ifelse((v8v3 + v9_intervals)<=maxtime,(v8v3 + v9_intervals),v8v3)
nomore <- ifelse((v8v3 + v9_intervals)>maxtime | nomore==1,1,0)
r9v3 <- time_point(v9v3)

v10_intervals <- ifelse(r9v3==0,rnorm(sum(r9v3==0),3,0.5),rnorm(sum(r9v3>0),5,0.5))
v10_intervals[v10_intervals<0] <- 0
v10v3 <- ifelse((v9v3 + v10_intervals)<=maxtime,(v9v3 + v10_intervals),v9v3)
nomore <- ifelse((v9v3 + v10_intervals)>maxtime | nomore==1,1,0)
r10v3 <- time_point(v10v3)

#create data set of observed
obsdat3 <- data.frame(r0v3, r1v3, r2v3, r3v3, r4v3, r5v3, r6v3, r7v3, r8v3, r9v3, r10v3, v1v3, v2v3, v3v3, v4v3, v5v3, v6v3, v7v3, v8v3, v9v3, v10v3)

#estimate for baseline HPV+, baseline HPV-, and new HPV+ (where the previous recorded visit was HPV-)
simsum <- data.frame()
for (j in 1:nsim){
  print(j)
  lower <- (j-1)*n+1
  upper <- j*n
  #dataset j  
  xsam3 <- data.frame(obsdat3[lower:upper,])
  
  #Relative to first visit, time interval in which precancer occurred
  #precancers present at time of first visit occur in time interval (-0.01,0]
  #time to last precancer negative assessment
  lower3 <- ifelse(xsam3$r10v3<3,xsam3$v10v3,
                   ifelse(xsam3$r9v3<3,xsam3$v9v3, 
                          ifelse(xsam3$r8v3<3,xsam3$v8v3,
                                 ifelse(xsam3$r7v3<3,xsam3$v7v3,
                                        ifelse(xsam3$r6v3<3,xsam3$v6v3,
                                               ifelse(xsam3$r5v3<3,xsam3$v5v3,
                                                      ifelse(xsam3$r4v3<3,xsam3$v4v3,
                                                             ifelse(xsam3$r3v3<3,xsam3$v3v3,
                                                                    ifelse(xsam3$r2v3<3,xsam3$v2v3,
                                                                           ifelse(xsam3$r1v3<3,xsam3$v1v3,
                                                                                  ifelse(xsam3$r0v3<3,0,-0.01)))))))))))
  
  #time to precancer diagnosis
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
  drisk_pos0 <- drisk_pos$sigma[max(which(drisk_pos$intmap[2,]<=0))]
  
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
  
  #Time of Detecting Incident HPV
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
  
  #Relative to time of detecting incident HPV, time interval in which precancer occurred
  #precancers present at time of first HPV detection occur in time interval (-0.01,0]
  #time to last precancer negative assessment
  dnlower <-  ifelse(is.na(dnewhpv_t)==1 | dnewhpv_t==0,NA,  #these were never HPV+ or were HPV+ at time 0
                     ifelse(xsam3$r10v3 < 3 & dnewhpv_t <= xsam3$v10v3, xsam3$v10v3-dnewhpv_t,
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

  #time to precancer detection
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
  
  #fit non-parametric estimator using EM-ICM
  H <- cbind(dnlower[is.na(dnlower)==0], dnupper[is.na(dnupper)==0])
  dnfit <- EMICM(H)
  
  #risk for new hpv+ under follow-up avoidance
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
  
  saveres <- data.frame(drisk_neg0, drisk_neg0.5,drisk_neg1,drisk_neg1.5,drisk_neg2,drisk_neg2.5,drisk_neg3,drisk_neg3.5,drisk_neg4,drisk_neg4.5,drisk_neg5,drisk_neg5.5,drisk_neg6,drisk_neg6.5,drisk_neg7,drisk_neg7.5,drisk_neg8,drisk_neg8.5,drisk_neg9,drisk_neg9.5,drisk_neg10,
                        drisk_pos0, drisk_pos0.5,drisk_pos1,drisk_pos1.5,drisk_pos2,drisk_pos2.5,drisk_pos3,drisk_pos3.5,drisk_pos4,drisk_pos4.5,drisk_pos5,drisk_pos5.5,drisk_pos6,drisk_pos6.5,drisk_pos7,drisk_pos7.5,drisk_pos8,drisk_pos8.5,drisk_pos9,drisk_pos9.5,drisk_pos10,
                        dnrisk_0, dnrisk_0.5,dnrisk_1,dnrisk_1.5,dnrisk_2,dnrisk_2.5,dnrisk_3,dnrisk_3.5,dnrisk_4,dnrisk_4.5,dnrisk_5,dnrisk_5.5,dnrisk_6,dnrisk_6.5,dnrisk_7,dnrisk_7.5,dnrisk_8,dnrisk_8.5,dnrisk_9,dnrisk_9.5,dnrisk_10)
  
  simsum <- rbind(simsum, saveres)
}

simsum[simsum==1] <- NA
simmean <- colMeans(simsum, na.rm=TRUE)
save(simsum,file="~/Desktop/HPVCPC_Exp3.RData")
