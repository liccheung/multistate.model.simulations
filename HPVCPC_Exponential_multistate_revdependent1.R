#multistate models under follow-up advoidance
sim.fun <- function(n.sim,output1,seed){
library(msm)
library(survival)
#increases number of infections to 3
#only 1 type
#General simulation settings
nsim <- n.sim #number of simulation datasets 
n <- 10000 #number of samples in each dataset
num <- nsim*n
set.seed(seed)

#LCC: Parameter values - set up so time represents years
p1 <- .1
p2 <- .5
p3 <- .055 #LCC: Didem's paper - 5.5% of HPV+ have prevalent CIN3+
#LCC: Didem's paper -4% of HPV- have HPV ~3 years later
lambda1 <- 0.055
#LCC: clearance parameters from Sally's paper
shape21 <- 1 #0.702
scale21 <- 1.5
#LCC: Didem's paper - approximately 2% of HPV+ are already cervical precancer at first detection
#LCC: Didem's paper - 3.8% progressed to cervical precancer in 5 years
shape22 <- 1
scale22 <- 60
#mean duration is scale22*gamma(1+1/shape22) is approximately 9.6 years from HPV acquisition to clearance

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

#dependent visits - 5 years for HPV+ result in very few visits before year 10, 
#try 3 years and add more variance - may need smaller
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
id <- seq(num)
obsdat <- rbind(cbind(id, 0, r0v3),
                cbind(id,v1v3,r1v3)[v1v3<=maxtime & r0v3<3,],
                cbind(id,v2v3,r2v3)[v2v3<=maxtime & r1v3<3,],
                cbind(id,v3v3,r3v3)[v3v3<=maxtime & r2v3<3,],
                cbind(id,v4v3,r4v3)[v4v3<=maxtime & r3v3<3,],
                cbind(id,v5v3,r5v3)[v5v3<=maxtime & r4v3<3,],
                cbind(id,v6v3,r6v3)[v6v3<=maxtime & r5v3<3,],
                cbind(id,v7v3,r7v3)[v7v3<=maxtime & r6v3<3,],
                cbind(id,v8v3,r8v3)[v8v3<=maxtime & r7v3<3,],
                cbind(id,v9v3,r9v3)[v9v3<=maxtime & r8v3<3,],
                cbind(id,v10v3,r10v3)[v10v3<=maxtime & r9v3<3,])

obsdat <- data.frame(obsdat)
colnames(obsdat) <- c("id","years","state")
obsdat <- obsdat[order(obsdat[,1],obsdat[,2]),]
obsdat$state[obsdat$state %in% c(0,1)] <- obsdat$state[obsdat$state %in% c(0,1)] + 1
#specify Markov model
Q <- rbind(c(0,lambda1,0),
           c(1/scale21,0,1/scale22),
           c(0,0,0))
theta <- c(lambda1,1/scale21,1/scale22)

#estimate for baseline HPV+, baseline HPV-, and new HPV+ (where the previous recorded visit was HPV-)
hpvcc.msm <- msm(state ~ years, subject=id, data=obsdat, qmatrix=Q, opt.method="optim")

p <- pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3]/(pmatrix.msm(hpvcc.msm,t1=0,t=10)[1,2]+pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3])
res <- c(hpvcc.msm$estimates.t,
         hpvcc.msm$QmatricesSE$baseline[1,2],hpvcc.msm$QmatricesSE$baseline[2,1],hpvcc.msm$QmatricesSE$baseline[2,3],
         hpvcc.msm$QmatricesL$baseline[1,2]<=0.055 & hpvcc.msm$QmatricesU$baseline[1,2]>=0.055,
         hpvcc.msm$QmatricesL$baseline[2,1]<=(1/1.5) & hpvcc.msm$QmatricesU$baseline[2,1]>=(1/1.5),
         hpvcc.msm$QmatricesL$baseline[2,3]<=(1/60) & hpvcc.msm$QmatricesU$baseline[2,3]>=(1/60),
         pmatrix.msm(hpvcc.msm,t1=0,t=0.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=1)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=1.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=2)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=2.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=3)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=3.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=4)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=4.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=5.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=6)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=6.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=7)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=7.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=8)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=8.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=9)[1,3],
         pmatrix.msm(hpvcc.msm,t1=0,t=9.5)[1,3],pmatrix.msm(hpvcc.msm,t1=0,t=10)[1,3],
         sum(h1==3)/sum(h1>=1),
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=0.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=1)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=1.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=2)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=2.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=3)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=3.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=4)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=4.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=5.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=6)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=6.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=7)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=7.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=8)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=8.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=9)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=9.5)[2,3],
         sum(h1==3)/sum(h1>=1)+(1-sum(h1==3)/sum(h1>=1))*pmatrix.msm(hpvcc.msm,t1=0,t=10)[2,3],
         p,
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=0.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=1)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=1.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=2)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=2.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=3)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=3.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=4)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=4.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=5.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=6)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=6.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=7)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=7.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=8)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=8.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=9)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=9.5)[2,3],
         p+(1-p)*pmatrix.msm(hpvcc.msm,t1=0,t=10)[2,3])

saveRDS(res, file=output1)
}
