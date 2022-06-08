#multistate models under mixed-cases interval censoring (random intervals)
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

v1<- ifelse(n_visits>=1, intervals[1:a1],0) #for people with 1+ visits, interval btwn initial visit and 1st follow up 
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

#create data set of observed
id <- seq(num)
obsdat <- rbind(cbind(id,0,r0v2),
                cbind(id,v1,r1v2)[r0v2<3 & v1>0,],
                cbind(id,v2,r2v2)[r1v2<3 & v2>v1,],
                cbind(id,v3,r3v2)[r2v2<3 & v3>v2,],
                cbind(id,v4,r4v2)[r3v2<3 & v4>v3,],
                cbind(id,v5,r5v2)[r4v2<3 & v5>v4,],
                cbind(id,v6,r6v2)[r5v2<3 & v6>v5,],
                cbind(id,v7,r7v2)[r6v2<3 & v7>v6,],
                cbind(id,v8,r8v2)[r7v2<3 & v8>v7,],
                cbind(id,v9,r9v2)[r8v2<3 & v9>v8,],
                cbind(id,v10,r10v2)[r9v2<3 & v10>v9,])

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

#prevalent precancer risk if new HPV infection detected at 5 year visit
p <- pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3]/(pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,2]+pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3])
res <- c(hpvcc.msm$estimates.t,
         hpvcc.msm$QmatricesSE$baseline[1,2],hpvcc.msm$QmatricesSE$baseline[2,1],hpvcc.msm$QmatricesSE$baseline[2,3],
         hpvcc.msm$QmatricesL$baseline[1,2]<=0.055 & hpvcc.msm$QmatricesU$baseline[1,2]>=0.055,
         hpvcc.msm$QmatricesL$baseline[2,1]<=(1/1.5) & hpvcc.msm$QmatricesU$baseline[2,1]>=(1/1.5),
         hpvcc.msm$QmatricesL$baseline[2,3]<=(1/60) & hpvcc.msm$QmatricesU$baseline[2,3]>=(1/60),
         #risk following HPV-negative at first visit
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
         #risk following HPV-positive at first visit
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
         #risk following newly detected HPV-positive at a visit in 5 years
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
