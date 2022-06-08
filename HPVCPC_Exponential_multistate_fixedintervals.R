#multistate models under fixed intervals
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
r8 <- time_point(rep(8,num)) #for checking
r8.5 <- time_point(rep(8.5,num))
r9 <- time_point(rep(9,num))
r9.5 <- time_point(rep(9.5,num))
r10 <- time_point(rep(10,num))

#create data set of observed
id <- seq(num)
obsdat <- rbind(cbind(id, 0, r0),
                cbind(id,0.5,r0.5)[r0<3,],
                cbind(id,1,r1)[r0.5<3,],
                cbind(id,1.5,r1.5)[r1<3,],
                cbind(id,2,r2)[r1.5<3,],
                cbind(id,2.5,r2.5)[r2<3,],
                cbind(id,3,r3)[r2.5<3,],
                cbind(id,3.5,r3.5)[r3<3,],
                cbind(id,4,r4)[r3.5<3,],
                cbind(id,4.5,r4.5)[r4<3,],
                cbind(id,5,r5)[r4.5<3,],
                cbind(id,5.5,r5.5)[r5<3,],
                cbind(id,6,r6)[r5.5<3,],
                cbind(id,6.5,r6.5)[r6<3,],
                cbind(id,7,r7)[r6.5<3,],
                cbind(id,7.5,r7.5)[r7<3,],
                cbind(id,8,r8)[r7.5<3,],
                cbind(id,8.5,r8.5)[r8<3,],
                cbind(id,9,r9)[r8.5<3,],
                cbind(id,9.5,r9.5)[r9<3,],
                cbind(id,10,r10)[r9.5<3,])

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

#prevalent precancer if new HPV infection detected at year 5
p <- pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3]/(pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,2]+pmatrix.msm(hpvcc.msm,t1=0,t=5)[1,3])
res <- c(hpvcc.msm$estimates.t,
         hpvcc.msm$QmatricesSE$baseline[1,2],hpvcc.msm$QmatricesSE$baseline[2,1],hpvcc.msm$QmatricesSE$baseline[2,3],
         hpvcc.msm$QmatricesL$baseline[1,2]<=0.055 & hpvcc.msm$QmatricesU$baseline[1,2]>=0.055,
         hpvcc.msm$QmatricesL$baseline[2,1]<=(1/1.5) & hpvcc.msm$QmatricesU$baseline[2,1]>=(1/1.5),
         hpvcc.msm$QmatricesL$baseline[2,3]<=(1/60) & hpvcc.msm$QmatricesU$baseline[2,3]>=(1/60),
         #precancer risk if HPV-negative at first visit
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
         #precancer risk if HPV-positive at first visit
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
         #precancer risk if newly HPV-positive at a visit in 5 years
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
