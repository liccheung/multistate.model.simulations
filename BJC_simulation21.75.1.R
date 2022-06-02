#misspecified Markov model (truth is precursor to clinical cancer having Gamma form)
#normally distributed covariate
#n.sim - number of simulation data sets to run
#output1 - filename of output
#seed - seed number for simulation run
sim.fun2 <- function(n.sim,output1,seed){
#no measurement error
#early detection has no effect on subsequent clinical cancer/death
nsim <- n.sim #number of simulation datasets
n <- 10000 #number of samples in each dataset
num <- nsim*n  
set.seed(seed)

#set parameters
#(a1,b1) are gamma parameters for time to cancer precursor
b11 <- -5.5
b12 <- 0.5
a1 <- 1
#time from cancer precursor to clinical cancer
b21 <- -3+log(2)
b22 <- 0.5
a2 <- 2   
#time to death from healthy
b31 <- -4
b32 <- 0.2
a3 <- 1   
#time to death from cancer precursor
b41 <- -4
b42 <- 0.2
a4 <- 1 
#40% survival after clinical detection with lung cancer
b51 <- 0
b52 <- 0
a5 <- 1
theta <- c(b11,b31,b21,b41,b51,b12,b32,b22,b42,b52)
#Simulate times
x <- rnorm(num,0,1)
t1 <- rgamma(num,shape=a1,rate=exp(b11)*exp(b12*x)) #time from healthy (age 45) to cancer precursor
t2 <- rgamma(num,shape=a2,rate=exp(b21)*exp(b22*x)) #sojourn time
t3 <- rgamma(num,shape=a3,rate=exp(b31)*exp(b32*x)) #time from healthy to death
t4 <- rgamma(num,shape=a4,rate=exp(b41)*exp(b42*x)) #time from cancer precursor to death
t5 <- rgamma(num,shape=a5,rate=exp(b51)*exp(b52*x)) #time from clinical cancer to death
#set to ~25% overall deaths after 15 years according to LYFS-CT mortality model for smokers, age 40-84
t_death <- ifelse(t3<=t1,t3,
                  ifelse(((t1+t4)<=(t1+t2)),(t1+t4),(t1+t2+t5)))
#14% with cancer before death
t_can <- ifelse((t1<t3) & ((t1+t2)<(t1+t4)),(t1+t2),NA)
t_precan <- ifelse(t1<t3,t1,NA)
v0 <- 0 #runif(num,min=0,max=25) #enrollment date and first screen for screening arm
n_visits <- sample.int(15, num, replace= TRUE)-1  #15 #uncomment for no dropout
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
a11 <- sum(n_visits>=11)
a12 <- sum(n_visits>=12)
a13 <- sum(n_visits>=13)
a14 <- sum(n_visits>=14)
a15 <- sum(n_visits>=15)

v1<- ifelse(n_visits>=1, intervals[1:a1],NA) #for people with 1+ visits, interval btwn initial visit and 1st follow up
v2 <- v1 + ifelse(n_visits>=2, intervals[(a1+1):(a1+a2)],NA) #interval between visits for people with 1-2 follow ups
v3 <- v2 + ifelse(n_visits>=3, intervals[(a1+a2+1):(a1+a2+a3)],NA)
v4 <- v3 + ifelse(n_visits>=4, intervals[(a1+a2+a3+1):(a1+a2+a3+a4)],NA)
v5 <- v4 + ifelse(n_visits>=5, intervals[(a1+a2+a3+a4+1):(a1+a2+a3+a4+a5)],NA)
v6 <- v5 + ifelse(n_visits>=6, intervals[(a1+a2+a3+a4+a5+1):(a1+a2+a3+a4+a5+a6)],NA)
v7 <- v6 + ifelse(n_visits>=7, intervals[(a1+a2+a3+a4+a5+a6+1):(a1+a2+a3+a4+a5+a6+a7)],NA)
v8 <- v7 + ifelse(n_visits>=8, intervals[(a1+a2+a3+a4+a5+a6+a7+1):(a1+a2+a3+a4+a5+a6+a7+a8)],NA)
v9 <- v8 + ifelse(n_visits>=9, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9)],NA)
v10 <- v9 + ifelse(n_visits>=10, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)],NA)
v11 <- v10 + ifelse(n_visits>=11, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11)],NA)
v12 <- v11 + ifelse(n_visits>=12, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)],NA)
v13 <- v12 + ifelse(n_visits>=13, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13)],NA)
v14 <- v13 + ifelse(n_visits>=14, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14)],NA)
v15 <- v14 + ifelse(n_visits>=15, intervals[(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+1):(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15)],NA)

cens <- v0+15 #study administratively censored after 15 years
vr <- function(vt){
  return(ifelse(cens<=vt | t_death<=vt | (is.na(t_can)==0 & t_can<=vt),NA,
                ifelse(is.na(t_precan)==0 & t_precan<=vt,2,1)))
}
r1 <- vr(v1)
r2 <- vr(v2)
r3 <- vr(v3)
r4 <- vr(v4)
r5 <- vr(v5)
r6 <- vr(v6)
r7 <- vr(v7)
r8 <- vr(v8)
r9 <- vr(v9)
r10 <- vr(v10)
r11 <- vr(v11)
r12 <- vr(v12)
r13 <- vr(v13)
r14 <- vr(v14)
r15 <- vr(v15)

#observed data - person id, time, x, state, obstype
id <- seq(num)

obsdat <- rbind(cbind(id,0,x,1,1),
                cbind(id,v1,x,r1,1)[is.na(r1)==0,],
                cbind(id,v2,x,r2,1)[is.na(r2)==0,],
                cbind(id,v3,x,r3,1)[is.na(r3)==0,],
                cbind(id,v4,x,r4,1)[is.na(r4)==0,],
                cbind(id,v5,x,r5,1)[is.na(r5)==0,],
                cbind(id,v6,x,r6,1)[is.na(r6)==0,],
                cbind(id,v7,x,r7,1)[is.na(r7)==0,],
                cbind(id,v8,x,r8,1)[is.na(r8)==0,],
                cbind(id,v9,x,r9,1)[is.na(r9)==0,],
                cbind(id,v10,x,r10,1)[is.na(r10)==0,],
                cbind(id,v11,x,r11,1)[is.na(r11)==0,],
                cbind(id,v12,x,r12,1)[is.na(r12)==0,],
                cbind(id,v13,x,r13,1)[is.na(r13)==0,],
                cbind(id,v14,x,r14,1)[is.na(r14)==0,],
                cbind(id,v15,x,r15,1)[is.na(r15)==0,],
                cbind(id,t_can,x,3,2)[is.na(t_can)==0 & t_can<=cens,],
                cbind(id,t_death,x,4,3)[t_death<=cens,],
                cbind(id,cens,x,ifelse(is.na(t_can)==0 & t_can<=cens,3,
                                       ifelse(pmax(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,na.rm=TRUE)==2,2,1.5)),3)[t_death>cens,])

obsdat <- data.frame(obsdat)
colnames(obsdat) <- c("id","years","x","state","obstype")
obsdat <- obsdat[order(obsdat[,1],obsdat[,2]),]

theta <- c(b11,b31,b21,b51,b12,b32)
loglik <- function(parest,xsam,loweri,upperi){
  liki <- function(z){
    perdat <- xsam[xsam$id==z,]
    #6 parameters - b11, b31=b41, b21, b51, b12=b22, b32=b42, b52=0
    q12 <- exp(parest[1]+parest[5]*perdat$x[1])
    q14 <- exp(parest[2]+parest[6]*perdat$x[1])
    q23 <- exp(parest[3]+parest[5]*perdat$x[1])
    q24 <- exp(parest[2]+parest[6]*perdat$x[1])    
    q34 <- exp(parest[4])
    
    lstate1 <- perdat$years[max(which(perdat$state==1))]
    fstate2 <- ifelse(sum(perdat$state==2)>0,perdat$years[min(which(perdat$state==2))],NA)
    lstate2 <- ifelse(sum(perdat$state==2)>0,perdat$years[max(which(perdat$state==2))],NA)
    fstate3 <- ifelse(sum(perdat$state==3)>0,perdat$years[min(which(perdat$state==3))],NA)
    lstate3 <- ifelse(sum(perdat$state==3)>0,perdat$years[max(which(perdat$state==3))],NA)
    fstate4 <- ifelse(sum(perdat$state==4)>0,perdat$years[min(which(perdat$state==4))],NA)
    fstate1.5 <- ifelse(sum(perdat$state==1.5)>0,perdat$years[min(which(perdat$state==1.5))],NA)   
    p11 <- exp(-1*(q12+q14)*lstate1)
    p12 <- ifelse(is.na(fstate2)==0,q12*(exp(-1*(q23+q24)*(fstate2-lstate1))-exp(-1*(q12+q14)*(fstate2-lstate1)))/(q12+q14-q23-q24),1)
    #either a true 1->1 followed by instant 4 transition or a 1->2 followed by instant 4 transition
    p14 <- ifelse(is.na(fstate4)==0 & is.na(fstate2)==1 & is.na(fstate3)==1,exp(-1*(q12+q14)*(fstate4-lstate1))*q14 +
             q12*(exp(-1*(q23+q24)*(fstate4-lstate1))-exp(-1*(q12+q14)*(fstate4-lstate1)))/(q12+q14-q23-q24)*q24,1)
    p22 <- ifelse(is.na(fstate2)==0,exp(-1*(q23+q24)*(lstate2-fstate2)),1)
    #clinical cancer is an exact time: 2->2 followed by instant 3 transition
    p23 <- ifelse(is.na(fstate3)==0 & is.na(fstate2)==0,exp(-1*(q23+q24)*(fstate3-lstate2))*q23,1)
    #clinical cancer is an exact time: 1->2 transition followed by instant 2->3 transition
    p13 <- ifelse(is.na(fstate3)==0 & is.na(fstate2)==1,q12*(exp(-1*(q23+q24)*(fstate3-lstate1))-exp(-1*(q12+q14)*(fstate3-lstate1)))/(q12+q14-q23-q24)*q23,1)
    #death from preclinical cancer is an exact time: 2->2 followed by instant 4 transition
    p24 <- ifelse(is.na(fstate2)==0 & is.na(fstate3)==1 & is.na(fstate4)==0,exp(-1*(q23+q24)*(fstate4-lstate2))*q24,1)
    p33 <- ifelse(is.na(fstate3)==0 & is.na(fstate4)==1,exp(-1*q34*(lstate3-fstate3)),1)
    #death from clinical cancer is an exact time: 3->3 followed by instant 4 transition
    p34 <- ifelse(is.na(fstate3)==0 & is.na(fstate4)==0,exp(-1*q34*(fstate4-lstate3))*q34,1)
    p1.5 <- ifelse(is.na(fstate1.5)==0,exp(-1*(q12+q14)*(fstate1.5-lstate1)) +
                    q12*(exp(-1*(q23+q24)*(fstate1.5-lstate1))-exp(-1*(q12+q14)*(fstate1.5-lstate1)))/(q12+q14-q23-q24),1)
    lik <- p11*p12*p13*p14*p22*p23*p24*p33*p34*p1.5
    return(lik)                                          
  }  
  lli <- sapply(seq(loweri,upperi),liki)
  print(c(-1*sum(log(lli)),parest))
  return(-1*sum(log(lli)))
}

simsum <- data.frame()
for (j in 1:nsim){
  print(j)
  loweri <- (j-1)*n+1
  upperi <- j*n
  lc <- data.frame(obsdat[which(obsdat$id %in% seq(loweri,upperi)),])
  lik_sim <-  optim(par = theta, fn = loglik, xsam = lc, loweri=loweri, upperi=upperi,
                    method="BFGS", control=list(trace=TRUE), hessian=TRUE)
  lcl <- lik_sim$par - 1.96*sqrt(diag(solve(lik_sim$hessian)))
  ucl <- lik_sim$par + 1.96*sqrt(diag(solve(lik_sim$hessian)))
  #confidence intervals for sojourn times in healthy, cancer precursor, and clinical cancer states
  ST1_lcl <- 1/sum(exp(lik_sim$par[1:2])) - 1.96*1/sum(exp(lik_sim$par[1:2]))^2*
    sqrt(exp(2*lik_sim$par[1])*solve(lik_sim$hessian)[1,1]+exp(2*lik_sim$par[2])*solve(lik_sim$hessian)[2,2]+2*exp(sum(lik_sim$par[1:2]))*solve(lik_sim$hessian)[1,2])
  ST1_ucl <- 1/sum(exp(lik_sim$par[1:2])) + 1.96*1/sum(exp(lik_sim$par[1:2]))^2*
    sqrt(exp(2*lik_sim$par[1])*solve(lik_sim$hessian)[1,1]+exp(2*lik_sim$par[2])*solve(lik_sim$hessian)[2,2]+2*exp(sum(lik_sim$par[1:2]))*solve(lik_sim$hessian)[1,2])
  ST2_lcl <- 1/sum(exp(lik_sim$par[2:3])) - 1.96*1/sum(exp(lik_sim$par[2:3]))^2*
    sqrt(exp(2*lik_sim$par[3])*solve(lik_sim$hessian)[3,3]+exp(2*lik_sim$par[2])*solve(lik_sim$hessian)[2,2]+2*exp(sum(lik_sim$par[2:3]))*solve(lik_sim$hessian)[2,3])
  ST2_ucl <- 1/sum(exp(lik_sim$par[2:3])) + 1.96*1/sum(exp(lik_sim$par[2:3]))^2*
    sqrt(exp(2*lik_sim$par[3])*solve(lik_sim$hessian)[3,3]+exp(2*lik_sim$par[2])*solve(lik_sim$hessian)[2,2]+2*exp(sum(lik_sim$par[2:3]))*solve(lik_sim$hessian)[2,3])
  ST3_lcl <- exp(-1*lik_sim$par[4]) - 1.96*exp(-1*lik_sim$par[4])*sqrt(solve(lik_sim$hessian)[4,4])
  ST3_ucl <- exp(-1*lik_sim$par[4]) + 1.96*exp(-1*lik_sim$par[4])*sqrt(solve(lik_sim$hessian)[4,4])
  ressum <- c(lik_sim$par,1/sum(exp(lik_sim$par[1:2])),1/sum(exp(lik_sim$par[2:3])),exp(-1*lik_sim$par[4]),
              diag(sqrt(solve(lik_sim$hessian))),
              1/sum(exp(lik_sim$par[1:2]))^2*
              sqrt(exp(2*lik_sim$par[1])*solve(lik_sim$hessian)[1,1]+exp(2*lik_sim$par[2])*solve(lik_sim$hessian)[2,2]+2*exp(sum(lik_sim$par[1:2]))*solve(lik_sim$hessian)[1,2]),
              1/sum(exp(lik_sim$par[2:3]))^2*
              sqrt(exp(2*lik_sim$par[3])*solve(lik_sim$hessian)[3,3]+exp(2*lik_sim$par[2])*solve(lik_sim$hessian)[2,2]+2*exp(sum(lik_sim$par[2:3]))*solve(lik_sim$hessian)[2,3]),
              exp(-1*lik_sim$par[4])*sqrt(solve(lik_sim$hessian)[4,4]),
              lik_sim$par[1:3]-theta[1:3],
              exp(lik_sim$par[4])-exp(theta[4]),
              lik_sim$par[5:6]-theta[5:6],
              lcl[1:3]<=theta[1:3] & ucl[1:3]>=theta[1:3],
              (exp(lik_sim$par[4])-1.96*exp(lik_sim$par[4])*sqrt(solve(lik_sim$hessian)[4,4])) <= exp(theta[4]) &
                 (exp(lik_sim$par[4])+1.96*exp(lik_sim$par[4])*sqrt(solve(lik_sim$hessian)[4,4])) >= exp(theta[4]),
              lcl[5:6]<=theta[5:6] & ucl[5:6]>=theta[5:6],              
              1/sum(exp(lik_sim$par[1:2]))-1/sum(exp(theta[1:2])),
              ST1_lcl <= 1/sum(exp(theta[1:2])) & ST1_ucl >= 1/sum(exp(theta[1:2])),
              1/sum(exp(lik_sim$par[2:3]))-(2*exp(theta[3]) + exp(theta[2]))/(exp(theta[3]) + exp(theta[2]))^2,
              ST2_lcl <= (2*exp(theta[3]) + exp(theta[2]))/(exp(theta[3]) + exp(theta[2]))^2 & 
                ST2_ucl >= (2*exp(theta[3]) + exp(theta[2]))/(exp(theta[3]) + exp(theta[2]))^2,
              exp(-1*lik_sim$par[4]) - exp(-1*theta[4]),
              ST3_lcl <= exp(-1*theta[4]) & ST3_ucl >= exp(-1*theta[4]))
  simsum <- rbind(simsum,ressum)
}
#colnames(simsum) <- c("A.1.2 Bias","A.1.4 Bias","A.2.3 Bias","exp(A.3.4 Bias)","A.1.2 X Bias","A.1.4 X Bias",
#                      "A.1.2 Cov","A.1.4 Cov","A.2.3 Cov","exp(A.3.4) Cov","A.1.2 X Cov","A.1.4 X Cov",
#                      "A.1.ST Bias","A.1.ST Cov","A.2.ST Bias","A.2.ST Cov","A.3.ST Bias","A.3.ST Cov")

saveRDS(simsum, file=output1)
}

