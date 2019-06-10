############################################
##   example code for density mediation   ##
##   developed by: Katrina Devick         ##
##   last updated: 9 June 2019            ##
############################################


##########################
####  **** NOTE ****  ####
##########################

## in order for this code to work properly, you need to:
## 1) download the file LDDPdensity.f that is on the GitHub repository with this file (https://github.com/kdevick/densitymediation)
## 2) replace the current LDDPdensity.f file in the DPpackage with this file
## 3) recompile and load this new DPpackage 


#### load required libraries 
library(devtools)
build("~/DPpackage")
install("~/DPpackage")
library(DPpackage, lib.loc = "~/")
library(sfsmisc)
library(utils)
library(survival)
library(MASS)
library(gtools)
library(plyr)
library(BSGW)




#################################################### 
####   format CanCORS data to use in analyses   ####
#################################################### 

load("~/Dropbox/cancors.RData")

## center your continous mediator variable if including a quadratic effect of the mediator on survival
cancors$bmi.c <- scale(cancors$bmi.dx, scale=FALSE)
cancors$bmisq.c <- (cancors$bmi.c)^2

data.comp <- data.frame(na.omit(cbind(cancors$truncsurv, cancors$fiveyrcens, cancors$race,cancors$bmi.dx,cancors$female,cancors$agecat1,cancors$agecat3,cancors$stage2,cancors$stage3,cancors$stage4,cancors$income1,cancors$income3,cancors$income45,cancors$bmi.dx.cat,cancors$bmisq.c,cancors$bmi.c)))
colnames(data.comp) <- c("truncsurv","fiveyrcens","race","bmi.dx","female","agecat1","agecat3","stage2","stage3","stage4","income1","income3","income45","bmi.dx.cat","bmisq.c","bmi.c")



#################################################### 
####         covariate patterns in data         ####
#################################################### 

## create a matrix of complete covariate data (not to include the exposure -- race)
data.complete <- data.frame(na.omit(cbind(cancors$bmi.dx,cancors$female,cancors$agecat1,cancors$agecat3,cancors$stage2,cancors$stage3,cancors$stage4,cancors$income1,cancors$income3,cancors$income45)))
colnames(data.complete) <- c("bmi.dx","female","agecat1","agecat3","stage2","stage3","stage4","income1","income3","income45")
n.complete <- dim(data.complete)[1] ## use this n so that probablities used in calc sum to 1
x<-ddply(data.complete,c("female","agecat1","agecat3","stage2","stage3","stage4","income1","income3","income45"),function(x) count=nrow(x)/n.complete)

## create a matrix of covariate patterns and empirical probabilities 
covpatterns <- as.matrix(x[,1:(dim(x)[2]-1)])
prob.covpatterns <- as.matrix(x[,dim(x)[2]])




#################################################### 
####      fit the Bayesian disparity model      ####
#################################################### 

## disp model with race*age interactions 
iter <- 320000
burnin <- 20000
set.seed(1112017)
outcome.disp.bayes.BEST <- bsgw(Surv(truncsurv +0.0001,fiveyrcens,type="right") ~ race+race*agecat1+race*agecat3+female+agecat1+agecat3+stage2+stage3+stage4+income1+income3+income45,data=data.comp, init="survreg"
                                , ordweib=TRUE, scale=0, control=bsgw.control(scalex=FALSE, iter=iter, burnin=burnin, sd.thresh=1e-4 , nskip=1000))
save(outcome.disp.bayes.BEST,file="outcome_disp_bayes_BEST.RData")


## load model and obtain coefficients
load("outcome_disp_bayes_BEST.RData")

## a vector with the indices of which posterior samples to keep after thin/burn in
burnthin <- seq(20001,320000,by=100)

## format output to get coefficient of interest
sum.disp.BEST <- summary(outcome.disp.bayes.BEST)
disp.scale.BEST <- mean(sum.disp.BEST$survreg.scale$median)
disp.sam.BEST <- outcome.disp.bayes.BEST$smp$beta[burnthin,]
### 
thetas.disp.BEST <- -apply(disp.sam.BEST,2,mean)*disp.scale.BEST
disp.bayes.BEST <- exp(thetas.disp.BEST[2])
quantile(exp(-disp.sam.BEST[,2]*disp.scale.BEST), probs=c(0.025,0.975))
save(thetas.disp.BEST,sum.disp.BEST,file="thetas_disp_bayes_BEST.RData")
load("thetas_disp_bayes_BEST.RData")



#################################################### 
####       fit the Bayesian outcome model       ####
####       (now including BMI/BMI^2 effects)    ####
#################################################### 

## outcome model including race*age interactions
iter <- 320000
burnin <- 20000

set.seed(1112017)
outcome.bayes.bmisq.BEST <-bsgw(Surv(truncsurv +0.0001,fiveyrcens,type="right") ~ race+bmi.c+bmisq.c+race*bmi.c+race*bmisq.c+female+agecat1+agecat3+stage2+stage3+stage4+income1+income3+income45+
                                  race*agecat1+race*agecat3, data=data.comp, init="survreg",
                                ordweib=TRUE, scale=0, control=bsgw.control(scalex=FALSE, iter=iter, burnin=burnin, sd.thresh=1e-4 , nskip=1000))
save(outcome.bayes.bmisq.BEST,file="outcome_bayes_bmisq_BEST.RData")



## load model and obtain coefficients
load("outcome_bayes_bmisq_BEST.RData")

## a vector with the indices of which posterior samples to keep after thin/burn in
burnthin <- seq(20001,320000,by=100)

## format output to get coefficient of interest and save
sum.bmisq.BEST <- summary(outcome.bayes.bmisq.BEST)
bmisq.scale.BEST <-  mean(sum.bmisq.BEST$survreg.scale$median)
outcome.bmisq.sam.BEST <- outcome.bayes.bmisq.BEST$smp$beta[burnthin,]
outcome.bmisq.shapesam.BEST <- outcome.bayes.bmisq.BEST$smp$betas[burnthin,]
thetas.bayes.bmisq.BEST <- -apply(outcome.bmisq.sam.BEST,2,mean)*bmisq.scale.BEST
thetassam.bmisq.BEST <- -outcome.bmisq.sam.BEST*bmisq.scale.BEST
save(thetassam.bmisq.BEST,sum.bmisq.BEST, file="thetassam_bmisq_BEST.RData")




#################################################### 
####    fit density regression for mediator     ####
#################################################### 

mcmc <- list(nburn = 10000, nsave = 3000, nskip = 50, ndisplay = 200) 


nrec <- length(cancors$bmi.c)
W.temp <- cbind(cancors$stage1,cancors$bmi.c,rep(1, nrec), cancors$race, cancors$female, cancors$agecat1, 
                cancors$agecat3, cancors$stage2, cancors$stage3, cancors$stage4, cancors$income1, cancors$income3, cancors$income45)
colnames(W.temp) <- c("stage1","y","int", "x", "female", "agecat1", "agecat3", "stage2", "stage3", "stage4","income1","income3","income45")
W.temp <- na.omit(W.temp) ### LDDPdensity does not run with NA in data, stage1 is included because of missingness
y <- W.temp[,2]
W <- W.temp[,-(1:2)]


## prior specification
S0 <- 1000 * solve(t(W) %*% W)
m0 <- solve(t(W) %*% W) %*% t(W) %*% y
prior <- list(a0 = 10, b0 = 1, m0 = m0, S0 = S0, tau1 = 6.01,
              taus1 = 6.01, taus2 = 2.01, nu = 9, psiinv = solve(S0))

## covariate patterns to calculate densities for (we only need to estimate the densities in the whites) 
Wpred <- cbind(rep(1,dim(covpatterns)[1]),rep(0,dim(covpatterns)[1]),as.matrix(covpatterns))
colnames(Wpred)[1:2] <- c("int","x")

#mcmc <- list(nburn = 50, nsave = 50, nskip = 1, ndisplay = 10)
work.dir <- "/Users/katrinadevick" 

set.seed(777)
denreg <- LDDPdensity(formula = y ~ W - 1, zpred = Wpred, ngrid = 200, 
                      compute.band = TRUE, type.band = "HPD", prior = prior, mcmc = mcmc,
                      state = NULL, status = TRUE,work.dir=work.dir)

denssam <- extractcdensity.LDDPdensity(denreg)
cdfsam <- extractccdf.LDDPdensity(denreg)

save(denreg, denssam, cdfsam, file="denreg_bmi_c_200.RData") 




#################################################### 
####    calc residual disparity for den reg     ####
#################################################### 


load("denreg_bmi_c_200.RData")
load("thetassam_bmisq_BEST.RData")

## cbind the posterior theta samples and CDFs 
colnames(thetassam.bmisq.BEST) <- rownames(sum.bmisq.BEST$coefficients$beta)
posteriorsam.bmisq.BEST <- cbind(cdfsam[[1]],thetassam.bmisq.BEST)

## the function densitymediation.bmisq.BEST is currently written for the above outcome model, it will need be to be tweaked 
## if the dimension of the polyinomial effect of the mediator on the outcome or interactions included in the model change  
denmed.bmisq.BEST.10000 <- densitymediation.bmisq.BEST(den.model=denreg, covpatterns=covpatterns, prob.covpatterns=prob.covpatterns,
                                                       posteriorsam=posteriorsam.bmisq.BEST, R=10000, seed=102717)
save(denmed.bmisq.BEST.10000,file="den_med_bmisq_BEST_10000.RData")




#################################################### 
####    calc conditional residual disparity     ####
#################################################### 


load("den_med_bmisq_BEST_10000.RData")

## indicators for each age group
age1 <- as.numeric(covpatterns[,"agecat1"]==1)
age2 <- as.numeric(covpatterns[,"agecat1"]==0 & covpatterns[,"agecat3"]==0)
age3 <- as.numeric(covpatterns[,"agecat3"]==1)

## calc conditional prob of each covariate pattern given a particular age group
age1.prob <- age1*prob.covpatterns/sum(age1*prob.covpatterns)
age2.prob <- age2*prob.covpatterns/sum(age2*prob.covpatterns)
age3.prob <- age3*prob.covpatterns/sum(age3*prob.covpatterns)

## obtain posterior samples of conditional RD 
RDsam.age1 <- denmed.bmisq.BEST.1000$RDsam %*% age1.prob
RDsam.age2 <- denmed.bmisq.BEST.1000$RDsam %*% age2.prob
RDsam.age3 <- denmed.bmisq.BEST.1000$RDsam %*% age3.prob

## numbers for "Disparity" and "RD Density" rows in table 2 from paper
xtable(cbind(c(disp.bayes.age1,disp.bayes.age2,disp.bayes.age3),
             rbind(quantile(RDsam.age1,probs=c(.5,0.025,.975)),
                   quantile(RDsam.age2,probs=c(.5,0.025,.975),na.rm=TRUE),
                   quantile(RDsam.age3,probs=c(.5,0.025,.975)))),digits=3)

