
###############################################
## R function to extract the samples for the
## consitional densities from a LDDP object
###############################################

extractcdensity.LDDPdensity <- function(object)
{
  if(is(object, "LDDPdensity"))
  {
    work.dir <- object$work.dir
    if(!is.null(work.dir))
    {
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    
    nsave <- object$mcmc$nsave
    ngrid <- object$ngrid
    npred <- object$npred
    
    denspw <- matrix(0,nrow=npred,ncol=ngrid)
    
    denssam <- matrix(0,nrow=nsave,ncol=(npred*ngrid))
    
    foo <- .Fortran("readlddpdens2",
                    nsave     = as.integer(nsave), 
                    npred     = as.integer(npred),
                    ngrid     = as.integer(ngrid),
                    denspw    = as.double(denspw),
                    denssam   = as.double(denssam),
                    PACKAGE="DPpackage")
    
    denssam <- matrix(foo$denssam,nrow=nsave,ncol=npred*ngrid)
    
    if(!is.null(work.dir))
    {
      setwd(old.dir)
    }
    
    out <- list(denssam=denssam)
    return(out)
  }
}





###############################################
## R function to extract the samples for the
## conditional cdf from a LDDP object
###############################################

extractccdf.LDDPdensity <- function(object)
{
  if(is(object, "LDDPdensity"))
  {
    
    work.dir <- object$work.dir
    if(!is.null(work.dir))
    {
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    
    
    nsave <- object$mcmc$nsave
    ngrid <- object$ngrid
    npred <- object$npred
    
    cdfpw <- matrix(0,nrow=npred,ncol=ngrid)
    
    cdfsam <- matrix(0,nrow=nsave,ncol=(npred*ngrid))
    
    foo <- .Fortran("readlddpcdf2",
                    nsave     = as.integer(nsave), 
                    npred     = as.integer(npred),
                    ngrid     = as.integer(ngrid),
                    cdfpw     = as.double(cdfpw),
                    cdfsam    = as.double(cdfsam),
                    PACKAGE="DPpackage")
    
    cdfsam <- matrix(foo$cdfsam,nrow=nsave,ncol=npred*ngrid)
    
    if(!is.null(work.dir))
    {
      setwd(old.dir)
    }
    
    out <- list(cdfsam=cdfsam)
    return(out)
  }
}



###############################################
## Example
###############################################

library(DPpackage)
setwd("/Users/katrinadevick")

# Simulated data

nobs <- 500
y1   <-rnorm(nobs, 3,.8)
y21 <- rnorm(nobs,1.5, 0.8)
y22 <- rnorm(nobs,4.0, 0.6)
u <- runif(nobs)
y2 <- ifelse(u<0.6,y21,y22)
y <- c(y1,y2)

# design matrix including a single factor
trt <- c(rep(0,nobs),rep(1,nobs))

# design matrix for posterior predictive 
zpred <- rbind(c(1,0),c(1,1))  

# Prior information

S0 <- diag(100,2)
m0 <- rep(0,2)
psiinv <- diag(1,2)

prior <- list(a0=10,b0=1,nu=4,m0=m0,S0=S0,psiinv=psiinv,
              tau1=6.01,taus1=6.01,taus2=2.01)

# Initial state
state <- NULL

# MCMC parameters
nburn <- 5
nsave <- 40
nskip <- 0
ndisplay <- 20
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)

# Working directory (please change this accordingly)
work.dir <- "/Users/katrinadevick"

# Fitting the model
fit1 <- LDDPdensity(y~trt,prior=prior,mcmc=mcmc,state=state,status=TRUE,
                    ngrid=200,zpred=zpred,work.dir=work.dir)

# Using the R function to extract the samples

denssam <- extractcdensity.LDDPdensity(fit1)
cdfsam <- extractccdf.LDDPdensity(fit1)

