###########################################
##   source code for density mediation   ##
##   developed by: Katrina Devick        ##
##   last updated: 8 May 2019            ##
###########################################


invcdf <- function(p,cdf,grid){ 
  ### linear interpolation to find grid value from inverse cdf sampling
  index <- which(abs(cdf-p)==min(abs(cdf-p)))
  if(index == 1){cdf.0 <- 0} else {cdf.0 <- cdf[index-1]}
  if(is.na(cdf[index+1])) {cdf.1 <- 1} else {cdf.1 <- cdf[index+1]} ### using this since values 
  
  y <- p-cdf[index]
  if(y >= 0){
    slope <- (cdf.1-cdf[index])/(diff(grid)[1])
  }else {
    slope <- (cdf[index]-cdf.0)/(diff(grid)[1])
  }
  toreturn <- grid[index]+y/slope
  # toreturn
  return(toreturn)
}



int.bmisq <- function(cdf.row,grid,R,theta21,theta22,theta31,theta32,ngrid,a){ 
  p <- runif(R)  
  randomsamples <- sapply(p,invcdf,cdf=cdf.row,grid=grid)
  
  expestimates <- exp((theta21+theta31*a)*(randomsamples) + (theta22+theta32*a)*(randomsamples)^2)  
  return(mean(expestimates)) ### this is the estimate of the int in the residual disparity 
  #### calculation for a particular covariate pattern.. need to run this serperate for the 
  #### blacks and whites to calcule resdisp 
}





firstterm.BEST <- function(covpatterns.row, theta1, theta6){
  age <- covpatterns.row[c("agecat1","agecat3")]
  
  term <- exp(theta1+theta6%*%age)
  return(term)
}




RDforMCMCsam.bmisq.BEST <- function(row, ngrid=ngrid, npred=npred, R=R, grid=grid, covpatterns=covpatterns){
  cdf <- matrix(row[1:(ngrid*npred)],nrow=npred,ncol=ngrid,byrow=TRUE)
  theta <- row[(ngrid*npred+1):length(row)] 
  theta1 <- theta["race"]
  theta21 <- theta["bmi.c"]
  theta22 <- theta["bmisq.c"]
  theta31 <- theta["race:bmi.c"]
  theta32 <- theta["race:bmisq.c"]
  theta6  <- theta[c("race:agecat1","race:agecat3")]
  
  int.black.shifted <- apply(cdf,1,int.bmisq,grid=grid,R=R,theta21=theta21,theta22=theta22,theta31=theta31,theta32=theta32,a=1) 
  int.white         <- apply(cdf,1,int.bmisq,grid=grid,R=R,theta21=theta21,theta22=theta22,theta31=theta31,theta32=theta32,a=0)
  expterm           <- apply(covpatterns,1,firstterm.BEST,theta1=theta1,theta6=theta6)
  
  resdisp <- expterm*(int.black.shifted/int.white)
  return(resdisp)
}




sumRD <- function(x){
  est <- mean(x)
  med <- median(x)
  CI <- quantile(x,probs=c(0.025,0.975),na.rm=TRUE)
  
  toreturn <- c(est,med,CI)
  names(toreturn) <- c("RD Mean", "RD Median", "2.5% Percentile", "97.5% Percentile")
  return(toreturn)
}






#### covpatterns is a npred x ncov matrix
#### prob.covpatterns is a npred x 1 matrix
densitymediation.bmisq.BEST <- function(den.model, covpatterns, prob.covpatterns, posteriorsam, R, seed){
  print("beg")
  set.seed(seed)
  ngrid  <- den.model$ngrid
  npred  <- den.model$npred
  nsave  <- den.model$mcmc$nsave
  grid   <- den.model$grid
  den.m  <- den.model$densp.m
  cdf.m  <- den.model$cdfp.m
  
  print("int")
  RDsam <- t(apply(posteriorsam,1,RDforMCMCsam.bmisq.BEST,ngrid=ngrid, npred=npred, R=R, grid=grid, covpatterns=covpatterns)) 
  print("rdsam")
  RDmargsam <- RDsam %*% prob.covpatterns ##### 3000x1 matrix
  
  print("marg rd")
  summaryRD <- apply(RDsam,2,sumRD)
  summaryRDmarg <- sumRD(RDmargsam)
  
  print("toreturn")
  toreturn <- list(RDsam,RDmargsam,summaryRD,summaryRDmarg)
  names(toreturn) <- c("RDsam","margRDsam","summaryRD","summaryRDmarg") 
  return(toreturn)
}




###############################################
## R function to extract the samples for the
## consitional densities from a LDDP object
## function developed by Alejandro Jara
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
## function developed by Alejandro Jara
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

