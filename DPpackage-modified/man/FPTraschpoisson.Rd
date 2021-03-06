\name{FPTraschpoisson}
\alias{FPTraschpoisson}
\alias{FPTraschpoisson.default}

\title{Bayesian analysis for a Finite Polya Tree Rasch Poisson model}
\description{
    This function generates a posterior density sample for a 
    Rasch Poisson model, using a Finite Polya Tree or a Mixture of Finite Polya Tree 
    prior for the distribution of the random effects. 
}
   
\usage{

FPTraschpoisson(y,prior,mcmc,offset,state,status,
                grid=seq(-10,10,length=1000),data=sys.frame(sys.parent()),
                compute.band=FALSE)
     
}

\arguments{
    \item{y}{        a matrix giving the data for which the Rasch Poisson Model
                     is to be fitted.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Finite Polya tree
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing), \code{mub} and \code{Sb} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, \code{mu}
                     giving the mean of the normal baseline distribution
                     (is must be specified if \code{mub} and \code{Sb} are missing),
                     \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution of variance of the normal baseline distribution,
                     \code{sigma} giving the standard deviation of the normal baseline distribution
                     (is must be specified if \code{tau1} and \code{tau2} are missing),
                     \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the difficulty
                     parameters, and \code{M} giving the finite level to be considered for the 
                     Finite Polya tree.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, and \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{offset}{   this can be used to specify an a priori known component to
                     be included in the linear predictor during the fitting.}

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{grid}{     grid points where the density estimate is 
                     evaluated. The default is seq(-10,10,length=1000).}

    \item{data}{     data frame.}       

    \item{compute.band}{  logical variable indicating whether the confidence band for the density 
                          and CDF must be computed.}       
    
}

\details{
  This generic function fits a semiparametric Rasch Poisson model as in
  San Martin et al. (2011), where the linear predictor is modeled as follows:
  \deqn{\eta_{ij} = \theta_i - \beta_j, i=1,\ldots,n, j=1,\ldots,k}{etaij = thetai - betaj, i=1,\ldots,n, j=1,\ldots,k}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha,\mu,\sigma^2 \sim FPT^M(\Pi^{\mu,\sigma^2},\textit{A})}{G | alpha,mu,sigma2 ~ FPT^M(Pi^{mu,sigma2},\textit{A})}  
  where, the the PT is centered around a \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution, by
  taking each \eqn{m} level of the partition \eqn{\Pi^{\mu,\sigma^2}}{Pi^{mu,sigma2}} to coincide 
  with the \eqn{k/2^m, k=0,\ldots,2^m} quantile of the \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution.
  The family \eqn{\textit{A}=\{\alpha_e: e \in E^{*}\}}{\textit{A}=\{alphae: e \in E^{*}\}}, where \eqn{E^{*}=\bigcup_{m=0}^{\infty} E^m}{E^{*}=\bigcup_{m=0}^{+infty} E^m} 
  and \eqn{E^m} is the \eqn{m}-fold product of \eqn{E=\{0,1\}},
  was specified as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}. 
  
  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}

  Each of the parameters of the baseline distribution, \eqn{\mu}{mu} and \eqn{\sigma}{sigma} 
  can be considered as random or fixed at some particular value. In the 
  first case, a Mixture of Polya Trees Process is considered as a prior
  for the distribution of the random effects. To let \eqn{\sigma^2}{sigma2} to be fixed at a particular
  value, set \eqn{\tau_1}{tau1} to NULL in the prior specification. 
  To let \eqn{\mu}{mu} to be fixed at a particular
  value, set \eqn{\mu_b}{mub} to NULL in the prior specification.

  In the computational implementation of the model, a Metropolis-Hastings step is used to sample
  the full conditional of the difficulty parameters. The full conditionals for abilities and PT
  parameters are sampled using slice sampling. We refer to Jara, Hanson and Lesaffre (2009) for more details and for 
  the description regarding sampling functionals of PTs. 

}

\value{
  An object of class \code{FPTraschpoisson} representing the Rasch Poisson
  model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{beta}, \code{mu}, \code{sigma2}, and the precision parameter 
  \code{alpha}.
  
  The function \code{DPrandom} can be used to extract the posterior mean of the 
  random effects.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{alpha}{ giving the value of the precision parameter.} 
  
  \item{b}{ a vector of dimension nsubjects giving the value of the random effects
  for each subject.} 

  \item{beta}{ giving the value of the difficulty parameters.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.}
  
  \item{sigma2}{ giving the variance of the normal baseline distributions.}
}

\seealso{
\code{\link{DPrandom}}, \code{\link{DPraschpoisson}}
}

\references{
Hanson, T., Johnson, W. (2002) Modeling regression error with a Mixture
  of Polya Trees. Journal of the American Statistical Association, 
  97: 1020 - 1033.

Jara, A., Hanson, T.,  Lesaffre, E. (2009) Robustifying Generalized Linear 
   Mixed Models using a New Class of Mixture of Multivariate Polya Trees. Journal of
   Computational and Graphical Statistics, 18(4): 838-860.

Lavine, M. (1992) Some aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 20: 1222-11235.

Lavine, M. (1994) More aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 22: 1161-1176.

San Martin, E., Jara, A., Rolin, J.-M., and Mouchart, M. (2011) On the Bayesian 
  nonparametric generalization of IRT-type models. Psychometrika (To appear).
  
}

\examples{
\dontrun{
    ####################################
    # A simulated Data Set
    ####################################
      nsubject <- 200
      nitem <- 10
      
      y <- matrix(0,nrow=nsubject,ncol=nitem)
      
      ind <- rbinom(nsubject,1,0.5)
      theta <- ind*rnorm(nsubject,1,0.25)+(1-ind)*rnorm(nsubject,3,0.25)
      beta <- c(0,seq(-1,1,length=nitem-1))
      true.density <- function(grid)
      {
         0.5*dnorm(grid,1,0.25)+0.5*dnorm(grid,3,0.25) 
      }  

      for(i in 1:nsubject)
      {
         for(j in 1:nitem)
         {
            eta<-theta[i]-beta[j]         
            mean<-exp(eta)
            y[i,j]<-rpois(1,mean)
         }
      }
 
    # Prior information

      beta0 <- rep(0,nitem-1)
      Sbeta0 <- diag(1000,nitem-1)

      prior <- list(alpha=1,
                    tau1=6.01,
                    tau2=2.01,
                    mub=0,
                    Sb=100,
                    beta0=beta0,
                    Sbeta0=Sbeta0,
                    M=5)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 5000
      nskip <- 0
      ndisplay <- 100
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fit the model
      fit1 <- FPTraschpoisson(y=y,prior=prior,mcmc=mcmc,
                              state=state,status=TRUE,
                              grid=seq(-1,5,0.01),
                              compute.band=TRUE)

    # Density estimate (along with HPD band) and truth
      plot(fit1$grid,fit1$dens.u,lwd=2,col="blue",type="l",lty=2,
           xlab=expression(theta),ylab="density")
      lines(fit1$grid,fit1$dens,lwd=2,col="blue")
      lines(fit1$grid,fit1$dens.l,lwd=2,col="blue",lty=2)
      lines(fit1$grid,true.density(fit1$grid),col="red")

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

    # Extract random effects
      DPrandom(fit1)
      plot(DPrandom(fit1))
      DPcaterpillar(DPrandom(fit1))

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
