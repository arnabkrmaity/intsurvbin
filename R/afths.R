#' Function to implement the horseshoe shrinkage prior in Bayesian linear regression
#'
#'
#' This function employs the algorithm proposed in Bhattacharya et al. (2015).
#' The global-local scale parameters are updated via a slice sampling scheme given
#'  in the online supplement of Polson et al. (2014). Two different algorithms are
#'  used to compute posterior samples of the \eqn{p*1} vector of regression coefficients \eqn{\beta}.
#'  The method proposed in Bhattacharya et al. (2015) is used when \eqn{p>n}, and the
#'  algorithm provided in Rue (2001) is used for the  case \eqn{p<=n}. The function
#'  includes options for full hierarchical Bayes versions with hyperpriors on all
#'  parameters, or empirical Bayes versions where some parameters are taken equal
#'  to a user-selected value.
#'
#'  The model is:
#' \deqn{y=X\beta+\epsilon, \epsilon \sim N(0,\sigma^2)}
#'
#' The full Bayes version of the horseshoe, with hyperpriors on both \eqn{\tau} and \eqn{\sigma^2} is:
#'
#' \deqn{\beta_j \sim N(0,\sigma^2 \lambda_j^2 \tau^2)}
#' \deqn{\lambda_j \sim Half-Cauchy(0,1), \tau \sim Half-Cauchy (0,1)}
#' \deqn{\sigma^2 \sim 1/\sigma^2}
#'
#' There is an option for a truncated Half-Cauchy prior (truncated to [1/p, 1]) on \eqn{\tau}.
#' Empirical Bayes versions are available as well, where \eqn{\tau} and/or
#' \eqn{\sigma^2} are taken equal to fixed values, possibly estimated using the data.
#'
#' @references Bhattacharya, A., Chakraborty, A. and Mallick, B.K. (2015), Fast Sampling
#' with Gaussian Scale-Mixture priors in High-Dimensional Regression.
#'
#' Polson, N.G., Scott, J.G. and Windle, J. (2014) The Bayesian Bridge.
#' Journal of Royal Statistical Society, B, 76(4), 713-733.
#'
#' Rue, H. (2001). Fast sampling of Gaussian Markov random fields. Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 63, 325--338.
#'
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010), The Horseshoe
#'  Estimator for Sparse Signals. Biometrika 97(2), 465--480.
#'
#'@param ct survival response, a \eqn{n*2} matrix with first column as response and second column as right censored indicator,
#'1 is event time and 0 is right censored.
#'@param X Matrix of covariates, dimension \eqn{n*p}.
#'@param method.tau Method for handling \eqn{\tau}. Select "truncatedCauchy" for full
#' Bayes with the Cauchy prior truncated to [1/p, 1], "halfCauchy" for full Bayes with
#' the half-Cauchy prior, or "fixed" to use a fixed value (an empirical Bayes estimate,
#' for example).
#'@param tau  Use this argument to pass the (estimated) value of \eqn{\tau} in case "fixed"
#' is selected for method.tau. Not necessary when method.tau is equal to"halfCauchy" or
#' "truncatedCauchy". The default (tau = 1) is not suitable for most purposes and should be replaced.
#'@param method.sigma Select "Jeffreys" for full Bayes with Jeffrey's prior on the error
#'variance \eqn{\sigma^2}, or "fixed" to use a fixed value (an empirical Bayes
#'estimate, for example).
#'@param Sigma2 A fixed value for the error variance \eqn{\sigma^2}. Not necessary
#'when method.sigma is equal to "Jeffreys". Use this argument to pass the (estimated)
#'value of Sigma2 in case "fixed" is selected for method.sigma. The default (Sigma2 = 1)
#'is not suitable for most purposes and should be replaced.
#'@param burn Number of burn-in MCMC samples. Default is 1000.
#'@param nmc Number of posterior draws to be saved. Default is 5000.
#'@param thin Thinning parameter of the chain. Default is 1 (no thinning).
#'@param alpha Level for the credible intervals. For example, alpha = 0.05 results in
#'95\% credible intervals.
#'@param Xtest test design matrix.
#'@param cttest test survival response.

#'
#'@return \item{BetaHat}{Posterior mean of Beta, a \eqn{p} by 1 vector.}
#' \item{LeftCI}{The left bounds of the credible intervals.}
#' \item{RightCI}{The right bounds of the credible intervals.}
#' \item{BetaMedian}{Posterior median of Beta, a \eqn{p} by 1 vector.}
#' \item{Sigma2Hat}{Posterior mean of error variance \eqn{\sigma^2}. If method.sigma =
#' "fixed" is used, this value will be equal to the user-selected value of Sigma2
#' passed to the function.}
#' \item{TauHat}{Posterior mean of global scale parameter tau, a positive scalar.
#' If method.tau = "fixed" is used, this value will be equal to the user-selected value
#' of tau passed to the function.}
#' \item{BetaSamples}{Posterior samples of Beta.}
#' \item{TauSamples}{Posterior samples of tau.}
#' \item{Sigma2Samples}{Posterior samples of Sigma2.}
#'
#'
#' @examples
#' \dontrun{#In this example, there are no relevant predictors
#' #20 observations, 30 predictors (betas)
#' y <- rnorm(20)
#' X <- matrix(rnorm(20*30) , 20)
#' res <- horseshoe(y, X, method.tau = "truncatedCauchy", method.sigma = "Jeffreys")
#'
#' plot(y, X%*%res$BetaHat) #plot predicted values against the observed data
#' res$TauHat #posterior mean of tau
#' HS.var.select(res, y, method = "intervals") #selected betas
#' #Ideally, none of the betas is selected (all zeros)
#' #Plot the credible intervals
#' library(Hmisc)
#' xYplot(Cbind(res$BetaHat, res$LeftCI, res$RightCI) ~ 1:30)
#'}
#'
#' \dontrun{ #The horseshoe applied to the sparse normal means problem
#' # (note that HS.normal.means is much faster in this case)
#' X <- diag(100)
#' beta <- c(rep(0, 80), rep(8, 20))
#' y <- beta + rnorm(100)
#' res2 <- horseshoe(y, X, method.tau = "truncatedCauchy", method.sigma = "Jeffreys")
#' #Plot predicted values against the observed data (signals in blue)
#' plot(y, X%*%res2$BetaHat, col = c(rep("black", 80), rep("blue", 20)))
#' res2$TauHat #posterior mean of tau
#' HS.var.select(res2, y, method = "intervals") #selected betas
#' #Ideally, the final 20 predictors are selected
#' #Plot the credible intervals
#' library(Hmisc)
#' xYplot(Cbind(res2$BetaHat, res2$LeftCI, res2$RightCI) ~ 1:100)
#'}
#'
#' @export


# 15 June 2016
# Use Polya Gamma data augmentation approach instead of Albert and Chib (1993) 


afths <- function(ct, X, method.tau = c("fixed", "truncatedCauchy","halfCauchy"), tau = 1,
                  method.sigma = c("fixed", "Jeffreys"), Sigma2 = 1,
                  burn = 1000, nmc = 5000, thin = 1, alpha = 0.05,
                  Xtest = NULL, cttest = NULL)
{
  
  method.tau = match.arg(method.tau)
  
  method.sigma = match.arg(method.sigma)
  
  ptm=proc.time()
  niter <- burn+nmc
  effsamp=(niter -burn)/thin
  
  
  # X <- as.matrix(bdiag(X, X))
  n=nrow(X)
  p <- ncol(X)
  if(is.null(Xtest))
  {
    Xtest <- X
    ntest <- n
    # cttest<- ct
  } else {
    ntest <- nrow(Xtest)
  }
  
  
 
  time         <- ct[, 1]
  status       <- ct[, 2]
  censored.id  <- which(status == 0)
  n.censored   <- length(censored.id)  # number of censored observations
  X.censored   <- X[censored.id, ]
  X.observed   <- X[-censored.id, ]
  y.s <- logtime <- log(time)   # for coding convenience, since the whole code is written with y
  y.s.censored <- y.s[censored.id]
  y.s.observed <- y.s[-censored.id]
  
  
  timetest         <- cttest[, 1]
  statustest       <- cttest[, 2]
  censored.idtest  <- which(statustest == 0)
  n.censoredtest   <- length(censored.idtest)  # number of censored observations
  y.stest <- logtimetest <- log(timetest)   # for coding convenience, since the whole code is written with y
  
  
  ## parameters ##
  beta.s   <- rep(0, p);
  lambda   <- rep(1, p);
  sigma_sq <- Sigma2;
  adapt.par  <- c(100, 20, 0.5, 0.75)
  prop.sigma <- 1
  
  
  ## output ##
  beta.sout          <- matrix(0, p, effsamp)
  lambdaout          <- matrix(0, p, effsamp)
  tauout             <- rep(0, effsamp)
  sigmaSqout         <- rep(1, effsamp)
  predsurvout        <- matrix(0, ntest, effsamp)
  logtimeout         <- matrix(0, ntest, effsamp)
  loglikelihood.sout <- rep(0, effsamp)
  likelihood.sout    <- matrix(0, n, effsamp)
  y.sout             <- matrix(0, n, effsamp)
  trace              <- array(dim = c(niter, length(sigma_sq)))
  
  
  
  ## start Gibb's sampling ##
  for(i in 1:niter)
  {
    
    ################################
    ####### survival part ##########
    ################################
    
    mean.impute <- X.censored %*% beta.s
    sd.impute   <- sqrt(sigma_sq)
    ## update censored data ##
    time.censored <- msm::rtnorm(n.censored, mean = mean.impute, sd = sd.impute, lower = y.s.censored)
    # truncated at log(time) for censored data
    y.s[censored.id] <- time.censored
    
    
    
    ## which algo to use ##
    if(p > n)
    {
      bs  = bayesreg.sample_beta(X, z = y.s, mvnrue = FALSE, b0 = 0, sigma2 = sigma_sq, tau2 = tau^2, 
                                 lambda2 = lambda^2, omega2  = matrix(1, n, 1), XtX = NA)
      beta.s <- bs$x
      
    } else {
      
      bs  = bayesreg.sample_beta(X, z = y.s, mvnrue = TRUE, b0 = 0, sigma2 = sigma_sq, tau2 = tau^2, 
                                 lambda2 = lambda^2, omega2  = matrix(1, n, 1), XtX = NA)
      beta.s <- bs$x
    }
    
    
    
    # ## update sigma_sq ##
    # if(method.sigma == "Jeffreys"){
    #   
    #   E_1=max(t(y.s-X%*%beta.s)%*%(y.s-X%*%beta.s),(1e-10))
    #   E_2=max(sum(beta.s^2/(tau*lambda)^2),(1e-10))
    #   
    #   sigma_sq= 1/stats::rgamma(1, (n + p)/2, scale = 2/(E_1+E_2))
    # }
    
    ## Update Sigma using Metropolis Hastings scheme
    sigma2.current <- sigma_sq
    theta.current  <- log(sigma2.current)
    trace[i, ]     <- theta.current
    
    # adjust sigma for proposal distribution (taken from Metro_Hastings() function)
    if (i > adapt.par[1] && i%%adapt.par[2] == 0 && i < (adapt.par[4] * niter))
    {
      len <- floor(i * adapt.par[3]):i
      t   <- trace[len, ]
      nl  <- length(len)
      p.sigma <- (nl - 1) * stats::var(t)/nl
      p.sigma <- MHadaptive::makePositiveDefinite(p.sigma)
      if (!(0 %in% p.sigma))
        prop.sigma <- p.sigma
    }
    
    theta.proposed  <- stats::rnorm(1, mean = theta.current, sd = sqrt(prop.sigma))
    sigma2.proposed <- exp(theta.proposed)
    y.s.tilde       <- y.s - X %*% beta.s
    kernel          <- min(1, exp(sum(stats::dnorm(y.s.tilde, sd = sqrt(sigma2.proposed), log = TRUE)) -
                                  sum(stats::dnorm(y.s.tilde, sd = sqrt(sigma2.current), log = TRUE)) +
                                  stats::dnorm(theta.proposed, log = TRUE) - 
                                  stats::dnorm(theta.current,  log = TRUE) + 
                                    abs(theta.proposed) - abs(theta.current)))
                                    # stats::dnorm(theta.current, log = TRUE) - 
                                    # stats::dnorm(theta.proposed, log = TRUE))) 
    if(stats::rbinom(1, size = 1, kernel) == 1)
    {
      theta.current <- theta.proposed
      sigma2.current <- sigma2.proposed
    }
    sigma_sq <- sigma2.current 
    
    
    
    
    beta <- c(beta.s)  # all \beta's together
    Beta <- matrix(beta, ncol = p, byrow = TRUE)
    Beta[1, ] <- Beta[1, ]/sqrt(sigma_sq)
    
    ## update lambda_j's in a block using slice sampling ##
    eta = 1/(lambda^2)
    upsi = stats::runif(p,0,1/(1+eta))
    tempps = apply(Beta^2, 2, sum)/(2*tau^2)
    ub = (1-upsi)/upsi
    # now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = stats::pgamma(ub, (1 + 1)/2, scale = 1/tempps)
    Fub[Fub < (1e-4)] = 1e-4;  # for numerical stability
    up = stats::runif(p,0,Fub)
    eta <- stats::qgamma(up, (1 + 1)/2, scale=1/tempps)
    lambda = 1/sqrt(eta);
    
    ## update tau ##
    ## Only if prior on tau is used
    if(method.tau == "halfCauchy"){
      tempt <- sum(apply(Beta^2, 2, sum)/(2*lambda^2))
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt = (1-utau)/utau
      Fubt = stats::pgamma(ubt,(p+1)/2,scale=1/tempt)
      Fubt = max(Fubt,1e-8) # for numerical stability
      ut = stats::runif(1,0,Fubt)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = 1/sqrt(et)
    }#end if
    
    if(method.tau == "truncatedCauchy"){
      tempt <- sum(apply(Beta^2, 2, sum)/(2*lambda^2))
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt_1=1
      ubt_2 = min((1-utau)/utau,p^2)
      Fubt_1 = stats::pgamma(ubt_1,(p+1)/2,scale=1/tempt)
      Fubt_2 = stats::pgamma(ubt_2,(p+1)/2,scale=1/tempt)
      #Fubt = max(Fubt,1e-8) # for numerical stability
      ut = stats::runif(1,Fubt_1,Fubt_2)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = 1/sqrt(et)
    }
    
    
    mean <- Xtest %*% beta.s
    sd   <- sqrt(sigma_sq)
    predictive.survivor <- stats::pnorm(mean/sd, lower.tail = FALSE)
    logt <- Xtest %*% beta.s
    
    
    
    ## Prediction ##
    mean.stest          <- Xtest %*% beta.s
    sd.stest            <- sqrt(sigma_sq)
    predictive.survivor <- stats::pnorm((y.stest - mean.stest)/sd.stest, lower.tail = FALSE)
    logt                <- mean.stest
    
    ## Following is required for DIC computation
    loglikelihood.s <- sum(c(stats::dnorm(y.s.observed, mean = X.observed %*% beta.s, sd = sqrt(sigma_sq),
                                   log = TRUE),
                             log(1 - stats::pnorm(y.s.censored, mean = X.censored %*% beta.s,
                                           sd = sqrt(sigma_sq)))))
    # loglikelihood.s <- sum(c(stats::dnorm(y.s[-censored.id], mean = X.observed %*% beta.s, sd = sqrt(sigma_sq), 
    #                                log = TRUE),
    #                          dtnorm(y.s[censored.id], mean = X.censored %*% beta.s, sd = sqrt(sigma_sq),
    #                                 lower = y.s.censored, log = TRUE)))
    
    
    ## Following is required for DIC computation
    likelihood.s <- c(stats::dnorm(y.s.observed, mean = X.observed %*% beta.s, sd = sqrt(sigma_sq)),
                      stats::pnorm(y.s.censored, mean = X.censored %*% beta.s, sd = sqrt(sigma_sq),
                            lower.tail = FALSE))
    
    
    if (i%%500 == 0)
    {
      print(i)
    }
    
   
    
    if(i > burn && i%%thin== 0)
    {
      beta.sout[ ,(i-burn)/thin]     <- beta.s
      lambdaout[ ,(i-burn)/thin]     <- lambda
      tauout[(i - burn)/thin]        <- tau
      sigmaSqout[(i - burn)/thin]    <- sigma_sq
      predsurvout[ ,(i - burn)/thin] <- predictive.survivor
      logtimeout[, (i - burn)/thin]  <- logt
      loglikelihood.sout[(i - burn)/thin] <- loglikelihood.s
      likelihood.sout[, (i - burn)/thin]  <- likelihood.s
      y.sout[, (i - burn)/thin] <- y.s
    }
  }
  
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  pMean.s   <- apply(beta.sout,1, mean)
  pMedian.s <- apply(beta.sout,1,stats::median)
  pLambda   <- apply(lambdaout, 1, mean)
  pSigma= mean(sigmaSqout)
  pTau=mean(tauout)
  pPS <- apply(predsurvout, 1, mean)
  pLogtime  <- apply(logtimeout, 1, mean)
  pLoglikelihood.s <- mean(loglikelihood.sout)
  pLikelihood.s    <- apply(likelihood.sout, 1, mean)
  py.s             <- apply(y.sout, 1, mean)
  
  
  loglikelihood.posterior.s <- sum(c(stats::dnorm(y.s.observed, mean = X.observed %*% pMean.s, sd = sqrt(pSigma), log = TRUE),
                                     log(1 - stats::pnorm(y.s.censored, mean = X.censored %*% pMean.s,
                                               sd = sqrt(pSigma)))))
  # loglikelihood.posterior.s <- sum(c(stats::dnorm(py.s[-censored.id], mean = X.observed %*% pMean.s, sd = sqrt(pSigma), 
  #                                log = TRUE),
  #                          stats::dtnorm(py.s[censored.id], mean = X.censored %*% pMean.s, sd = sqrt(pSigma),
  #                                 lower = y.s.censored, log = TRUE)))
  DIC.s <- -4 * pLoglikelihood.s + 2 * loglikelihood.posterior.s
  WAIC <- -2 * (sum(log(pLikelihood.s)) - 
                  2 * (sum(log(pLikelihood.s)) - pLoglikelihood.s))
  
  
  #construct credible sets
  left <- floor(alpha*effsamp/2)
  right <- ceiling((1-alpha/2)*effsamp)
  
  beta.sSort <- apply(beta.sout, 1, sort, decreasing = F)
  left.spoints <- beta.sSort[left, ]
  right.spoints <- beta.sSort[right, ]
  
  
  result=list("SurvivalHat" = pPS, "LogTimeHat" = pLogtime, "Beta.sHat"= pMean.s, 
              "LeftCI.s" = left.spoints, "LeftCI.s" = right.spoints,
              "Beta.sMedian" = pMedian.s, 
              "LambdaHat" = pLambda,
              "Sigma2Hat"=pSigma,"TauHat"=pTau, "Beta.sSamples" = beta.sout, 
              "TauSamples" = tauout, "Sigma2Samples" = sigmaSqout, "DIC.s" = DIC.s)
  return(result)
}




# ============================================================================================================================
# Sample the regression coefficients
bayesreg.sample_beta <- function(X, z, mvnrue, b0, sigma2, tau2, lambda2, omega2, XtX)
{
  alpha  = (z - b0)
  Lambda = sigma2 * tau2 * lambda2
  sigma  = sqrt(sigma2)
  
  # Use Rue's algorithm
  if (mvnrue)
  {
    # If XtX is not precomputed
    if (any(is.na(XtX)))
    {
      omega = sqrt(omega2)
      X0    = apply(X,2,function(x)(x/omega))
      bs    = bayesreg.fastmvg2_rue(X0/sigma, alpha/sigma/omega, Lambda)
    }
    
    # XtX is precomputed (Gaussian only)
    else {
      bs    = bayesreg.fastmvg2_rue(X/sigma, alpha/sigma, Lambda, XtX/sigma2)
    }
  }
  
  # Else use Bhat. algorithm
  else
  {
    omega = sqrt(omega2)
    X0    = apply(X,2,function(x)(x/omega))
    bs    = bayesreg.fastmvg_bhat(X0/sigma, alpha/sigma/omega, Lambda)
  }
  
  return(bs)
}


# ============================================================================================================================
# function to generate multivariate normal random variates using Rue's algorithm
bayesreg.fastmvg2_rue <- function(Phi, alpha, d, PtP = NA)
{
  Phi   = as.matrix(Phi)
  alpha = as.matrix(alpha)
  r     = list()
  
  # If PtP not precomputed
  if (any(is.na(PtP)))
  {
    PtP = t(Phi) %*% Phi
  }
  
  p     = ncol(Phi)
  if (length(d) > 1)
  {
    Dinv  = diag(as.vector(1/d))
  }
  else
  {
    Dinv   = 1/d
  }
  L     = t(chol(PtP + Dinv))
  v     = forwardsolve(L, t(Phi) %*% alpha)
  r$m   = backsolve(t(L), v)
  w     = backsolve(t(L), stats::rnorm(p,0,1))
  
  r$x   = r$m + w
  return(r)
}


# ============================================================================================================================
# function to generate multivariate normal random variates using Bhat. algorithm
bayesreg.fastmvg_bhat <- function(Phi, alpha, d)
{
  d     = as.matrix(d)
  p     = ncol(Phi)
  n     = nrow(Phi)
  r     = list()
  
  u     = as.matrix(stats::rnorm(p,0,1)) * sqrt(d)
  delta = as.matrix(stats::rnorm(n,0,1))
  v     = Phi %*% u + delta
  Dpt   = (apply(Phi, 1, function(x)(x*d)))
  W     = Phi %*% Dpt + diag(1,n)
  w     = solve(W,(alpha-v))
  
  r$x   = u + Dpt %*% w
  r$m   = r$x
  
  return(r)
}
