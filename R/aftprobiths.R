#' Function to implement the horseshoe shrinkage prior in integrated survival and binary regression
#'
#'
#'
#' @references Maity, A. K., Carroll, R. J., and Mallick, B. K. (2018),
#' Integration of Survival and Binary Data for Variable Selection and Prediction: A Bayesian Approach.
#'
#'
#'@param ct survival response, a \eqn{n*2} matrix with first column as response and second column as right censored indicator,
#'1 is event time and 0 is right censored.
#'@param z binary response, a \eqn{n*1} vector with numeric values 0 or 1.
#'@param X Matrix of covariates, dimension \eqn{n*p}.
#'@param burn Number of burn-in MCMC samples. Default is 1000.
#'@param nmc Number of posterior draws to be saved. Default is 5000.
#'@param thin Thinning parameter of the chain. Default is 1 (no thinning).
#'@param alpha Level for the credible intervals. For example, alpha = 0.05 results in
#'95\% credible intervals.
#'@param Xtest test design matrix.
#'@param cttest test survival response.
#'@param ztest test binary response.
#'
#'
#'
#'
#'@return \item{Beta.sHat}{Posterior mean of \eqn{\beta} for survival model, a \eqn{p} by 1 vector.}
#'\item{Beta.bHat}{Posterior mean of \eqn{\beta} for binary model, a \eqn{p} by 1 vector.}
#'\item{LeftCI.s}{The left bounds of the credible intervals for Beta.sHat.}
#'\item{RightCI.s}{The right bounds of the credible intervals for Beta.sHat.}
#'\item{LeftCI.b}{The left bounds of the credible intervals for Beta.bHat.}
#'\item{RightCI.b}{The right bounds of the credible intervals for Beta.bHat.}
#'\item{Beta.sMedian}{Posterior median of \eqn{beta} for survival model, a \eqn{p} by 1 vector.}
#'\item{Beta.bMedian}{Posterior median of \eqn{beta} for binary model, a \eqn{p} by 1 vector.}}
#'\item{SigmaHat}{Posterior mean of variance covariance matrix.}
#'\item{LambdaHat}{Posterior mean of \eqn{\lambda}, a \eqn{p*1} vector.}
#'\item{TauHat} = {Posterior mean of \eqn{\tau}, a \eqn{2*1} vector.}
#'\item{Beta.sSamples}{Posterior samples of \eqn{\beta} for survival model.}
#'\item{Beta.bSamples}{Posterior samples of \eqn{\beta} for binary model.}
#'\item{LambdaSamples}{Posterior samples of \eqn{\lambda}.}
#'\item{TauSamples}{Posterior samples of \eqn{\tau}.}
#'\item{SigmaSamples}{Posterior samples of variance covariance matrix.}
#'\item{DIC.s}{DIC for survival model.}
#'\item{DIC.b}{DIC for binary model.}
#'\item{SurvivalHat}{Predictive survival probability.}
#'\item{LogTimeHat}{Predictive log time.}
#'
#'
#'
#' @examples
#' \dontrun{ #
#'}
#'
#' @export



aftprobiths <- function(ct, z, X, burn = 1000, nmc = 5000, thin = 1, alpha = 0.05,
                         Xtest = NULL, cttest = NULL, ztest = NULL)
{

  niter=burn+nmc
  effsamp=(niter - burn)/thin

  # Survival response
  time         <- ct[, 1]
  status       <- ct[, 2]
  censored.id  <- which(status == 0)
  n.censored   <- length(censored.id)  # number of censored observations
  X.censored   <- X[censored.id, ]
  X.observed   <- X[-censored.id, ]
  y.s <- logtime <- log(time)   # for coding convenience, since the whole code is written with y
  y.s.censored <- y.s[censored.id]
  y.s.observed <- y.s[-censored.id]

  # Binary response
  w    <- z   # we switch to w because the remaining code is written as y so we set y = z, latent variable
  id1  <- which(w == 1)  # index of y where y = 1
  id0  <- which(w == 0)  # index of y where y = 0
  nid1 <- length(id1)    # number of samples where y = 1
  nid0 <- length(id0)    # number of samples where y = 0
  X1   <- X[id1, ]       # covariates where y = 1
  X0   <- X[id0, ]       # covariates where y = 0
  z1   <- z[id1]
  z0   <- z[id0]
  y.b  <- z              # for coding convenience

  p.s  <- ncol(X)
  n.s  <- nrow(X)
  Xbig <- kronecker(diag(1, 2), X)
  n    <- nrow(Xbig)
  p    <- ncol(Xbig)

  if(is.null(Xtest))
  {
    Xtest <- X
    ntest <- n
    cttest<- ct
    ztest <- z
  } else {
    ntest <- nrow(Xtest)
  }

  timetest         <- cttest[, 1]
  statustest       <- cttest[, 2]
  censored.idtest  <- which(statustest == 0)
  n.censoredtest   <- length(censored.idtest)  # number of censored observations
  Xtest.censored   <- Xtest[censored.idtest, ]
  y.stest <- logtimetest <- log(timetest)   # for coding convenience, since the whole code is written with y

  wtest    <- ztest   # we switch to w because the remaining code is written as y so we set y = z, latent variable
  id1test  <- which(wtest == 1)  # index of y where y = 1
  id0test  <- which(wtest == 0)  # index of y where y = 0
  nid1test <- length(id1test)    # number of samples where y = 1
  nid0test <- length(id0test)    # number of samples where y = 0
  X1test   <- Xtest[id1test, ]       # covariates where y = 1
  X0test   <- Xtest[id0test, ]       # covariates where y = 0
  z1test   <- ztest[id1test]
  z0test   <- ztest[id0test]
  y.btest  <- ztest              # for coding convenience



  ## parameters ##
  beta   <- rep(0, p);
  beta.s <- beta[1:p.s]        # survival coefficients
  beta.b <- beta[(p.s + 1):p]  # binary coefficients
  lambda <- rep(1, p.s);
  tau    <- rep(1, 2)
  V0     <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  v0     <- 2
  Sigma  <- V0
  theta.mean <- c(0, log(0.5))
  adapt.par  <- c(100, 20, 0.5, 0.75)
  prop.sigma <- diag(2)


  ## output ##
  betaout            <- matrix(0, p, effsamp)
  lambdaout          <- matrix(0, p.s, effsamp)
  tauout             <- matrix(0, 2, effsamp)
  Sigmaout           <- array(0, c(2, 2, effsamp))
  predsurvout        <- matrix(0, ntest, effsamp)
  logtimeout         <- matrix(0, ntest, effsamp)
  loglikelihood.sout <- rep(0, effsamp)
  loglikelihood.bout <- rep(0, effsamp)
  trace              <- array(dim = c(niter, length(theta.mean)))

  ## which algo to use ##
  if(p>n) {algo=1} else {algo=2}

  ## matrices ##
  I_n=diag(n)
  l0=rep(0,p)
  l1=rep(1,n)
  l2=rep(1,p)


  ## start Gibb's sampling ##
  for(i in 1:niter)
  {
    ## Update survival latent variable ##
    y.b.censored     <- y.b[censored.id]
    cory             <- Sigma[1, 2]/sqrt(Sigma[1, 1] * Sigma[2, 2])
    mean.impute.s    <- (X.censored %*% beta.s +
                           cory * sqrt(Sigma[1, 1]/Sigma[2, 2]) * (y.b.censored - X.censored %*% beta.b))
    sd.impute.s      <- sqrt(Sigma[1, 1] * (1 - cory^2))
    time.censored    <- msm::rtnorm(n.censored, mean = mean.impute.s, sd = sd.impute.s, lower = y.s.censored)
    y.s[censored.id] <- time.censored  # truncated at log(time) for censored data


    ## Update binary latent variable ##
    y.s1           <- y.s[id1]
    y.s0           <- y.s[id0]
    mean.impute.b1 <- X1 %*% beta.b + cory * sqrt(Sigma[2, 2]/Sigma[1, 1]) * (y.s1 - X1 %*% beta.s)
    mean.impute.b0 <- X0 %*% beta.b + cory * sqrt(Sigma[2, 2]/Sigma[1, 1]) * (y.s0 - X0 %*% beta.s)
    sd.impute.b    <- sqrt(Sigma[2, 2] * (1 - cory^2))
    y.b[id1]       <- msm::rtnorm(nid1, mean = mean.impute.b1, lower = 0)  # Equation (6) of Albert and Chib (1993)
    y.b[id0]       <- msm::rtnorm(nid0, mean = mean.impute.b0, upper = 0)

    ## Update binary latent variable for test dataset ##
    y.s1test           <- y.stest[id1test]
    y.s0test           <- y.stest[id0test]
    mean.impute.b1test <- X1test %*% beta.b + cory * sqrt(Sigma[2, 2]/Sigma[1, 1]) * (y.s1test - X1test %*% beta.s)
    mean.impute.b0test <- X0test %*% beta.b + cory * sqrt(Sigma[2, 2]/Sigma[1, 1]) * (y.s0test - X0test %*% beta.s)
    sd.impute.b    <- sqrt(Sigma[2, 2] * (1 - cory^2))
    y.btest[id1test]       <- msm::rtnorm(nid1test, mean = mean.impute.b1test, lower = 0)  # Equation (6) of Albert and Chib (1993)
    y.btest[id0test]       <- msm::rtnorm(nid0test, mean = mean.impute.b0test, upper = 0)


    ## update beta ##
    Sigmachol   <- chol(Sigma)
    Sigmalower  <- t(Sigmachol)  # lower traingular matrix
    y           <- c(y.s, y.b)
    ySigma      <- kronecker(diag(n.s), Sigmalower) %*% y
    XbigSigma   <- kronecker(diag(n.s), Sigmalower) %*% Xbig
    lambda.star <- c(tau[1]*lambda, tau[2]*lambda)

    if(algo==1)
    {
      U      <- as.numeric(lambda.star^2)*t(XbigSigma)
      ## step 1 ##
      u      <- stats::rnorm(l2, l0, lambda.star)
      v      <- XbigSigma%*%u + stats::rnorm(n)
      ## step 2 ##
      v_star <- solve((XbigSigma%*%U + I_n),(ySigma - v))
      beta   <- u + U%*%v_star
    } else if(algo==2)
    {
      Qstar <- t(XbigSigma)%*%XbigSigma
      L     <- chol((Qstar + diag(1/as.numeric(lambda.star^2), p, p)))
      v     <- solve(t(L), t(t(ySigma)%*%XbigSigma))
      mu    <- solve(L, v)
      u     <- solve(L, stats::rnorm(p))
      beta  <- mu + u
    }


    ## update lambda in a block using slice sampling ##
    beta.s <- beta[1:p.s]        # survival coefficients
    beta.b <- beta[(p.s + 1):p]  # binary coefficients
    Beta   <- cbind(beta.s/tau[1], beta.b/tau[2])
    eta    <- 1/lambda^2
    upsi   <- stats::runif(p.s, 0, 1/(1+eta))
    tempps <- apply(Beta^2, 1, sum)/2
    ub     <- (1-upsi)/upsi
    Fub    <- stats::pgamma(ub, (2 + 1)/2, scale = 1/tempps)
    Fub[Fub < (1e-4)] <- 1e-4;  # for numerical stability
    up     <- stats::runif(p.s, 0, Fub)
    eta    <- stats::qgamma(up, (2 + 1)/2, scale = 1/tempps)
    lambda <- 1/sqrt(eta);


    ## update tau ##
    Beta  <- cbind(beta.s/lambda, beta.b/lambda)
    tempt <- apply(Beta^2, 2, sum)/2
    et    <- 1/tau^2
    utau  <- stats::runif(2, 0,1 /(1 + et))
    ubt   <- (1-utau)/utau
    Fubt  <- stats::pgamma(ubt,(p.s + 1)/2,scale = 1/tempt)
    Fubt[Fubt < (1e-8)] = 1e-8;  # for numerical stability
    ut    <- stats::runif(2, 0,Fubt)
    et    <- stats::qgamma(ut, (p.s + 1)/2,scale = 1/tempt)
    tau   <- 1/sqrt(et)


    ## Update Sigma using Metropolis Hastings scheme
    Sigma.current <- Sigma
    theta.current <- c(log(Sigma.current[1, 1]), log(Sigma.current[1, 2]))
    trace[i, ]    <- theta.current

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

    theta.proposed <- tmvtnorm::rtmvnorm(1, mean = theta.current, sigma = prop.sigma,
                               lower = c(2 * theta.current[2], -Inf))
    Sigma.proposed <- matrix(c(exp(theta.proposed[1]), exp(theta.proposed[2]), exp(theta.proposed[2]), 1), nrow = 2)
    y.s.tilde      <- y.s - X %*% beta.s
    y.b.tilde      <- y.b - X %*% beta.b
    kernel         <- min(1, exp(sum(mvtnorm::dmvnorm(cbind(y.s.tilde, y.b.tilde), sigma = Sigma.proposed, log = TRUE)) -
                                   sum(mvtnorm::dmvnorm(cbind(y.s.tilde, y.b.tilde), sigma = Sigma.current, log = TRUE)) +
                                   mvtnorm::dmvnorm(theta.proposed, log = TRUE) -
                                   mvtnorm::dmvnorm(theta.current, log = TRUE) +
                                   abs(sum(theta.proposed)) - abs(sum(theta.current))))
    if(stats::rbinom(1, size = 1, kernel) == 1)
    {
      theta.current <- theta.proposed
      Sigma.current <- Sigma.proposed
    }
    Sigma <- Sigma.current


    ## Prediction ##
    # mean.stest          <- Xtest %*% beta.s + cory * sqrt(Sigma[1, 1]/Sigma[2, 2]) * (y.btest - Xtest %*% beta.b)
    mean.stest          <- Xtest %*% beta.s + stats::cor(y.s, y.b) * sqrt(Sigma[1, 1]/Sigma[2, 2]) * (y.btest - Xtest %*% beta.b)
    sd.stest            <- sqrt(Sigma[1, 1] * (1 - cory^2))
    # mean.stest          <- Xtest %*% beta.s
    # sd.stest            <- sqrt(Sigma[1, 1])
    predictive.survivor <- stats::pnorm((y.stest - mean.stest)/sd.stest, lower.tail = FALSE)
    logt                <- mean.stest

    ## Following is required for DIC computation
    loglikelihood.s <- sum(c(stats::dnorm(y.s.observed, mean = X.observed %*% beta.s, sd = sqrt(Sigma[1, 1]),
                                   log = TRUE),
                             log(1 - stats::pnorm(y.s.censored, mean = X.censored %*% beta.s,
                                           sd = sqrt(Sigma[1, 1])))))
    loglikelihood.b <- sum(stats::dbinom(z, size = 1, prob = stats::pnorm(X %*% beta.b), log = TRUE))



    if (i%%500 == 0)
    {
      print(i)
    }

    if(i > burn && i%%thin== 0)
    {
      betaout[, (i-burn)/thin]            <- beta
      lambdaout[ ,(i-burn)/thin]          <- lambda
      tauout[, (i-burn)/thin]             <- tau
      Sigmaout[, ,(i-burn)/thin]          <- Sigma
      predsurvout[ ,(i - burn)/thin]      <- predictive.survivor
      logtimeout[, (i - burn)/thin]       <- logt
      loglikelihood.sout[(i - burn)/thin] <- loglikelihood.s
      loglikelihood.bout[(i - burn)/thin] <- loglikelihood.b
    }
  }


  pMean            <- apply(betaout, 1, mean)
  pMedian          <- apply(betaout, 1, stats::median)
  pLambda          <- apply(lambdaout, 1, mean)
  pSigma           <- apply(Sigmaout, c(1, 2), mean)
  pTau             <- apply(tauout, 1, mean)
  pPS              <- apply(predsurvout, 1, mean)
  pLogtime         <- apply(logtimeout, 1, mean)
  pLoglikelihood.s <- mean(loglikelihood.sout)
  pLoglikelihood.b <- mean(loglikelihood.bout)



  ## Compute DIC ##
  pMean.s <- pMean[1:p.s]
  pMean.b <- pMean[(p.s + 1):p]

  loglikelihood.posterior.s <- sum(c(stats::dnorm(y.s.observed, mean = X.observed %*% pMean.s, sd = sqrt(pSigma[1, 1]),
                                           log = TRUE),
                                     log(1 - stats::pnorm(y.s.censored, mean = X.censored %*% pMean.s,
                                                   sd = sqrt(pSigma[1, 1])))))
  loglikelihood.posterior.b <- sum(stats::dbinom(z, size = 1, prob = stats::pnorm(X %*% pMean.b), log = TRUE))

  DIC.s <- -4 * pLoglikelihood.s + 2 * loglikelihood.posterior.s
  DIC.b <- -4 * pLoglikelihood.b + 2 * loglikelihood.posterior.b


  #construct credible sets
  left <- floor(alpha*effsamp/2)
  right <- ceiling((1-alpha/2)*effsamp)

  BetaSort <- apply(betaout, 1, sort, decreasing = F)
  left.points <- BetaSort[left, ]
  right.points <- BetaSort[right, ]



  result=list("Beta.sHat" = pMean.s, "Beta.bHat" = pMean.b,
              "LeftCI.s" = left.points[1:p.s], "RightCI.s" = right.points[1:p.s],
              "LeftCI.b" = left.points[(p.s + 1):p], "RightCI.b" = right.points[(p.s + 1):p],
              "Beta.sMedian" = pMedian[1:p.s], "Beta.bMedian" = pMedian[(p.s + 1):p],
              "SigmaHat" = pSigma, "LambdaHat" = pLambda, "TauHat" = pTau,
              "Beta.sSamples" = betaout[1:p.s, ], "Beta.bSamples" = betaout[(p.s + 1):p, ],
              "LambdaSamples" = lambdaout, "TauSamples" = tauout, "SigmaSamples" = Sigmaout,
              "DIC.s" = DIC.s, "DIC.b" = DIC.b,
              "SurvivalHat" = pPS, "LogTimeHat" = pLogtime)
  return(result)
}

