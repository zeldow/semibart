#'Semiparametric regression using BART. For continuous outcomes \eqn{y}, the model is
#'\eqn{y = \omega(x) + a \beta + \epsilon}, where \eqn{\epsilon \sim N(0,\sigma^2)}, 
#'\eqn{x} are some covariates, and \eqn{a} is a smaller subset of covariates (and possibly
#'interactions) that may be of immediate scientific interest. The 
#'covariates \eqn{a} represent the design matrix for variables and interactions that are 
#'modeled parametrically. The functional form of \eqn{\omega(x)} is unspecified and modeled
#'using Bayesian Additive Regression Trees (BART) (Chipman et al, 2010). To complete the 
#'model, we use a normal prior on \eqn{\beta} and an inverse chi square prior on 
#'\eqn{\sigma}. 
#'
#'For binary \eqn{y}, the model is \eqn{P(Y=1 | x, a) = F(\omega(x) + a\beta)}, where \eqn{F} 
#'denotes the standard normal cdf (probit link) and BART is used to model the nonparametric 
#'\eqn{\omega(x)}.
#'
#'The covariates in the parametric and nonparametric parts may overlap. That is, a covariate
#'included in \eqn{x} may also be included in the parametric part as an interaction with the
#'exposure variable if effect modification is of scientific interest. In this case, special 
#'care is recommended as a larger sample size or larger number of trees may be needed.
#'
#'For information on causal interpretations as structural mean model, see 
#'Vansteelandt et al (2014) and Zeldow et al (2016).
#'
#'
#'@references Chipman, H., George, E., and McCulloch R. (2010) 
#'  Bayesian Additive Regression Trees.
#'  \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298.
#'@references Vansteelandt, S, and Joffe, M. (2014) Structural nested models and g-estimation: 
#'  The partially realized promise. \emph{Statistical Science:} 707-731.
#'@references Zeldow, B, Lo Re, V, Roy, J. (2016) Bayesian semiparametric 
#'  regression and structural mean models with BART.
#'
#'
#'@param x.train Design matrix of values to be modeled with BART.
#'@param Design matrix of values to be modeled linearly.
#'@param y.train Vector of outcomes (continuous or binary). When binary, elements
#'  must be either 0 or 1.
#'@param sigest Estimate of regression error. If no value is supplied and sigest=NA,
#'  the least squares estimate is used. Must be a positive number. Ignored if y.train
#'  is binary.
#'@param sigdf Degrees of freedom on prior for error variance.
#'@param sigquant The quantile of the prior that the rough estimate (see sigest) is 
#'  placed at. The closer the quantile is to 1, the more aggresive the fit will be as 
#'  you are putting more prior weight on error standard deviations (\eqn{\sigma}{sigma}) 
#'  less than the rough estimate. Not used if y.train is binary.
#'@param k For numeric y, k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} 
#'  is away from +/-.5. The response (y.train) is internally scaled to range from -0.5 to 0.5. 
#'  For binary y, k is the number of prior standard deviations \eqn{f(x)} is away from +/-3. 
#'  In both cases, the bigger k is, the more conservative the fitting will be.
#'@param power Power parameter for prior on tree depth.
#'@param base Base parameter on prior on tree depth.
#'@param meanb Prior mean on regression coefficients. Length must equal # columns in a.train, 
#'  that is: length(meanb) == ncol(a.train).
#'@param sigb Prior standard deviation on regression coefficients. Prior is 
#'  \eqn{\beta \sim N(meanb, sigb^2 I)} where \eqn{I} is the identity matrix of 
#'  appropriate dimension.
#'@param ntree Number of trees to use for BART.
#'@param ndpost Number of MCMC iterations, including burn-in.
#'@param numcut Number of cutpoints for each variable in BART. Must be of length 1 or have 
#'  length ncol(x.train).
#'@param usequants Indicates whether to use observed quantiles for cutpoints or evenly 
#'  spaced cutpoints based on min and max for each column in x.train.
#'@param offset Offset for regression -- used only when outcome is binary.
#'@param binarylink Indicates whether to use probit or logit link for binary data. 
#'  Currently only the probit link is supported.
#'@param verbose Indicates whether or not user wants printed output to check progress of 
#'  MCMC algorithm.
#'@param printevery Indicates how often to print an update on completion of algorithm. 
#'  Default is to print a message every 100 iterations. Ignored if verbose = FALSE.
#'
#'
#'@return Returns a list containing a matrix of MCMC draws for regression parameters 
#'  (the dimension is ndpost x ncol(a.train)). When y.train is continuous also returns 
#'  vector of draws of the error variance. Retrieve the regression parameters using $beta 
#'  and $sigma, for the regression parameters and variance parameters, respectively.
#'
#'
#'@importFrom stats lm
#'
#'@export
#'
#'
#'@examples
#'set.seed(1)
#'n <- 200; nc <- 5
#'x <- matrix(rnorm(n * nc), nrow = n, ncol = nc)
#'a <- rbinom(n, 1, 0.5)
#'y <- 2 + 3 * x[ ,1] + 0.5 * x[ ,2] - 2 * x[ ,3] + 5 * x[ ,5] + 2 * a + rnorm(n)
#'sb <- semibart(x, as.matrix(a), y)
semibart = function(
  x.train,a.train, y.train,
  sigest=NA, sigdf=3, sigquant=.90, 
  k=2.0,
  power=2.0, base=.95,
  meanb=rep(0,ncol(a.train)), sigb=4,
  ntree=200,
  ndpost=1000,
  numcut = 100,
  usequants = FALSE,
  offset = 0,
  binarylink = "probit",
  verbose = TRUE,
  printevery = 100)
{
  binary=FALSE

  if(is.factor(y.train)) {
    if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
    binary = TRUE
    y.train = as.numeric(y.train)-1
  } else {
    if((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0)) {
      cat('NOTE: assumming numeric response is binary\n')
      binary = TRUE
      if ( binarylink != "probit" & binarylink != "logit" & binarylink != "logistic" ) stop("for binary outcomes, must choose probit or logit link (default = probit)")
    }
  }
  
  if(mode(numcut)!="numeric") stop("input numcut must be a numeric vector")
  if(length(numcut)==1) numcut = rep(numcut,ncol(x.train))
  if(length(numcut) != ncol(x.train)) stop("length of numcut must equal number of columns of x.train")
  numcut = as.integer(numcut)
  if(min(numcut)<1) stop("numcut must be >= 1")
  
  if((mode(k)!="numeric") || (k<0)) stop("input k must be a positive number")
  
  if((mode(sigquant)!="numeric") || (sigquant<0)) stop("input sigquant must be a positive number")
  if((mode(ntree)!="numeric") || (ntree<0)) stop("input ntree must be a positive number")
  if((mode(ndpost)!="numeric") || (ndpost<0)) stop("input ndpost must be a positive number")
  
  if(typeof(usequants) != "logical") stop("input usequants must be logical")  
  if(typeof(verbose) != "logical") stop("input verbose must be logical")

  probitlink <- 0;
  if (binarylink == "logit" || binarylink == "logistic") probitlink <- 1
  rgy = range(y.train)
  y = -.5 + (y.train-rgy[1])/(rgy[2]-rgy[1])

  # if sigest=NA, fit a lm to training data to get the value of sigest...
  # sigest is on the scale of the transformed y, so we do the lm after the scaling above...
  if(!binary) {
    if (is.na(sigest)) {
      templm = stats::lm(y~x.train+a.train-1) 
      sigest = summary(templm)$sigma
    } else {
      sigest = sigest/(rgy[2]-rgy[1]) #put input sigma estimate on transformed scale
    }
  } else {
    sigest=1
  }
  
  if(!binary) {
    offset = -1000.0
  }
  
  bartres_ <- semibart_cpp(x.train,a.train,y,as.double(sigest),
                    as.integer(sigdf),as.double(sigquant),as.double(k),
                    as.double(power),as.double(base),
                    as.double(meanb),as.double(sigb),
                    as.integer(ntree),as.integer(ndpost),
                    as.integer(numcut),as.integer(usequants),
                    as.double(offset), as.integer(probitlink),
                    as.integer(verbose),as.integer(printevery))
  
  # now read in the results...
  if(!binary) {
    bartres_$sigmaReps = bartres_$sigmaReps*(rgy[2]-rgy[1])
    bartres_$betaReps = bartres_$betaReps*(rgy[2]-rgy[1])
  }

         
  if(binary) {
    retval = list(
      beta = bartres_$betaReps
      )
  } else {
    retval = list(
      sigma = bartres_$sigmaReps,
      beta = bartres_$betaReps
      )
  }
  class(retval) = 'semibart'
  return(invisible(retval))
}
