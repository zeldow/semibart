#'Semiparametric regression using BART
#'
#'@param x.train Design matrix of values to be modeled with BART
#'@param a.train Design matrix of values to be modeled linearly
#'@param y.train Vector of outcome (continuous or binary)
#'@param sigest Estimate of regression error standard deviation
#'@param sigdf Prior value on error term
#'@param sigquant Prior value on error term
#'@param k Controls how aggressive the fit is
#'@param power Prior on tree depth
#'@param base Prior on tree depth
#'@param meanb Prior mean on regression coefficients -- length must equal # columns in a.train
#'@param sigb Prior standard deviation on regression coefficients ~ sigb^2 I
#'@param ntree Number of trees to use for BART
#'@param ndpost Number of MCMC iterations
#'@param numcut Number of cutpoints for each variable in BART. Must be of length 1 or have length ncol(x.train)
#'@param usequants Indicates whether to use observed quantiles for cutpoints or evenly spaced cutpoints based on min and max
#'@param offset Offset for regression -- used only when outcome is binary
#'@param binarylink Indicates whether to use probit or logit link for binary data
#'@param verbose Indicates whether or not user wants printed output
#'@param printevery Printevery
#'@return Returns matrix of regression parameters (ndpost x ncol(a.train)). When y.train is continuous also returns vector of draws of regression parameter.
semibart = function(
  x.train,a.train, y.train,
  sigest=NA, sigdf=3, sigquant=.90, 
  k=2.0,
  power=2.0, base=.95,
  meanb=rep(0,ncol(a.train)), sigb=4,
  ntree=200,
  ndpost=1000,
  numcut = 100,
  usequants = 1,
  offset = 0,
  binarylink = "probit",
  verbose = 1,
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
  
  probitlink <- 0;
  if (binarylink == "logit" || binarylink == "logistic") probitlink <- 1
  rgy = range(y.train)
  y = -.5 + (y.train-rgy[1])/(rgy[2]-rgy[1])

  # if sigest=NA, fit a lm to training data to get the value of sigest...
  # sigest is on the scale of the transformed y, so we do the lm after the scaling above...
  if(!binary) {
    if (is.na(sigest)) {
      templm = lm(y~x.train+a.train-1) 
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
