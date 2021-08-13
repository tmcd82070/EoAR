#' Make posterior predictions
#' 
#' @param object A fitted \code{eoar} object.
#' @param newdata A data frame of covariates for new observations.
#' @param newoffset A vector offsets for new observations.
#' @param beta.params A data frame containing, at a minimum, two columns 
#'     named \code{$alpha} and \code{$beta}. These are the alpha and beta 
#'     parameters of a Beta distribution that determine the g-values in 
#'     each "cell". Length of \code{$alpha} and \code{$beta} is either one, 
#'     in which case g is assumed constant over "cells", or one element 
#'     per "cell".
#' @param nyears A scalar for the number of years to make predictions for.
#' @param meanLambda A logical indicating whether the sample lambda according
#'     to the estimated lambda standard deviation (FALSE, defualt) or to
#'     ignore sampling of lambda to predict using the mean expectation (TRUE).  
#' @return A list of posterior predictive distributions.
#' @export
#' @author Andrew Tredennick
PostPredictEoAR <- function(object, 
                            newdata,
                            newoffset, 
                            beta.params = NULL, 
                            nyears = 1,
                            meanLambda = FALSE) {
  X <- model.matrix(object, newdata)
  betaMat <- t(as.matrix(object$out[,object$coef.labels]))
  sigmaVec <- as.numeric(as.matrix(object$out[,"sigmaLambda"]))
  if(meanLambda == TRUE) {
    sigmaVec[] <- 0
  }
  
  n <- nrow(X)
  k <- ncol(betaMat)
  
  # loop over observation rows
  lambdaExp <- matrix(data = 0, nrow = n, ncol = k)
  M <- matrix(data = 0, nrow = n, ncol = k)
  MOneYear <- matrix(data = 0, nrow = n, ncol = k)
  C <- matrix(data = 0, nrow = n, ncol = k)
  for(i in 1:nrow(X)) {
    lambdaMu <- X[i, ] %*% betaMat
    lambdaLog <- rnorm(k, lambdaMu[1, ], sigmaVec)
    lambdaExp[i, ] <- exp(lambdaLog)
    lambdaScl <- lambdaExp[i, ] * nyears
    M[i, ] <- rpois(k, lambdaScl * newoffset[i])
    MOneYear[i, ] <- rpois(k, lambdaExp[i, ] * newoffset[i])
    
    if(!is.null(beta.params)){
      g <- rbeta(1, beta.params$alpha[i], beta.params$beta[i])
      C[i, ] <- rbinom(k, M[i, ], g)
    }
  }
  
  if(!is.null(beta.params)){
    Csummary <- t(apply(C, 1, quantile, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  }
  
  if(is.null(beta.params)){
    C <- NULL
    Csummary <- NULL
  }
  
  # Compute fleet-wide mortality
  MfleetOneYear <- apply(MOneYear, 2, sum)
  Mfleet <- apply(M, 2, sum)
  
  Lsummary <- t(apply(lambdaExp, 1, quantile, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  Msummary <- t(apply(M, 1, quantile, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  MOnesummary <- t(apply(MOneYear, 1, quantile, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  Mfleetsummary <- quantile(Mfleet, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
  MfleetOneYearsummary <- quantile(MfleetOneYear, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
  
  
  out <- list(
    Posteriors = list(
      postLambda = lambdaExp,
      postMortality = M,
      postMortalityOneYear = MOneYear,
      postMFleet = Mfleet,
      postMFleetOneYear = MfleetOneYear,
      postCarcasses = C
    ),
    Summaries = list(
      Lambda = Lsummary,
      Mortality = Msummary,
      MortalityOneYear = MOnesummary,
      MFleet = Mfleetsummary,
      MfleetOneYear = MfleetOneYearsummary,
      Carcasses = Csummary
    )
  )
  return(out)
}