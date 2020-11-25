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
#' @return A list of posterior predictive distributions.
#' @export
#' @author Andrew Tredennick
PostPredictEoAR <- function(object, 
                            newdata,
                            newoffset, 
                            beta.params = NULL, 
                            nyears = 1) {
  X <- model.matrix(object, newdata)
  betaMat <- t(as.matrix(object$out[,object$coef.labels]))
  
  n <- nrow(X)
  k <- ncol(betaMat)
  
  # loop over observation rows
  lambdaExp <- matrix(data = 0, nrow = n, ncol = k)
  M <- matrix(data = 0, nrow = n, ncol = k)
  C <- matrix(data = 0, nrow = n, ncol = k)
  for(i in 1:nrow(X)) {
    lambda <- X[i, ] %*% betaMat + log(newoffset[i])
    lambdaExp[i, ] <- exp(lambda)
    lambdaScl <- lambdaExp[i, ] * nyears
    M[i, ] <- rpois(k, lambdaScl)
    
    if(!is.null(beta.params)){
      g <- rbeta(1, beta.params$alpha[i], beta.params$beta[i])
      C[i, ] <- rbinom(k, M[i, ], g)
    }
  }
  
  if(is.null(beta.params)){
    C <- NULL
  }
  
  out <- list(
    postLambda = lambdaExp,
    postMortality = M,
    postCarcasses = C
  )
  return(out)
}