#' @export
#'
#' @title coef.eoar - Coefficient extractor for EoAR models.
#'
#' @description Extracts the log-linear model coefficients from an
#' \code{eoar} model object.  See \code{\link{eoar}}.
#'
#' @param obj An object of class \code{eoar}.
#'
#' @return Estimates of the log-linear model's coefficients.
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{eoar}}, \code{\link{labels.eoar}}
#'
#' @examples
#' # A 3 year study of 7 sites. 21 "cells". lambda change = 20/year
#' set.seed(9430834) # fixes Y and g of this example, but not the RNG's used in chains
#' ns <- 3
#'
#' ny <- 7
#' g <- data.frame(
#'  alpha = rnorm(ns*ny,70,2),
#'  beta = rnorm(ns*ny,700,25)
#' )
#' Y <- rbinom(ns*ny, c(rep(20,ny), rep(40,ny), rep(60,ny)), g$alpha/(g$alpha+g$beta))
#'
#' df <- data.frame(year=factor(c(rep("2015",ny),rep("2016",ny),rep("2017",ny))),
#'    Year=c(rep(1,ny),rep(2,ny),rep(3,ny)))
#'
#' # Uninformed eoar (use low number of iterations because it's and example)
#' eoar.1 <- eoar(Y~year, g, df, nburn = 1000, niters= 50*10, nthins = 10 )
#'
#' coef(eoar.1)
coef.eoar <- function(obj){
  obj$estimates[obj$coef.labels,"Estimate"]
}
