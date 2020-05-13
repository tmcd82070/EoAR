#' @name print.eoa
#'
#' @title print.eoar - Print a Evidence of Absence Regression model
#'
#' @description Print method for EoAR models produced by \code{eoar},
#' which are of class \code{eoar}.
#'
#' @param x An estimated eoar object from \code{eoar}.
#'
#' @param \dots Included for compatibility with other print methods.  Ignored here.
#'
#' @return The input value of \code{obj} is invisibly returned.
#' @author Trent McDonald, WEST Inc. \email{tmcdonald@west-inc.com}
#'
#' @seealso \code{\link{summary.eoar}}
#' @examples
#'
#' @keywords models
#' @export

print.eoar <- function( x, ... ){

  summary(x)

}
