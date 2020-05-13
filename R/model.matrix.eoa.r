#' @title model.matrix.eoa
#'
#' @description Model matrix method for eoar objects produced
#' by the \code{eoar} function.
#'
#' @param object An object of class \code{eoar}. In particular,
#' this object must have a $terms component.
#'
#' @param \dots Included for compatability with other methods
#'
#' @author Trent McDonald
#'
#' @return The design matrix for an \code{eoar} object with all the factor
#' levels and interactions expanded.
#'
#' @export
#'
model.matrix.eoar <- function(object, data, ...){

  Terms <- delete.response(object$terms)
  if(missing(data) || is.null(data)){
    data <- object$data
  }
  m <- model.frame(Terms, data)
  model.matrix.default(Terms, m)

}
