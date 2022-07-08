#' @title Simulated Dataset with One Endogenous Continuous Regressor
#' @description A dataset with one endogenous regressor \code{x}, an instrument \code{z}
#'  used to build \code{x}, an intercept and a dependent variable, \code{y}.
#'  The true parameter values for the coefficients are: \code{b0 = 1} for the intercept
#'  and \code{a1 = 2.5} for \code{x}. The correlation between \code{z} and \code{x} is 0.7
#' @name dataLatentIV
#' @usage data("dataLatentIV")
#' @format A data frame with 2500 observations on 3 variables:
#' \describe{
#' \item{\code{y}}{a numeric vector representing the dependent variable.}
#' \item{\code{x}}{a numeric vector representing the endogenous variable.}
#' \item{\code{z}}{a numeric vector used in the construction of the endogenous variable, x.}
#' }
#' @docType data
#' @author Yanping Liu \email{yanping.liu@@iqvia.com}
"dataLatentIV"
