#' @title  Function to run simple Latent Instrument Variable model via MLE
#'
#' @description  Fits linear models with one endogenous regressor and no additional explanatory variables using the latent instrumental variable approach
#'
#' @param formula A symbolic description of the model to be fitted. Of class "formula".
#' @param endogenous.regressor str the name of the endogenous variable in the input data
#' @param data A data.frame containing the data of all parts specified in the formula parameter.
#' @param iters number of iterations for the optimization.
#' @param parameters A named vector containing a set of parameters to use in the first optimization iteration.
#' The names have to correspond exactly to the names of the components specified in the formula parameter.
#' If not provided, a linear model is fitted to derive them.
#'
#' @return A named vector of all coefficients resulting from model fitting.
#'
#' @seealso \code{\link[optimx]{optimx}} for possible elements of parameter \code{optimx.arg}
#'
#' @references   Ebbes, P., Wedel,M., Böckenholt, U., and Steerneman, A. G. M. (2005). 'Solving and Testing for Regressor-Error
#' (in)Dependence When no Instrumental Variables are Available: With New Evidence for the Effect of Education on Income'.
#' Quantitative Marketing and Economics, 3:365--392.
#'
#' Ebbes P., Wedel M., Böckenholt U. (2009). “Frugal IV Alternatives to Identify the Parameter for an Endogenous Regressor.” Journal of Applied Econometrics, 24(3), 446–468.
#'
#' @examples
#' \donttest{
#' data("dataLatentIV")
#'
#' # function call without any initial parameter values
#' livcoef <- LIV(y~x, "x", dataLatentIV, NULL)
#' livcoef
#'
#'
#' # function call with initial parameter values given by the user
#' fit<-kmeans(dataLatentIV$x, 2)
#' starts<-c(1, 2.5, fit$centers[1], fit$centers[2], 0.5, 1, 0.5, 1)
#' livcoef2<-LIV(y~x, "x", dataLatentIV, NULL, starts)
#' livcoef2
#'
#' }
#'
#' @importFrom stats lm kmeans model.frame
#' @import optimx
#' @export


LIV<-function(formula, endogenous.regressor, data, iters = NULL, parameters=c()){
  #Start parameters will be:
  LIVstart<-numeric(8);
  #Default parameters if nothing is given
  if(is.null(parameters)){
    message("Defaults")
    #Regression for starting parameters
    startlm<-lm(formula, data) #intercept and coefficient
    for(i in 1:2){
      LIVstart[i]<-startlm$coefficients[i]
    }
    #LIVstart[3]<-mean(data[,endogenous.regressor]); LIVstart[4]<-mean(data[,endogenous.regressor])+sd(data[,endogenous.regressor])
    fit<-kmeans(data[,endogenous.regressor], 2)
    LIVstart[3]<-fit$centers[1]; LIVstart[4]<-fit$centers[2]
    LIVstart[5]<-0.5; LIVstart[6]<-1; LIVstart[7]<-0.5; LIVstart[8]<-1}
  else { LIVstart<-parameters }

  names(LIVstart)<-c("Intercept", "Endogenous.regressor", "pi1", "pi2", "lambda1", "Var(E)", "Cov(EV)", "Var(V)")
  #Optimize log-likelihood function
  res.optimx <- tryCatch(expr =
                           optimx::optimx(par = LIVstart,
                                          fn  = LIV_LL,
                                          data = model.frame(formula = formula, data = data),
                                          method = "Nelder-Mead",
                                          hessian = FALSE,
                                          itnmax  = iters,
                                          control = list(kkt=FALSE, trace = 0,
                                                         dowarn = FALSE)),
                         error   = function(e){ return(e)})

  if(is(res.optimx, "error"))
    stop("Failed to optimize the log-likelihood function with error \'", res.optimx$message,
         "\'. Please revise your start parameter and data.", call. = FALSE)
  #Coefficients
  LIVcoefficients<-res.optimx[1:8]; names(LIVcoefficients)<-names(LIVstart)
  return(LIVcoefficients)
}
