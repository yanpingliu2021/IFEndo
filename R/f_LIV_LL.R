#' @title Latent Instrumental Variables (LIV) Log-likelihood function to optimize
#'
#' @description A helper function to construct the log-likelihood function for simple LIV model
#'
#' @param starting.parameters the 8 parameters ("Intercept", "Endogenous regressor", "pi1", "pi2", "lambda1",
#' "Var(E)", "Cov(EV)", "Var(V)") to estimate as outlined in Ebbes 2004
#' @param data the input data.frame
#' @return The log-likelihood function
#' @export
#' @import mvtnorm

########################
# Frequentist LIV
########################
LIV_LL <- function(starting.parameters, data){

  #As outlined in Ebbes 2004, there are 8 parameters to estimate for simple LIV model: Betas, pi's,lambda, Sigma

  #Get starting beta coefficients (these come from lm() if no starting values are specified)
  b00<-starting.parameters[1]
  b01<-starting.parameters[2]

  # Probability for group membership; this must be bounded between 0 and 1 as log(-x) will produce NaN's
  #Give hard constraints
  lambda<-ifelse(starting.parameters[5] < 0, 0,
                 ifelse(starting.parameters[5] > 1, 1, starting.parameters[5]))
  sigma<-matrix(c(starting.parameters[6], starting.parameters[7],
                  starting.parameters[7],starting.parameters[8]),2)
  #Varcov matrix reduced form from Ebbes 2005
  s2e     <- sigma[1,1]
  s2v     <- sigma[2,2]
  sev     <- sigma[1,2]
  varcov      <- matrix(0,2,2)
  varcov[1,1] <- b01*b01*s2v+2*b01*sev+s2e
  varcov[2,1] <- b01*s2v+sev
  varcov[1,2] <- varcov[2,1]
  varcov[2,2] <- s2v
  #Group 1 mean vector from Ebbes 2004
  pi1 <- starting.parameters[3]
  mu1 <- matrix(data = c(b00+b01*pi1, pi1), nrow=2, ncol=1)
  #Group 2 mean vector from Ebbes 2004
  pi2 <- starting.parameters[4]
  mu2 <- matrix(c(b00+b01*pi2, pi2), nrow=2,ncol=1)
  #Log Likelihood to be maximized (we multiply by -1 and minimize)
  log.pdf1 <- mvtnorm::dmvnorm(data, mean=mu1, sigma=varcov, log = TRUE)
  log.pdf2 <- mvtnorm::dmvnorm(data, mean=mu2, sigma=varcov, log = TRUE)
  max.AB   <- pmax(log(lambda)+log.pdf1, log(1 - lambda)+log.pdf2)
  logLL.lse<- sum(max.AB+log(lambda*exp(log.pdf1-max.AB) + (1-lambda)*exp(log.pdf2-max.AB)))
  return(-1*logLL.lse)
}
