#' Fit the proportional rate model with intermittently measured time-dependent covariates under VCAR
#'
#' This function fit the proportional rate model with intermittently
#' measured time-dependent covariates under the visiting completely at random (VCAR) assumption.
#'
#' @param nonEventDF2 A data frame that contains the ID, non-event visit times and covariates (for the event model) measured at each non-event visit.
#' The first two columns must be \code{ID} and \code{Time}:
#' \code{ID} is the subject ID for each observation;
#' \code{Time} include time zero and the non-event visit times.
#' The other columns are covariate values measured at \code{Time}.
#'
#' @param eventDF A data frame that contains the ID, event visit times and covariates (for the event model) measured at each event visit.
#' The first two columns must be \code{ID} and \code{Time}:
#' \code{ID} is the subject ID for each observation;
#' \code{Time} include time zero and the event times.
#' The other columns are covariate values measured at \code{Time}.
#'
#'
#' @param tau Numeric value: data on [0,\code{tau}] will be used in the estimating equation in the second step.
#' @param h Numeric value: bandwidth for kernel estimation.
#'
#'
#' @return A vector of the coefficients
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)
#' ndf2 <- dat$nonEventDF2
#' edf <- dat$eventDF
#' fit <- reVCAR(ndf2, edf, tau = 4.5, h = 0.4)


reVCAR <- function(nonEventDF2, eventDF, tau, h)
{
  u.all <- eventDF$Time
  t.all <- nonEventDF2$Time

  zwu.all <- as.matrix(eventDF[-c(1:2)])
  zwt2.all <- as.matrix(nonEventDF2[t.all>0,-c(1:2)])

  P <- dim(zwu.all)[2]

  # With 2, there is no visit at time 0
  u2.all <- u.all[u.all>0]
  t2.all <- t.all[t.all>0]

  s1arr <- array(dim = c(P, length(u2.all), length(t2.all)))#matrix(nrow = length(u2.all), ncol = length(t2.all))
  s0mat <-  matrix(nrow = length(u2.all), ncol = length(t2.all))

  for(j in 1:length(u2.all))
  {
    uj <- u2.all[j]
    if(uj>h)
    {
      s0mat[j, ] <- K((uj-t2.all)/h)/h
    }else{
      s0mat[j, ] <- K((h-t2.all)/h)/h
    }
  }

  for(p in 1:P)
  {
    s1arr[p,,] <- t(t(s0mat)*zwt2.all[,p])
  }

  U <- function(b)
  {
    expZb <- as.vector(exp(zwt2.all%*%b))
    fun1 <- function(matp)
    {
      colSums(t(matp)*expZb)/colSums(t(s0mat)*expZb)
    }
    ee <- apply(s1arr,1,fun1)
    ee[is.na(ee)] <- 0
    colSums(zwu.all[u.all<=tau & u.all>0, ])- colSums(ee[u2.all<=tau,])
  }

  beta.est <- nleqslv::nleqslv(rep(0,P),U)$x
  beta.est
}

