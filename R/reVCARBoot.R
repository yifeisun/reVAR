#' Resampling for variance estimation using reVCAR
#'
#' This function gives the estimates using bootstrapped samples.
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
#' @param tau Numeric value: data on [0,\code{tau}] will be used in the estimating equation in the second step.
#' @param h Numeric value: bandwidth for kernel estimation.
#' @param B Integer value: number of bootstrap samples.
#'
#' @return A matrix of the coefficients from the bootstrap samples;
#' each row corresponds to the estimate on a bootstrapp sample
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)
#' ndf2 <- dat$nonEventDF2
#' edf <- dat$eventDF
#' fit <- reVCAR(ndf2, edf, tau = 4.5, h = 0.4)
#' bt <- reVCARBoot(ndf2, edf, tau = 4.5, h = 0.4, B = 500)
#' # Variance Estimation
#' apply(bt, 2, sd)


reVCARBoot <- function(nonEventDF2, eventDF, tau, h, B)
{
  est0 <- function(u.all, t.all, zwu.all, zwt.all, tau, baseline, n0)
  {
    zwt2.all <- as.matrix(zwt.all[t.all>0,])

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
        s0mat[j, ] <- K((uj-t2.all)/h)
      }else{
        s0mat[j, ] <- K((h-t2.all)/h)
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

    return( list(beta = beta.est) )
  }

  K <- function(x)
  {
    (abs(x)<1)*0.75*(1-(x)^2)
  }

  mu <- table(eventDF$ID)
  # mu-1 is the number of event visits for N subjects
  mt <- table(nonEventDF2$ID)
  # mt-1 is the number of non-event visits for N subjects

  N <- length(unique(eventDF$ID))
  P <- dim(eventDF)[2]-2

  res1 <- matrix(ncol = P, nrow = B)

  num.all1L <- split(1:length(eventDF$ID), eventDF$ID)
  num.all2L <- split(1:length(nonEventDF2$ID), nonEventDF2$ID)

  for(b in 1:B)
  {
    bind <- sample(1:N, N, replace = TRUE)

    bind1 <- unlist(num.all1L[bind])
    bind2 <- unlist(num.all2L[bind])

    nonEventDF2b <- nonEventDF2[bind2,]
    eventDFb <- eventDF[bind1,]

    id.all1b <- rep(1:N, mu[bind])
    id.all2b <- rep(1:N, mt[bind])

    u.allb <- eventDFb$Time
    t.allb <- nonEventDF2b$Time
    zwu.allb <- eventDFb[,-c(1:2)]
    zwt.allb <- nonEventDF2b[,-c(1:2)]

    fitb <- est0(u.allb, t.allb, zwu.allb, zwt.allb, tau, baseline, n0)

    res1[b,] <- fitb$beta
  }

  res1

}
