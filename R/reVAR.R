#' Fit the proportional rate model with intermittently measured time-dependent covariates under VAR
#'
#' This function fits the proportional rate model with intermittently
#' measured time-dependent covariates under the visiting at random (VAR) assumption.
#'
#'#' @param nonEventDF1 A data frame that is used for fitting the non-event time model, contains ID, event time information and covariates.
#' Each row records the covariates on a time interval.
#' The first four columns must be \code{ID}, \code{Start}, \code{End}, \code{Status}:
#' \code{ID} is the subject ID for each observation;
#' \code{Start} is the starting time of the interval;
#' \code{End} is the ending time of the interval;
#' \code{Status} is the status indicator, 1 = non-event visit, 0 = no.
#' The other columns are covariate values for the non-event model on the time interval (\code{Start}, \code{End}].
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
#' @param h Numeric value: bandwidth for kernel estimation in the estimating equations.
#' @param baseline Logical value: if TRUE, the cumulative baseline rate functions of the event and non-event processes are returned.
#'
#' @param h2 Numeric value: bandwidth for kernel estimation for the baseline rate function in the non-event visit model. Only needed if \code{baseline = TRUE}.
#' @param n0 Integer value: baseline cumulative rate functions at \code{n0} equally spaced points on [0,\code{tau}] are returned. Only needed if \code{baseline = TRUE}.
#'
#' @return A list containing the following components:
#'
#' \item{beta}{A vector of coefficients in the event model}
#'
#' \item{alpha}{A vector of coefficients in the non-event visit model}
#'
#' If \code{baseline = TRUE}, the function also returns the following components:
#'
#' \item{t0}{A vector of time points on which the baseline functions are evaluated}
#'
#' \item{M0}{A vector of the cumulative baseline rate function in the event model}
#'
#' \item{L0}{A vector of the cumulative baseline rate function in the non-event visit model}
#'
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)
#' ndf1 <- subset(dat$nonEventDF1_s2, select = -Status2)
#' ndf2 <- dat$nonEventDF2
#' edf <- dat$eventDF
#' fit <- reVAR(ndf1, ndf2, edf, tau = 4.5, h = 0.4)



reVAR <- function(nonEventDF1, nonEventDF2, eventDF,
                  tau, h,
                  baseline = FALSE, h2 = h, n0 = 100)
{
  u.all <- eventDF$Time
  t.all <- nonEventDF2$Time

  zwu.all <- as.matrix(eventDF[-c(1:2)])
  xt2.all <- as.matrix(nonEventDF1[nonEventDF1$Status==1,-c(1:4)])
  zwt2.all <- as.matrix(nonEventDF2[t.all>0,-c(1:2)])

  P <- dim(zwu.all)[2]
  # First step
  foright <- paste(names(nonEventDF1)[-c(1:4)],collapse = "+")
  fo <- paste("survival::Surv(Start, End, Status) ~ ",foright, sep = "")

  fit1 <- survival::coxph(as.formula(fo), data = nonEventDF1)
  alpha.est <- as.vector(coef(fit1))

  # Second step

  # With 2, there is no visit at time 0
  u2.all <- u.all[u.all>0]
  t2.all <- nonEventDF1$End[nonEventDF1$Status==1] # same as t.all[t.all>0]

  expXa <- exp(-xt2.all%*%(alpha.est))

  s1arr <- array(dim = c(P, length(u2.all), length(t2.all)))#matrix(nrow = length(u2.all), ncol = length(t2.all))
  s0mat <-  matrix(nrow = length(u2.all), ncol = length(t2.all))

  for(j in 1:length(u2.all))
  {
    uj <- u2.all[j]
    if(uj>h)
    {
      s0mat[j, ] <- K((uj-t2.all)/h)*expXa/h
    }else{
      s0mat[j, ] <- K((h-t2.all)/h)*expXa/h
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

  if(baseline == TRUE)
  {
    hz1 <- withCallingHandlers(survival::basehaz(fit1, centered = FALSE),
                               warning=function(w){invokeRestart("muffleWarning")})
    L1 <- c(0, hz1$hazard)
    dL1 <- diff(L1)
    t1 <- hz1$time

    l1u2.all <- u2.all
    for(i in 1:length(u2.all))
    {
      if(u2.all[i]<h)
      {
        l1u2.all[i] <- sum(K((h-t1)/h2)*dL1)/h2
      }else{
        l1u2.all[i] <- sum(K((u2.all[i]-t1)/h2)*dL1)/h2
      }
    }

    expZbeta <- as.vector(exp(zwt2.all%*%beta.est))
    s0u2.all <- colSums(t(s0mat)*expZbeta)

    mu2.all <- l1u2.all/s0u2.all
    mu2.all[is.na(mu2.all)] <- 0
    dM0u2.all <- mu2.all[order(u2.all)]

    t0 <- seq(from = 0, to = tau, length = n0)
    Lt <- function(tt)
    {
      sum(dL1[t1<=tt])
    }
    Mt <- function(tt)
    {
      sum(dM0u2.all[u2.all<=tt])
    }

    L0 <- sapply(t0, Lt)
    M0 <- sapply(t0, Mt)

    return( list(beta = beta.est,
                 alpha = alpha.est,
                 t0 = t0,
                 M0 = M0,
                 L0 = L0) )
  }else{
    return( list(beta = beta.est, alpha = alpha.est) )
  }
}



