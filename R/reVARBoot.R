#' Resampling for variance estimation using reVAR
#'
#' This function gives the estimates using bootstrapped samples.
#'
#' @param nonEventDF1 A data frame that is used for fitting the non-event time model, contains ID, event time information and covariates.
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
#' @param tau Numeric value: data on [0,\code{tau}] will be used
#' in the estimating equation in the second step.
#'
#' @param h Numeric value: bandwidth for kernel estimation in the estimating equations.
#' @param baseline Logical value: if TRUE, the cumulative baseline rate functions of the event and non-event processes are returned.
#'
#' @param h2 Numeric value: bandwidth for kernel estimation for the baseline rate function in the non-event visit model. Only needed if \code{baseline = TRUE}.
#' @param n0 Integer value: baseline cumulative rate functions at \code{n0} equally spaced points on [0,\code{tau}] are returned. Only needed if \code{baseline = TRUE}.
#'
#' @param B Integer value: number of bootstrap samples.
#'
#' @return A list containing the following components:
#'
#' \item{beta}{A matrix of the coefficients in the event model from the bootstrap samples;
#' each row corresponds to the estimate on a bootstrapp sample}
#'
#' \item{alpha}{A matrix of the coefficients in the non-event visit model from the bootstrap samples;
#' each row corresponds to the estimate on a bootstrapp sample}
#'
#' If \code{baseline = TRUE}, the function also returns the following components:
#'
#' \item{t0}{A vector of time points on which the baseline functions are evaluated}
#'
#' \item{M0}{A matrix of the cumulative baseline rate function in the event model;
#' each row corresponds to the estimate on a bootstrapp sample}
#'
#' \item{L0}{A matrix of the cumulative baseline rate function in the non-event visit model;
#' each row corresponds to the estimate on a bootstrapp sample}
#'
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
#' bt <- reVARBoot(ndf1, ndf2, edf, tau = 4.5, h = 0.4, B = 500)
#' # Variance Estimation for coefficients in the event model
#' apply(bt$beta, 2, sd)



reVARBoot <- function(nonEventDF1, nonEventDF2, eventDF,
                      tau, h,
                      baseline = FALSE, h2 = h, n0 = 100,
                      B = 500)
{
  est <- function(nonEventDF1, u.all, t.all, zwu.all, zwt.all, tau, baseline, n0)
  {
    xt2.all <- as.matrix(nonEventDF1[nonEventDF1$Status==1,-c(1:4)])
    zwt2.all <- as.matrix(zwt.all[t.all>0,])

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
        s0mat[j, ] <- K((uj-t2.all)/h)*expXa
      }else{
        s0mat[j, ] <- K((h-t2.all)/h)*expXa
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
          l1u2.all[i] <- sum(K((h-t1)/h2)*dL1)
        }else{
          l1u2.all[i] <- sum(K((u2.all[i]-t1)/h2)*dL1)
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
      return( list(beta = beta.est,
                   alpha = alpha.est) )
    }

  }

  K <- function(x)
  {
    (abs(x)<1)*0.75*(1-(x)^2)
  }

  mu <- table(eventDF$ID)
  # mu-1 is the number of event visits for N subjects
  mt <- table(nonEventDF2$ID)
  # mt-1 is the number of non-event visits for N subjects
  mdat <- table(nonEventDF1$ID)
  # mdat-1 is the number of visits for N subjects
  N <- length(unique(eventDF$ID))

  P <- dim(eventDF)[2]-2
  P2 <- dim(nonEventDF1)[2]-4

  res1 <- matrix(ncol = P, nrow = B)
  res2 <- matrix(ncol = P2, nrow = B)

  if(baseline == TRUE)
  {
    res3 <- matrix(ncol = n0, nrow = B)
    res4 <- matrix(ncol = n0, nrow = B)
  }

  num.all1L <- split(1:length(eventDF$ID), eventDF$ID)
  num.all2L <- split(1:length(nonEventDF2$ID), nonEventDF2$ID)
  num.all3L <- split(1:length(nonEventDF1$ID), nonEventDF1$ID)

  for(b in 1:B)
  {
    bind <- sample(1:N, N, replace = TRUE)

    bind1 <- unlist(num.all1L[bind])
    bind2 <- unlist(num.all2L[bind])
    bind3 <- unlist(num.all3L[bind])

    nonEventDF1b <- nonEventDF1[bind3,]
    nonEventDF2b <- nonEventDF2[bind2,]
    eventDFb <- eventDF[bind1,]

    ut.allb <- nonEventDF1b$Start
    ut.all2b <- nonEventDF1b$End
    ind.allb <- nonEventDF1b$Status
    ind2.allb <- nonEventDF1b$Status2
    xw.allb <- nonEventDF1b[, -c(1:4)]


    id.all1b <- rep(1:N, mu[bind])
    id.all2b <- rep(1:N, mt[bind])
    id.all3b <- rep(1:N, mdat[bind])


    u.allb <- eventDFb$Time
    t.allb <- nonEventDF2b$Time
    zwu.allb <- eventDFb[,-c(1:2)]
    zwt.allb <- nonEventDF2b[,-c(1:2)]


    datb <- data.frame(ID = id.all3b, Start = ut.allb,
                       End = ut.all2b, Status = ind.allb)
    datb <- cbind(datb, xw.allb)


    fitb <- est(datb, u.allb, t.allb, zwu.allb, zwt.allb, tau, baseline, n0)

    res1[b,] <- fitb$beta
    res2[b,] <- fitb$alpha

    if(baseline == TRUE)
    {
      res3[b,] <- fitb$M0
      res4[b,] <- fitb$L0
    }
  }

  if(baseline == TRUE)
  {
    return(list(beta = res1,
                alpha = res2,
                t0 = fitb$t0,
                M0 = res3,
                L0 = res4))
  }else{
    return(list(beta = res1,
                alpha = res2))
  }

}


