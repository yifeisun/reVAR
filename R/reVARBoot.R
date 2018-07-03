#' Resampling for variance estimation using reVAR
#'
#' This function gives the estimates using bootstrapped samples.
#'
#' @param nonEventDF1 A data  frame that is used for fitting the non-event time model, contains ID, event time information and covariates
#' @param nonEventDF2 A dataframe that contains the ID, non-event visit times and covariates measured at each non-event visit
#' @param eventDF A data frame that contains the ID, event visit times and covariates measured at each event visit
#' @param tau Data on [0,tau] will be used in the estimating equation in the second step
#' @param h Bandwidth for kernel estimation
#' @param B Number of bootstrap samples
#' @return A matrix of estimates from bootstrapped samples
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)
#' estBoot <- reVARBoot(dat$nonEventDF1,dat$nonEventDF2,dat$eventDF,4.5, h = 2*100^(-1/3), B = 500)
#' # Variance Estimation
#' apply(estBoot, 2, sd)


reVARBoot <- function(nonEventDF1, nonEventDF2, eventDF, tau, h, B)
{
  est <- function(nonEventDF1, u.all, t.all, zwu.all, zwt.all, tau)
  {
    xt2.all <- as.matrix(nonEventDF1[nonEventDF1$Status==1,-c(1:5)])
    zwt2.all <- as.matrix(zwt.all[t.all>0,])

    P <- dim(zwu.all)[2]
    # First step
    foright <- paste(names(nonEventDF1)[-c(1:5)],collapse = "+")
    fo <- paste("survival::Surv(Start, Stop, Status) ~ ",foright, sep = "")

    fit1 <- survival::coxph(as.formula(fo), data = nonEventDF1)
    alpha.est <- as.vector(coef(fit1))

    # Second step

    # With 2, there is no visit at time 0
    u2.all <- u.all[u.all>0]
    t2.all <- nonEventDF1$Stop[nonEventDF1$Status==1] # same as t.all[t.all>0]

    expXa <- exp(-xt2.all%*%(alpha.est))

    s1arr <- array(dim = c(2, length(u2.all), length(t2.all)))#matrix(nrow = length(u2.all), ncol = length(t2.all))
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
    list(beta = beta.est, alpha = alpha.est)
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

  res <- matrix(ncol = P, nrow = B)
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
    ut.all2b <- nonEventDF1b$Stop
    ind.allb <- nonEventDF1b$Status
    ind2.allb <- nonEventDF1b$Status2
    xw.allb <- nonEventDF1b[, -c(1:5)]


    id.all1b <- rep(1:N, mu[bind])
    id.all2b <- rep(1:N, mt[bind])
    id.all3b <- rep(1:N, mdat[bind])


    u.allb <- eventDFb$Time
    t.allb <- nonEventDF2b$Time
    zwu.allb <- eventDFb[,-c(1:2)]
    zwt.allb <- nonEventDF2b[,-c(1:2)]


    datb <- data.frame(ID = id.all3b, Start = ut.allb,
                       Stop = ut.all2b, Status = ind.allb,
                       Status2 = ind2.allb)
    datb <- cbind(datb, xw.allb)


    fitb <- est(datb, u.allb, t.allb, zwu.allb, zwt.allb, tau)

    res[b,] <- c(fitb$beta)
  }

  res
}

#set.seed(1)
#se2 <- reVARBoot(dat$nonEventDF1,dat$nonEventDF2,dat$eventDF,4.5, h = 2*N^(-1/3), B = 100)

