#' Fit the proportional rate model with intermittently measured time-dependent covariates under VAR
#'
#' This function fit the proportional rate model with intermittently
#' measured time-dependent covariates under VAR.
#' @param nonEventDF1 A data frame that is used for fitting the non-event time model, contains ID, event time information and covariates
#' @param nonEventDF2 A data frame that contains the ID, non-event visit times and covariates measured at each non-event visit
#' @param eventDF A data frame that contains the ID, event visit times and covariates measured at each event visit
#' @param tau Data on [0,tau] will be used in the estimating equation in the second step
#' @param h Bandwidth for kernel estimation
#' @return A list of parameter estimates in the first and second step
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)
#' fit <- reVAR(dat$nonEventDF1,dat$nonEventDF2,dat$eventDF, tau = 4.5, h = 2*100^(-1/3))



reVAR <- function(nonEventDF1, nonEventDF2, eventDF, tau, h)
{
  u.all <- eventDF$Time
  t.all <- nonEventDF2$Time

  zwu.all <- as.matrix(eventDF[-c(1:2)])
  xt2.all <- as.matrix(nonEventDF1[nonEventDF1$Status==1,-c(1:5)])
  zwt2.all <- as.matrix(nonEventDF2[t.all>0,-c(1:2)])

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

