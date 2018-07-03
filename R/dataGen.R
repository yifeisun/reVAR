#' Generate a simulated dataset
#'
#' This function generate a simulated dataset under the first scenario in the
#' simulation study of the paper.
#'
#' @param N Number of observations that will be generated
#' @return A list of simulated data, contains the data frames that will be used in the function reVAR
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)

genData <- function(N)
{
  genOneObs <- function(cen, xi, z0, w, gamma)
  {

    # For event visit
    Lam <- function(t, t0, gamma, si, w, z0, lv)
    {
      sit <- c(si[si<t], t)
      sit0 <- c(si[si<t0], t0)
      if(z0 == 1)
      {
        exp1 <- exp(beta1+beta2*w+gamma)
        exp2 <- exp(gamma+beta2*w)
        Lamt <- sum(diff(sit^2)*rep(c(exp1,exp2), length = length(sit)-1))
        Lamt0 <- sum(diff(sit0^2)*rep(c(exp1,exp2), length = length(sit0)-1))
      }else{
        exp2 <- exp(beta1+gamma+beta2*w)
        exp1 <- exp(gamma+beta2*w)
        Lamt <- sum(diff(sit^2)*rep(c(exp1,exp2), length = length(sit)-1))
        Lamt0 <- sum(diff(sit0^2)*rep(c(exp1,exp2), length = length(sit0)-1))
      }
      Lamt-Lamt0 + lv
    }
    # Generate recurrent event data
    # z0 is Z(0)
    # gamma is the frailty
    # si are the time points that Z(t) changes value
    # cen is censoring time
    # lv is log(V)
    # w is baseline covariate

    # For regular visit
    Lam2 <- function(t, t0, si, w, ui, zui, z0, lv)
    {
      zt0 <- ifelse(sum(si<=t0)%%2 == 0, 1-z0, z0)
      t0ui <- c(t0, ui[ui>t0 & ui<t], t)
      zt0ui <- c(zt0, zui[ui>t0 & ui<t])
      Lamt_t0 <- sum(l0*exp(alpha1*w+alpha2*zt0ui)*diff(t0ui))
      Lamt_t0 + lv
    }
    # ui are event times

    # Generate Z(t)
    si <- 0
    while(tail(si,1) < cen)
    {
      snew <- tail(si,1) + rexp(1, xi)
      if(snew < cen)
      {
        si <- c(si, snew)
      }else{ break }
    }

    # Generate event visits
    ui <- 0
    zui <- z0
    while(tail(ui,1) < cen)
    {
      v <- runif(1)
      uold <- tail(ui,1)
      unew <- uniroot(Lam, c(10^{-10}, 10^10),
                      lv = log(v), t0 = uold,
                      gamma = gamma, si = si,
                      z0 = z0, w = w)$root
      if(unew < cen)
      {
        ui <- c(ui, unew)
        znew <- ifelse(sum(si<=unew)%%2 == 0, 1-z0, z0)
        zui <- c(zui, znew)
      }else{ break }
    }


    # Generate regular visit (also using covariates observed at event visit)
    ti <- 0
    zti <- z0
    while(tail(ti,1) < cen)
    {
      v <- runif(1)
      told <- tail(ti,1)
      tnew  <- uniroot(Lam2, c(10^{-10}, 10^10),
                       lv = log(v), t0 = told, si = si,
                       ui = ui, zui = zui,
                       z0 = z0, w = w)$root

      if(tnew < cen)
      {
        ti <- c(ti, tnew)
        znew <- ifelse(sum(si<=tnew)%%2 == 0, 1-z0, z0)
        zti <- c(zti, znew)
      }else{ break }
    }
    list(zui = zui, ui = ui,
         zti = zti, ti = ti,
         z0 = z0, si = si)
  }

  # true parameters
  beta1 <- 0.5
  beta2 <- 0.5
  alpha1 <- .5
  alpha2 <- .5
  mugamma <- -1
  l0 <- exp(0)
  num <- 1


  C <- runif(N, 0, 5)
  # For generating Z(t)
  Xi <- rgamma(N, 5, 5)
  Z0 <- rbinom(N, 1, 0.2)
  Gamma <- rnorm(N, mean = mugamma, sd = 0.5)
  W <- runif(N, -.5, .5)

  u.all <- list()
  # time of all the event visits, including 0
  t.all <- list()
  # times of all the regular visits, including 0
  zu.all <- list()
  # covariates at event visits, including z0
  zt.all <- list()
  # covariates at regular visits, including z0
  id.all1 <- list()
  # id of event visits
  id.all2 <- list()
  # id of regular visits
  id.all3 <- list()
  # id of all the visits mixed together
  ut.all <- list()
  # times of all the visits mixed together, including 0, for Start in Surv()
  ut.all2 <- list()
  # time of all the visits mixed together, including C, for Stop in Surv()
  x.all <- list()
  # X at all the Start time points (ut.all)
  ind.all <- list()
  ind2.all <- list()
  # Whether the mixed event (including C) is REGULAR or not, for Status in Surv()

  wu.all <- list()
  wt.all <- list()
  # Baseline covariates
  x2.all <- list()
  # X2
  w.all3 <- list()
  # for dat.AG

  for(i in 1:N)
  {
    dat <- genOneObs(C[i],Xi[i],Z0[i],W[i],Gamma[i])
    z0 <- Z0[i]
    if(sum(dat$ui[-1]<=1e-6)>0)
    {
      dat$ui[dat$ui>0 & dat$ui<=1e-6] <- 1e-5
    }
    if(sum(dat$ti[-1]<=1e-6)>0)
    {
      dat$ti[dat$ti>0 & dat$ti<=1e-6] <- 1e-5
    }
    uti <- c(0, dat$ui[-1], dat$ti[-1])
    # event or non-event
    indi <- c(0, rep(0, length(dat$ui)-1), rep(1, length(dat$ti)-1))
    zuti <- c(z0, dat$zui[-1], dat$zti[-1])

    ind.all[[i]] <- c((indi[order(uti)])[-1], 0)
    ind2i <- c(0, rep(1, length(dat$ui)-1), rep(0, length(dat$ti)-1))
    ind2.all[[i]] <- c((ind2i[order(uti)])[-1], 0)
    uti2 <- uti[order(uti)]

    x.all[[i]] <- zuti[order(uti)]
    ut.all[[i]] <- uti2
    ut.all2[[i]] <- c(uti2[-1],C[i])
    id.all3[[i]] <- rep(i, length(uti2))
    w.all3[[i]] <- rep(W[i],length(uti2))

    u.all[[i]] <- dat$ui
    t.all[[i]] <- dat$ti
    zu.all[[i]] <- dat$zui
    zt.all[[i]] <- dat$zti
    id.all1[[i]] <- rep(i, length(dat$ui))
    id.all2[[i]] <- rep(i, length(dat$ti))
    wu.all[[i]] <- rep(W[i], length(dat$ui))
    wt.all[[i]] <- rep(W[i], length(dat$ti))
  }

  ut.allv <- unlist(ut.all)
  ut.all2v <- unlist(ut.all2)
  ind.allv <- unlist(ind.all)
  ind2.allv <- unlist(ind2.all)

  x.allv <- unlist(x.all)

  id.all1v <- unlist(id.all1)
  id.all2v <- unlist(id.all2)
  id.all3v <- unlist(id.all3)
  w.all3v <- unlist(w.all3)

  u.allv <- unlist(u.all)
  t.allv <- unlist(t.all)
  zu.allv <- unlist(zu.all)
  zt.allv <- unlist(zt.all)
  wu.allv <- unlist(wu.all)
  wt.allv <- unlist(wt.all)

  nonEventDF1 <- data.frame(ID = id.all3v, Start = ut.allv,
                            Stop = ut.all2v, Status = ind.allv,
                            Status2 = ind2.allv,
                            X = x.allv, W = w.all3v)
  # nonEventDF1 is the data frame to be passed into coxph for the estimation of the non-event visit model

  mu <- table(id.all1v)
  # mu-1 is the number of event visits for N subjects
  mt <- table(id.all2v)
  # mt-1 is the number of non-event visits for N subjects
  mdat <- table(nonEventDF1$ID)
  # mdat-1 is the number of visits for N subjects

  nonEventDF2 <- data.frame(ID = id.all2v,
                            Time = t.allv,
                            Z = zt.allv,
                            W = wt.allv)

  eventDF <- data.frame(ID = id.all1v,
                        Time = u.allv,
                        Z = zu.allv,
                        W = wu.allv)

  list(nonEventDF1 = nonEventDF1, nonEventDF2 = nonEventDF2, eventDF = eventDF)
}




