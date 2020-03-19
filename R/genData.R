#' Generate a simulated dataset
#'
#' This function generates a simulated dataset under scenarios of the
#' simulation studies in Sun et al. (2020+).
#'
#' @param N Integer value: number of observations that will be generated
#' @param sce Character value: The scenario number of simulation studies in Sun et al. (2020+).
#' Possible values are "VCAR1", "VAR2", "VAR3", "VCAR4", "VAR-M5", "VAR-M6", "VNAR7", "VNAR8", "VNAR9", "VNAR10".
#'
#' @return A list of simulated data, contains the three data frames that will be used in the function reVAR
#' \item{nonEventDF1_s2}{A data frame that is used for fitting the non-event time model, contains ID, event time information and covariates.
#' Each row records the covariates on a time interval.
#' \code{ID} is the subject ID for each observation;
#' \code{Start} is the starting time of the interval;
#' \code{End} is the ending time of the interval;
#' \code{Status} is the status indicator, 1 = non-event visit, 0 = no;
#' \code{Status2} is the status indicator, 1 = event visit, 0 = no.
#' The other columns are covariate values for the non-event model on the time interval (\code{Start}, \code{End}].}
#'
#'
#' \item{nonEventDF2}{A data frame that contains the ID, non-event visit times and covariates (for the event model) measured at each non-event visit.
#' \code{ID} is the subject ID for each observation;
#' \code{Time} are time zero and the non-event visit times.
#' The other columns are covariate values measured at \code{Time}.}
#'
#'
#' \item{eventDF}{A data frame that contains the ID, event visit times and covariates (for the event model) measured at each event visit.
#' \code{ID} is the subject ID for each observation;
#' \code{Time} are time zero and the event times.
#' The other columns are covariate values measured at \code{Time}.}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- genData(100)
#' str(dat)


genData <- function(N, sce = "VCAR1")
{
  getOneVNAR <- function(cen, xi, beta, alpha, gamma)
  {
    # Generate Z_1(t)
    si <- 0
    while(tail(si,1) < cen)
    {
      snew <- tail(si,1) + rexp(1, xi)
      if(snew < cen)
      {
        si <- c(si, snew)
      }else{ break }
    }

    z0 <- rbinom(1,1,0.5)

    tg <- (1:10000)/10000*5
    tg <- round(tg[tg<=cen],4)
    #as.numeric(cut(tg, c(si,Inf), right = FALSE))%%2

    # Z on time grid
    Z1tg <- (z0==1)*as.numeric(cut(tg, c(si,Inf), right = FALSE))%%2 +
      (z0==0)*(1-as.numeric(cut(tg, c(si,Inf), right = FALSE))%%2)

    pic <- pi
    cz2 <- runif(1,0,2*pi)
    cw <- runif(1,0,2*pi)

    Z2tg <- sin(tg*pic+cz2)
    Wtg <- sin(tg*pic+cw)

    Z3 <- runif(1,-.5,.5)

    # event process
    mutg <- tg*exp(beta[1]*Z1tg+beta[2]*Z2tg+beta[3]*Z3+gamma[1]*Wtg - 1)/10000*5

    ui <- NHPoisson::simNHP.fun(mutg)$posNH/10000*5
    ui <- c(0,ui)
    ui <- unique(ui)
    Z1ui <- (z0==1)*as.numeric(cut(ui, c(si,Inf), right = FALSE))%%2 +
      (z0==0)*(1-as.numeric(cut(ui, c(si,Inf), right = FALSE))%%2)

    Z2ui <- sin(ui*pic+cz2)
    Z3ui <- rep(Z3, length(ui))
    #Z1ui[1] <- z0

    ti <- 0
    # the time of last visit
    utlast <- 0
    while(utlast < cen)
    {
      uilast <- tail(ui[ui<=utlast],1)
      uinext <- ui[ui>utlast][1]
      if(is.na(uinext))
      {
        uinext <- cen
      }
      Z1last <- (z0==1)*as.numeric(cut(utlast, c(si,Inf), right = FALSE))%%2 +
        (z0==0)*(1-as.numeric(cut(utlast, c(si,Inf), right = FALSE))%%2)

      Z2last <- sin(utlast*pic+cz2)
      tgsub <- tg[tg<uinext-1e-6 & tg>utlast+1e-6]

      Z1tgsub <- (z0==1)*as.numeric(cut(tgsub, c(si,Inf), right = FALSE))%%2 +
        (z0==0)*(1-as.numeric(cut(tgsub, c(si,Inf), right = FALSE))%%2)
      Z2tgsub <- sin(tgsub*pic+cz2)
      Wtgsub <- sin(tgsub*pic+cw)

      lambdatgsub <- 1*exp(alpha[1]*Z1last+alpha[2]*Z2last+alpha[3]*Z3+alpha[4]*Z1tgsub+
                             alpha[5]*Z2tgsub+gamma[2]*Wtgsub)/10000*5
      nhp <- NHPoisson::simNHP.fun(lambdatgsub)$posNH
      if(length(nhp)>0)
      {# tinew < next ui
        tinew <- utlast + nhp[1]/10000*5
        #if(tinew > cen-1e-6 | tinew > uinext-1e-6) browser()
        ti <- c(ti, tinew)
        utlast <- tinew
      }else{
        # no non-event visit
        utlast <- uinext
      }
    }

    Z3ti <- rep(Z3, length(ti))
    Z1ti <- (z0==1)*as.numeric(cut(ti, c(si,Inf), right = FALSE))%%2 +
      (z0==0)*(1-as.numeric(cut(ti, c(si,Inf), right = FALSE))%%2)

    Z2ti <- sin(ti*pic+cz2)
    #Z1ti[1] <- z0

    list(ui = ui, zui = cbind(Z1ui,Z2ui,Z3ui),
         ti = ti, zti = cbind(Z1ti,Z2ti,Z3ti))
  }

  beta <- c(-1,1,-1)
  if(sce == "VCAR1")
  {
    alpha <- c(0,0,0,0,0)
    gamma <- c(0,0)
  }else if(sce == "VAR2")
  {
    alpha <- c(-1,1,-1,0,0)/2
    gamma <- c(0,0)
  }else if(sce == "VAR3")
  {
    alpha <- c(-1,1,-1,0,0)
    gamma <- c(0,0)
  }else if(sce == "VCAR4")
  {
    alpha <- c(0,0,0,0,0)
    gamma <- c(1,1)
  }else if(sce == "VAR-M5")
  {
    alpha <- c(-1,1,-1,0,0)/2
    gamma <- c(1,1)
  }else if(sce == "VAR-M6")
  {
    alpha <- c(-1,1,-1,0,0)
    gamma <- c(1,1)
  }else if(sce == "VNAR7")
  {
    alpha <- c(-1,1/2,-1,0,1/2)
    # X1,X2,Z3,Z1,Z2
    gamma <- c(1,1)
  }else if(sce == "VNAR8")
  {
    alpha <- c(-1/2,1,-1,-1/2,0)
    gamma <- c(1,1)
  }else if(sce == "VNAR9")
  {
    alpha <- c(-1/2,1/2,-1,-1/2,1/2)
    gamma <- c(1,1)
  }else{
    alpha <- c(0,0,-1,-1/2,1/2)
    gamma <- c(1,1)
  }

  #############################################
  # Generate data
  #############################################
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

  C <- runif(N,0,5)
  xi <- rgamma(N,5,5)

  for(i in 1:N)
  {
    dat <- getOneVNAR(C[i],xi[i], beta, alpha, gamma)
    uti <- c(0, dat$ui[-1], dat$ti[-1])
    # event or non-event
    indi <- c(0, rep(0, length(dat$ui)-1), rep(1, length(dat$ti)-1))
    zuti <- rbind(dat$zui[1,], dat$zui[-1,], dat$zti[-1,])

    ind.all[[i]] <- c((indi[order(uti)])[-1], 0)
    ind2i <- c(0, rep(1, length(dat$ui)-1), rep(0, length(dat$ti)-1))
    ind2.all[[i]] <- c((ind2i[order(uti)])[-1], 0)
    uti2 <- uti[order(uti)]

    x.all[[i]] <- zuti[order(uti),]
    ut.all[[i]] <- uti2
    ut.all2[[i]] <- c(uti2[-1],C[i])
    id.all3[[i]] <- rep(i, length(uti2))

    u.all[[i]] <- dat$ui
    t.all[[i]] <- dat$ti
    zu.all[[i]] <- dat$zui
    zt.all[[i]] <- dat$zti
    id.all1[[i]] <- rep(i, length(dat$ui))
    id.all2[[i]] <- rep(i, length(dat$ti))
    #print(i)
  }

  #############################################
  # Prepare data frames for reVAR
  #############################################

  ut.allv <- unlist(ut.all)
  ut.all2v <- unlist(ut.all2)
  ind.allv <- unlist(ind.all)
  ind2.allv <- unlist(ind2.all)

  x.allm <- Reduce(rbind, x.all)
  colnames(x.allm) <- c("X1","X2","X3")
  id.all1v <- unlist(id.all1)
  id.all2v <- unlist(id.all2)
  id.all3v <- unlist(id.all3)


  u.allv <- unlist(u.all)
  t.allv <- unlist(t.all)
  zu.allm <- Reduce(rbind, zu.all)
  zt.allm <- Reduce(rbind, zt.all)


  nonEventDF1_s2 <- data.frame(ID = id.all3v,
                               Start = ut.allv,
                               End = ut.all2v,
                               Status = ind.allv,
                               Status2 = ind2.allv,
                               x.allm)
  # nonEventDF1 is the data frame to be passed into coxph for the estimation of the non-event visit model
  # Status2 is not used in the method proposed in Sun et al.

  nonEventDF2 <- data.frame(ID = id.all2v,
                            Time = t.allv,
                            zt.allm)

  eventDF <- data.frame(ID = id.all1v,
                        Time = u.allv,
                        zu.allm)

  rownames(nonEventDF1_s2) <- NULL
  rownames(nonEventDF2) <- NULL
  rownames(eventDF) <- NULL

  list(nonEventDF1_s2 = nonEventDF1_s2,
       nonEventDF2 = nonEventDF2,
       eventDF = eventDF)
}


