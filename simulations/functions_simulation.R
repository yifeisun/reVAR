###############################################################
#         Functions used in the simulation studies
###############################################################

# genData2() is slightly different from genData(). The difference is that
# nonEventDF1 returned by genData2() contains an additional column of event
# indicator, which is used in the LCCF method

genData2 <- function(N, sce = "VCAR1")
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

    ui <- simNHP.fun(mutg)$posNH/10000*5
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
  # time of all the visits mixed together, including C, for End in Surv()
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


  nonEventDF1 <- data.frame(ID = id.all3v,
                            Start = ut.allv,
                            End = ut.all2v,
                            Status = ind.allv,
                            Status2 = ind2.allv,
                            x.allm)
  # nonEventDF1 is the data frame to be passed into coxph for the estimation of the non-event visit model


  nonEventDF2 <- data.frame(ID = id.all2v,
                            Time = t.allv,
                            zt.allm)

  eventDF <- data.frame(ID = id.all1v,
                        Time = u.allv,
                        zu.allm)

  rownames(nonEventDF1) <- NULL
  rownames(nonEventDF2) <- NULL
  rownames(eventDF) <- NULL

  list(nonEventDF1 = nonEventDF1,
       nonEventDF2 = nonEventDF2,
       eventDF = eventDF)
}

est <- function(nonEventDF1, u.all, t.all, zu.all, zt.all, tau, h)
{
  K <- function(x)
  {
    (abs(x)<1)*0.75*(1-(x)^2)
  }

  xt2.all <- as.matrix(nonEventDF1[nonEventDF1$Status==1,-c(1:5)])
  zt2.all <- as.matrix(zt.all[t.all>0,])

  P <- dim(zu.all)[2]
  # First step
  foright <- paste(names(nonEventDF1)[-c(1:5)],collapse = "+")
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
    s1arr[p,,] <- t(t(s0mat)*zt2.all[,p])
  }

  U <- function(b)
  {
    expZb <- as.vector(exp(zt2.all%*%b))
    fun1 <- function(matp)
    {
      colSums(t(matp)*expZb)/colSums(t(s0mat)*expZb)
    }
    ee <- apply(s1arr,1,fun1)
    ee[is.na(ee)] <- 0
    colSums(zu.all[u.all<=tau & u.all>0, ])- colSums(ee[u2.all<=tau,])
  }

  beta.est <- nleqslv::nleqslv(rep(0,P),U)$x
  list(beta = beta.est, alpha = alpha.est)
}

est0 <- function(nonEventDF1, u.all, t.all, zu.all, zt.all, tau, h)
{
  K <- function(x)
  {
    (abs(x)<1)*0.75*(1-(x)^2)
  }

  xt2.all <- as.matrix(nonEventDF1[nonEventDF1$Status==1,-c(1:5)])
  zt2.all <- as.matrix(zt.all[t.all>0,])

  P <- dim(zu.all)[2]
  # Second step

  # With 2, there is no visit at time 0
  u2.all <- u.all[u.all>0]
  t2.all <- nonEventDF1$End[nonEventDF1$Status==1] # same as t.all[t.all>0]



  s1arr <- array(dim = c(P, length(u2.all), length(t2.all)))
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
    s1arr[p,,] <- t(t(s0mat)*zt2.all[,p])
  }

  U <- function(b)
  {
    expZb <- as.vector(exp(zt2.all%*%b))
    fun1 <- function(matp)
    {
      colSums(t(matp)*expZb)/colSums(t(s0mat)*expZb)
    }
    ee <- apply(s1arr,1,fun1)
    ee[is.na(ee)] <- 0
    colSums(zu.all[u.all<=tau & u.all>0, ])- colSums(ee[u2.all<=tau,])
  }

  beta.est <- nleqslv::nleqslv(rep(0,P),U)$x
  list(beta = beta.est)
}

reVARBoot2 <- function(nonEventDF1, nonEventDF2, eventDF, tau, h, B = 500)
{
  mu <- table(eventDF$ID)
  # mu-1 is the number of event visits for N subjects
  mt <- table(nonEventDF2$ID)
  # mt-1 is the number of non-event visits for N subjects
  mdat <- table(nonEventDF1$ID)
  # mdat-1 is the number of visits for N subjects
  N <- length(unique(eventDF$ID))
  P <- dim(eventDF)[2]-2

  res <- matrix(ncol = P*2, nrow = B)
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
    x.allb <- nonEventDF1b[, -c(1:5)]


    id.all1b <- rep(1:N, mu[bind])
    id.all2b <- rep(1:N, mt[bind])
    id.all3b <- rep(1:N, mdat[bind])


    u.allb <- eventDFb$Time
    t.allb <- nonEventDF2b$Time
    zu.allb <- eventDFb[,-c(1:2)]
    zt.allb <- nonEventDF2b[,-c(1:2)]


    datb <- data.frame(ID = id.all3b, Start = ut.allb,
                       End = ut.all2b, Status = ind.allb,
                       Status2 = ind2.allb)
    datb <- cbind(datb, x.allb)


    fitb <- est(datb, u.allb, t.allb, zu.allb, zt.allb, tau, h)
    fitb0 <- est0(datb, u.allb, t.allb, zu.allb, zt.allb, tau, h)

    res[b,] <- c(fitb$beta,fitb0$beta)
  }

  res
}

oneiter <- function(sce = "VCAR1", B = 500, N, tau, h)
{
  dat <- genData2(N, sce)
  nonEventDF1 <- dat$nonEventDF1
  nonEventDF2 <- dat$nonEventDF2
  eventDF <- dat$eventDF

  u.all <- eventDF$Time
  t.all <- nonEventDF2$Time
  zu.all <- eventDF[,-c(1:2)]
  zt.all <- nonEventDF2[,-c(1:2)]


  fit <- est(nonEventDF1, u.all, t.all, zu.all, zt.all, tau, h)
  fit0 <- est0(nonEventDF1, u.all, t.all, zu.all, zt.all, tau, h)

  fit.LCCF <- coxph(Surv(Start, End, Status2) ~ X1 + X2 + X3 + cluster(ID),
                    data = dat$nonEventDF1)


  bt <- reVARBoot2(nonEventDF1, nonEventDF2, eventDF, tau, h, B)
  se <- apply(bt, 2, sd)

  beta <- c(-1,1,-1)

  est1 <- fit$beta
  L1 <- fit$beta[1] - 1.96*se[1]
  U1 <- fit$beta[1] + 1.96*se[1]
  L2 <- fit$beta[2] - 1.96*se[2]
  U2 <- fit$beta[2] + 1.96*se[2]
  L3 <- fit$beta[3] - 1.96*se[3]
  U3 <- fit$beta[3] + 1.96*se[3]
  cr <- c(L1<= beta[1] & U1>=beta[1],
          L2<= beta[2] & U2>=beta[2],
          L3<= beta[3] & U3>=beta[3])

  est2 <- fit0$beta
  L10 <- fit0$beta[1] - 1.96*se[4]
  U10 <- fit0$beta[1] + 1.96*se[4]
  L20 <- fit0$beta[2] - 1.96*se[5]
  U20 <- fit0$beta[2] + 1.96*se[5]
  L30 <- fit0$beta[3] - 1.96*se[6]
  U30 <- fit0$beta[3] + 1.96*se[6]
  cr0 <- c(L10<= beta[1] & U10>=beta[2],
           L20<= beta[2] & U20>=beta[3],
           L30<= beta[3] & U30>=beta[3])

  s.LCCF <- summary(fit.LCCF)$coef

  L1.LCCF <- s.LCCF[1,1] - 1.96*s.LCCF[1,4]
  U1.LCCF <- s.LCCF[1,1] + 1.96*s.LCCF[1,4]
  L2.LCCF <- s.LCCF[2,1] - 1.96*s.LCCF[2,4]
  U2.LCCF <- s.LCCF[2,1] + 1.96*s.LCCF[2,4]
  L3.LCCF <- s.LCCF[3,1] - 1.96*s.LCCF[3,4]
  U3.LCCF <- s.LCCF[3,1] + 1.96*s.LCCF[3,4]
  seLCCF <- c(s.LCCF[1,4], s.LCCF[2,4], s.LCCF[3,4])
  crLCCF <- c(L1.LCCF<= beta[1] & U1.LCCF>=beta[1],
              L2.LCCF<= beta[2] & U2.LCCF>=beta[2],
              L3.LCCF<= beta[3] & U3.LCCF>=beta[3])

  estLCCF <- coef(fit.LCCF)

  c(est1,est2,estLCCF,se,seLCCF,cr,cr0,crLCCF)
  # returns the point estimate/SE/coverage(yes or no) of 95% CI of the proposed method/Li et al. (2016)/LCCF
}

K <- function(x)
{
  (abs(x)<1)*0.75*(1-(x)^2)
}

