source("./simulations/functions_simulation.R")

library(doParallel)
library(doRNG)

cl <- makeCluster(12, outfile = "") # for use on a cluster with 12 nodes
registerDoParallel(cl)
registerDoRNG(1)


# The following code generate the results under Scenario I
# For other scenarios, change "VCAR1" to other values

err <- foreach(j = 1:1000,
               .combine = rbind,
               .packages = c("survival", "NHPoisson", "nleqslv")) %dorng%
{
  tryCatch(oneiter("VCAR1", B = 500, N = 200, tau = 4.5, h = 0.43),
           error = function(e) rep(NA, 27))
}


# bias
# the proposed method
beta <- c(-1,1,-1)
colMeans(err[,1:3]) - beta
# Li et al. (2016)
colMeans(err[,4:6]) - beta
# LCCF
colMeans(err[,7:9]) - beta

# SE
# the proposed method
apply(err[,1:3], 2, sd)
# Li et al. (2016)
apply(err[,4:6], 2, sd)
# LCCF
apply(err[,7:9], 2, sd)

# average estimated standard errors
# the proposed methods
colMeans(err[,10:12])
# Li et al. (2016)
colMeans(err[,13:15])
# LCCF
colMeans(err[,16:18])

# coverage probability
# the proposed method
colMeans(err[,19:21])
# Li et al. (2016)
colMeans(err[,22:24])
# LCCF
colMeans(err[,25:27])



