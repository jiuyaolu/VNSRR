setwd("/Users/jiuyaolu/Desktop/rerandomization/VNSRR/")

library(MASS)
library(tidyverse)
library(parallel)


set.seed(2023)

n_vec = c(30, 100, 500)
R2_vec = c(0.5)
rp_vec = c(0.5)

n_max = max(n_vec)
p_max = round(n_max * max(rp_vec))
X_full = mclapply(1:1000, function(i) {
  mvrnorm(n_max, 
          mu = rep(0, p_max), 
          Sigma = diag(p_max),
          empirical = T)
}, mc.cores = 10)

for (R2 in R2_vec) {
  for (n in n_vec) {
    for (rp in rp_vec) {
      p = round(n * rp)
      X = mclapply(1:1000, function(i) {
        tmp = X_full[[i]][1:n, 1:p]
        tmp <- scale(tmp, TRUE, FALSE)
        tmp <- tmp %*% svd(tmp, nu = 0)$v
        tmp <- scale(tmp, FALSE, TRUE)
        tmp
      }, mc.cores = 10)
      
      Y0 = mclapply(1:1000, function(i) {
        fx = rowSums(X[[i]])
        var_e = (1 - R2) / R2 * var(fx)
        as.vector(fx + mvrnorm(n, mu = 0, Sigma = var_e, empirical = T))
      }, mc.cores = 10)
      Y1 = mclapply(1:1000, function(i) {
        Y0[[i]] + 0.3 * sd(Y0[[i]])
      }, mc.cores = 10)
      save(X, Y0, Y1, file = sprintf("data/sim_simple/sim_simple_n%s_p%s.RData",n,p))
    }
  }
}
