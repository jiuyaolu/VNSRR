setwd("/Users/jiuyaolu/Desktop/rerandomization/VNSRR/")

library(MASS)
library(tidyverse)
library(parallel)


set.seed(2023)

ngroup_vec = c(2)
R2_vec = c(0.5)
p_vec = c(50, 250)

n_max = max(ngroup_vec) * max(p_vec) * 2
p_max = max(p_vec)
X_full = mclapply(1:1000, function(i) {
  mvrnorm(n_max, 
          mu = rep(0, p_max), 
          Sigma = diag(p_max),
          empirical = F)
}, mc.cores = 10)

for (R2 in R2_vec) {
  for (ngroup in ngroup_vec) {
    for (p in p_vec) {
      n = ngroup * p * 2
      
      X = mclapply(1:1000, function(i) {
        X_full[[i]][1:n, 1:p]
      }, mc.cores = 10)
      
      Y0 = mclapply(1:1000, function(i) {
        fx = rowSums(X[[i]])
        var_e = (1 - R2) / R2 * p
        as.vector(fx + mvrnorm(n, mu = 0, Sigma = var_e, empirical = F))
      }, mc.cores = 10)
      Y1 = mclapply(1:1000, function(i) {
        Y0[[i]] + 0.3 * sqrt(p / R2)
      }, mc.cores = 10)
      save(X, Y0, Y1, file = sprintf("data/sim_sequential/sim_sequential_n%s_p%s.RData",n,p))
    }
  }
}
