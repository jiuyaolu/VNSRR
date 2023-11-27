# compute Mahalanobis distance --------------------------------------------

compute_M_set <- function(Wt, X, prop_cov_inv = NULL) {
  if (is.null(prop_cov_inv)) {
    n <- nrow(X)
    n1 <- length(Wt)
    prop_cov_inv <- n1 * (1 - n1 / n) * solve(cov(X))
  }
  dif <- colMeans(X[Wt, ]) - colMeans(X[-Wt, ])
  as.vector(t(dif) %*% prop_cov_inv %*% dif)
}

compute_Mstar_set <- function(Wt, Wstar, i, j, Mt, H, h) {
  Mt + (2 * sum(H[j, Wstar]) - H[j, j]) - 
    (2 * sum(H[i, Wt]) - H[i, i]) + h[i] - h[j]
}

compute_M <- function(W, X, prop_cov_inv = NULL) {
  Wt <- which(W == 1)
  compute_M_set(Wt, X, prop_cov_inv)
}


# # test
# # three functions output same results
# # speed: compute_Mstar_set > compute_M_set > compute_M
# n <- 100
# p <- 10
# X <- matrix(rnorm(n * p), n, p)
# n1 <- round(n / 2)
# cov_inv <- solve(cov(X))
# prop_cov_inv <- n1 * (1 - n1 / n) * cov_inv
# H <- X %*% cov_inv %*% t(X) / (n1 * (1 - n1 / n))
# h <- 2 * (n1 / n) * rowSums(H)
# Wt <- sample(1:n, n1) # original assignment
# Mt <- compute_M_set(Wt, X, prop_cov_inv)
# Wstar <- Wt
# i <- sample(Wt, 1)
# j <- sample(setdiff(1:n, Wt), 1)
# Wstar[which(Wstar == i)] <- j # pair-switched assignment
# W <- rep(0, n)
# W[Wstar] <- 1 # pair-switched assignment with binary representation
# bench::mark(
#   compute_M_set(Wstar, X, prop_cov_inv),
#   compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h),
#   compute_M(W, X, prop_cov_inv),
#   compute_M(W, X)
# )


# rerandomization ---------------------------------------------------------

rerand_one <- function(X, n, n1, prop_cov_inv, H, h,
                       a, pair_switch, temper,
                       max_iter, report_iter) {
  # start point
  Wt <- sample(1:n, n1)
  if (a == Inf) {
    cont <- FALSE
    Mt <- NA
    Wbest <- Wt
    Mbest <- Mt
  } else {
    Mt <- compute_M_set(Wt, X, prop_cov_inv)
    n_iter <- 1
    n_jump <- 1
    Wbest <- Wt
    Mbest <- Mt
    if (Mt <= a) {
      cont <- FALSE
    } else {
      cont <- TRUE
    }
  }
  
  # main ----
  if (!report_iter) {
    if (max_iter == Inf) {
      # 1 rerand, not limit number of iterations
      if (!pair_switch) {
        # 1.1 independent rerand
        while (cont) {
          Wt <- sample(1:n, n1)
          Mt <- compute_M_set(Wt, X, prop_cov_inv)
          if (Mt <= a) cont <- FALSE
        }
      } else {
        # 1.2 pair-switching rerand
        while (cont) {
          Wstar <- Wt
          i <- sample(Wt, 1)
          j <- sample(setdiff(1:n, Wt), 1)
          Wstar[which(Wstar == i)] <- j
          Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
          if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
            Wt <- Wstar
            Mt <- Mstar
          }
          if (Mt <= a) cont <- FALSE
        }
      }
    } else {
      # 2 rerand, limit number of iterations 
      if (!pair_switch) {
        # 2.1 independent rerand
        if (cont) {
          for (iter in 2:max_iter) {
            Wt <- sample(1:n, n1)
            Mt <- compute_M_set(Wt, X, prop_cov_inv)
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) break
          }
        }
      } else {
        # 2.2 pair-switching rerand
        if (cont) {
          for (iter in 2:max_iter) {
            Wstar <- Wt
            i <- sample(Wt, 1)
            j <- sample(setdiff(1:n, Wt), 1)
            Wstar[which(Wstar == i)] <- j
            Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
            if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
              Wt <- Wstar
              Mt <- Mstar
            }
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) break
          }
        }
      }
      # output the best rather than the last
      Wt <- Wbest
      Mt <- Mbest
    }
    # output binary representation of Wt
    Wout <- rep(0, n)
    Wout[Wt] <- 1
    return(Wout)
  } else {
    # A record iterations (for developer) ----
    if (max_iter == Inf) {
      # A1 rerand, not limit number of iterations
      if (!pair_switch) {
        # A1.1 independent rerand
        while (cont) {
          Wt <- sample(1:n, n1)
          Mt <- compute_M_set(Wt, X, prop_cov_inv)
          n_iter <- n_iter + 1
          if (Mt <= a) cont <- FALSE
        }
        n_jump <- n_iter
      } else {
        # A1.2 pair-switching rerand
        while (cont) {
          Wstar <- Wt
          i <- sample(Wt, 1)
          j <- sample(setdiff(1:n, Wt), 1)
          Wstar[which(Wstar == i)] <- j
          Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
          n_iter <- n_iter + 1
          if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
            Wt <- Wstar
            Mt <- Mstar
            n_jump <- n_jump + 1
          }
          if (Mt <= a) cont <- FALSE
        }
      }
    } else {
      # A2 rerand, limit number of iterations 
      if (!pair_switch) {
        # A2.1 independent rerand
        if (cont) {
          n_jump <- n_iter <- max_iter
          for (iter in 2:max_iter) {
            Wt <- sample(1:n, n1)
            Mt <- compute_M_set(Wt, X, prop_cov_inv)
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) {
              n_jump <- n_iter <- iter
              break
            }
          }
        }
      } else {
        # A2.2 pair-switching rerand
        if (cont) {
          n_iter <- max_iter
          for (iter in 2:max_iter) {
            Wstar <- Wt
            i <- sample(Wt, 1)
            j <- sample(setdiff(1:n, Wt), 1)
            Wstar[which(Wstar == i)] <- j
            Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
            if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
              Wt <- Wstar
              Mt <- Mstar
              n_jump <- n_jump + 1
            }
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) {
              n_iter <- iter
              break
            }
          }
        }
      }
      # output the best rather than the last
      Wt <- Wbest
      Mt <- Mbest
    }
    # output n_iter & n_jump
    data.frame(n_jump = n_jump, n_iter = n_iter)
  }
}

rerand <- function(X, n1 = round(nrow(X) / 2), n_fisher = 1, 
                   p_a = 0.001, pair_switch = FALSE, temper = 1 / 10,
                   max_iter = Inf, parallel = FALSE, n_cores = 1, 
                   report_iter = FALSE) {
  n <- nrow(X)
  
  if (p_a == 1) {
    cov_inv <- NULL
    prop_cov_inv <- NULL
    H <- h <- NULL
  } else {
    cov_inv <- solve(cov(X))
    prop_cov_inv <- n1 * (1 - n1 / n) * cov_inv
    if (pair_switch) {
      H <- X %*% cov_inv %*% t(X) / (n1 * (1 - n1 / n))
      h <- 2 * (n1 / n) * rowSums(H)
    } else {
      H <- h <- NULL
    }
  }
  
  a <- qchisq(p = p_a, df = ncol(X))
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      rerand_one(X, n, n1, prop_cov_inv, H, h,
                 a, pair_switch, temper,
                 max_iter, report_iter)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      rerand_one(X, n, n1, prop_cov_inv, H, h,
                 a, pair_switch, temper,
                 max_iter, report_iter)
    })
  }
  if (n_fisher == 1) {
    unlist(out)
  } else {
    out
  }
}


# greedy pair switching ---------------------------------------------------

gps_one <- function(X, n, n1, prop_cov_inv, H, h, report_iter) {
  # start point
  Wt <- sample(1:n, n1)
  Mt <- compute_M_set(Wt, X, prop_cov_inv)
  n_jump <- 1
  cont <- TRUE
  
  # main ----
  if (!report_iter) {
    while (cont) {
      Wstar_all <- expand.grid(Wt, (1:n)[-Wt])
      Mstar_all <- apply(Wstar_all, 1, function(ij) {
        Wstar <- Wt
        Wstar[which(Wstar == ij[1])] <- ij[2]
        compute_Mstar_set(Wt, Wstar, ij[1], ij[2], Mt, H, h)
      })
      Mstar <- min(Mstar_all)
      ijstar <- Wstar_all[which.min(Mstar_all), ]
      if (Mstar < Mt) {
        Wt[which(Wt == ijstar[, 1])] <- ijstar[, 2]
        Mt <- Mstar
      } else {
        cont <- FALSE
      }
    }
    # output binary representation of local optimal Wt 
    Wout <- rep(0, n)
    Wout[Wt] <- 1
    return(Wout)
  } else {
    # A record iterations (for developer) ----
    while (cont) {
      Wstar_all <- expand.grid(Wt, (1:n)[-Wt])
      Mstar_all <- apply(Wstar_all, 1, function(ij) {
        Wstar <- Wt
        Wstar[which(Wstar == ij[1])] <- ij[2]
        compute_Mstar_set(Wt, Wstar, ij[1], ij[2], Mt, H, h)
      })
      Mstar <- min(Mstar_all)
      ijstar <- Wstar_all[which.min(Mstar_all), ]
      if (Mstar < Mt) {
        Wt[which(Wt == ijstar[, 1])] <- ijstar[, 2]
        Mt <- Mstar
        n_jump <- n_jump + 1
      } else {
        cont <- FALSE
      }
    }
    # output n_iter & n_jump
    n_iter <- n_jump * n1 * (n - n1)
    data.frame(n_jump = n_jump, n_iter = n_iter)
  }
}

gps <- function(X, n1 = round(nrow(X) / 2), n_fisher = 1, 
                parallel = FALSE, n_cores = 1, 
                report_iter = FALSE) {
  n <- nrow(X)
  cov_inv <- solve(cov(X))
  prop_cov_inv <- n1 * (1 - n1 / n) * cov_inv
  H <- X %*% cov_inv %*% t(X) / (n1 * (1 - n1 / n))
  h <- 2 * (n1 / n) * rowSums(H)
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      gps_one(X, n, n1, prop_cov_inv, H, h, report_iter)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      gps_one(X, n, n1, prop_cov_inv, H, h, report_iter)
    })
  }
  if (n_fisher == 1) {
    unlist(out)
  } else {
    out
  }
}


# # test
# # speed: pair-switching rerand > independent rerand > gps
# # number of iter: pair-switching rerand < independent rerand
# n <- 100
# p <- 10
# X <- matrix(rnorm(n * p), n, p)
# system.time(rerand(X, n_fisher = 10, pair_switch = F))
# system.time(rerand(X, n_fisher = 10, pair_switch = T))
# system.time(gps(X, n_fisher = 10))
# rerand(X, pair_switch = F, report_iter = T)
# rerand(X, pair_switch = T, report_iter = T)
# rerand(X, pair_switch = F, report_iter = T, max_iter = 10)
# rerand(X, pair_switch = F, report_iter = T, max_iter = 100)
# rerand(X, pair_switch = T, report_iter = T, max_iter = 10)
# rerand(X, pair_switch = T, report_iter = T, max_iter = 100)
# gps(X, report_iter = T)
# 
# bench::mark(
#   rerand(X, n_fisher = 10, pair_switch = T),
#   rerand(X, n_fisher = 10, pair_switch = T, report_iter = T),
#   check = F, max_iterations = 1
# )



# sequential rerandomization ---------------------------------------------------------
# W_fix: binary representation
# Wt_fix: set representation


rerand_part <- function(X_fix, X_wait, n1_wait, p_a, W_fix = NULL, Wt_fix = NULL, 
                        M_fix = NULL, prop_cov_inv_fix = NULL, 
                        prop_cov_inv = NULL, H = NULL, h = NULL,
                        pair_switch = FALSE, temper = 1 / 10, renew_iter = Inf, 
                        max_iter = Inf, 
                        output_set = FALSE, report_iter = FALSE) {
  if (is.null(Wt_fix) & !is.null(W_fix)) {
    Wt_fix <- which(W_fix == 1)
  }
  if (is.null(X_fix)) {
    X <- X_wait
    p <- ncol(X)
    n <- nrow(X)
    n_fix <- 0
    n_wait <- n
    n1 <- n1_wait
    M_fix <- 0
  } else {
    X <- rbind(X_fix, X_wait)
    p <- ncol(X)
    n <- nrow(X)
    n_fix <- nrow(X_fix)
    n_wait <- nrow(X_wait)
    n1 <- length(Wt_fix) + n1_wait
    if (is.null(M_fix)) {
      M_fix <- compute_M_set(Wt_fix, X_fix, prop_cov_inv_fix)
    }
  }
  ncp_chi <- n_fix / n_wait * M_fix
  a <- qchisq(p = p_a, df = p, ncp = ncp_chi) * n_wait / n
  
  if (is.null(prop_cov_inv)) {
    prop_cov_inv <- n1 * (1 - n1 / n) * solve(cov(X))
  }
  if (pair_switch) {
    if (is.null(H)) {
      H <- X %*% cov_inv %*% t(X) / (n1 * (1 - n1 / n))
    }
    if (is.null(h)) {
      h <- 2 * (n1 / n) * rowSums(H)
    }
  }
  
  # start point
  Wt_wait <- (n_fix + 1):n
  Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
  Mt <- compute_M_set(Wt, X, prop_cov_inv)
  n_iter <- 1
  n_jump <- 1
  if (Mt <= a) {
    cont <- FALSE
    Wbest <- Wt
    Mbest <- Mt
  } else {
    cont <- TRUE
    Wbest <- Wt
    Mbest <- Mt
  }
  
  # main ----
  if (!report_iter) {
    if (max_iter == Inf) {
      # 1 rerand, not limit number of iterations
      if (!pair_switch) {
        # 1.1 independent rerand
        while (cont) {
          Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
          Mt <- compute_M_set(Wt, X, prop_cov_inv)
          if (Mt <= a) cont <- FALSE
        }
      } else if (renew_iter == Inf) {
        # 1.2 pair-switching rerand
        while (cont) {
          Wstar <- Wt
          i <- sample(intersect(Wt_wait, Wt), 1)
          j <- sample(setdiff(Wt_wait, Wt), 1)
          Wstar[which(Wstar == i)] <- j
          Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
          if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
            Wt <- Wstar
            Mt <- Mstar
          }
          if (Mt <= a) cont <- FALSE
        }
      } else {
        # 1.3 renewable pair-switching rerand
        while (cont) {
          if (n_iter %% renew_iter != 0) {
            Wstar <- Wt
            i <- sample(intersect(Wt_wait, Wt), 1)
            j <- sample(setdiff(Wt_wait, Wt), 1)
            Wstar[which(Wstar == i)] <- j
            Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
          } else {
            Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
            Mt <- compute_M_set(Wt, X, prop_cov_inv)
          }
          n_iter <- n_iter + 1
          if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
            Wt <- Wstar
            Mt <- Mstar
          }
          if (Mt <= a) cont <- FALSE
        }
      }
    } else {
      # 2 rerand, limit number of iterations
      if (!pair_switch) {
        # 2.1 independent rerand
        if (cont) {
          for (iter in 2:max_iter) {
            Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
            Mt <- compute_M_set(Wt, X, prop_cov_inv)
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) break
          }
        }
      } else if (renew_iter == Inf) {
        # 2.2 pair-switching rerand
        if (cont) {
          for (iter in 2:max_iter) {
            Wstar <- Wt
            i <- sample(intersect(Wt_wait, Wt), 1)
            j <- sample(setdiff(Wt_wait, Wt), 1)
            Wstar[which(Wstar == i)] <- j
            Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
            if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
              Wt <- Wstar
              Mt <- Mstar
            }
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) break
          }
        }
      } else {
        # 2.3 renewable pair-switching rerand
        if (cont) {
          for (iter in 2:max_iter) {
            if (iter %% renew_iter != 0) {
              Wstar <- Wt
              i <- sample(intersect(Wt_wait, Wt), 1)
              j <- sample(setdiff(Wt_wait, Wt), 1)
              Wstar[which(Wstar == i)] <- j
              Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
            } else {
              Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
              Mt <- compute_M_set(Wt, X, prop_cov_inv)
            }
            if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
              Wt <- Wstar
              Mt <- Mstar
            }
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) break
          }
        }
      }
      # output the best rather than the last
      Wt <- Wbest
      Mt <- Mbest
    }
    # output set/binary representation of Wt
    if (output_set) {
      return(list(W = Wt, M = Mt))
    } else {
      Wout <- rep(0, n)
      Wout[Wt] <- 1
      return(list(W = Wout, M = Mt))
    }
  } else {
    # A record iterations (for developer) ----
    if (max_iter == Inf) {
      # A1 rerand, not limit iterations
      if (!pair_switch) {
        # A1.1 independent rerand
        while (cont) {
          Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
          Mt <- compute_M_set(Wt, X, prop_cov_inv)
          n_iter <- n_iter + 1
          if (Mt <= a) cont <- FALSE
        }
        n_jump <- n_iter
      } else if (renew_iter == Inf) {
        # A1.2 pair-switching rerand
        while (cont) {
          Wstar <- Wt
          i <- sample(intersect(Wt_wait, Wt), 1)
          j <- sample(setdiff(Wt_wait, Wt), 1)
          Wstar[which(Wstar == i)] <- j
          Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
          n_iter <- n_iter + 1
          if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
            Wt <- Wstar
            Mt <- Mstar
            n_jump <- n_jump + 1
          }
          if (Mt <= a) cont <- FALSE
        }  
      } else if (renew_iter != Inf) {
        # A1.3 renewable pair-switching rerand
        while (cont) {
          if (n_iter %% renew_iter != 0) {
            Wstar <- Wt
            i <- sample(intersect(Wt_wait, Wt), 1)
            j <- sample(setdiff(Wt_wait, Wt), 1)
            Wstar[which(Wstar == i)] <- j
            Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
          } else {
            Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
            Mt <- compute_M_set(Wt, X, prop_cov_inv)
          }
          n_iter <- n_iter + 1
          if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
            Wt <- Wstar
            Mt <- Mstar
            n_jump <- n_jump + 1
          }
          if (Mt <= a) cont <- FALSE
        }
      }
    } else {
      # A2 rerand, limit iterations
      if (!pair_switch) {
        # A2.1 independent rerand
        if (cont) {
          n_jump <- n_iter <- max_iter
          for (iter in 2:max_iter) {
            Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
            Mt <- compute_M_set(Wt, X, prop_cov_inv)
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) {
              n_jump <- n_iter <- iter
              break
            }
          }
        }
      } else if (renew_iter == Inf) {
        # A2.2 pair-switching rerand
        if (cont) {
          n_iter <- max_iter
          for (iter in 2:max_iter) {
            Wstar <- Wt
            i <- sample(intersect(Wt_wait, Wt), 1)
            j <- sample(setdiff(Wt_wait, Wt), 1)
            Wstar[which(Wstar == i)] <- j
            Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
            if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
              Wt <- Wstar
              Mt <- Mstar
              n_jump <- n_jump + 1
            }
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) {
              n_iter <- iter
              break
            }
          }
        }
      } else if (renew_iter != Inf) {
        # A2.3 renewable pair-switching rerand
        if (cont) {
          n_iter <- max_iter
          for (iter in 2:max_iter) {
            if (iter %% renew_iter != 0) {
              Wstar <- Wt
              i <- sample(intersect(Wt_wait, Wt), 1)
              j <- sample(setdiff(Wt_wait, Wt), 1)
              Wstar[which(Wstar == i)] <- j
              Mstar <- compute_Mstar_set(Wt, Wstar, i, j, Mt, H, h)
            } else {
              Wt <- c(Wt_fix, sample(Wt_wait, n1_wait))
              Mt <- compute_M_set(Wt, X, prop_cov_inv)
            }
            if (rbinom(1, 1, min((Mt / Mstar)^(1 / temper), 1)) == 1) {
              Wt <- Wstar
              Mt <- Mstar
              n_jump <- n_jump + 1
            }
            if (Mt < Mbest) {
              Wbest <- Wt
              Mbest <- Mt
            }
            if (Mt <= a) {
              n_iter <- iter
              break
            }
          }
        }
      }
      # output the best rather than the last
      Wt <- Wbest
      Mt <- Mbest
    }
    # output W, M, n_iter, n_jump
    if (output_set) {
      return(list(W = Wt, M = Mt, n_jump = n_jump, n_iter = n_iter))
    } else {
      Wout <- rep(0, n)
      Wout[Wt] <- 1
      return(list(W = Wout, M = Mt, n_jump = n_jump, n_iter = n_iter))
    }
  }
}


seq_rerand_one <- function(X, n_cum, n1_seq, p_a_seq,
                           prop_cov_inv_seq, H_seq, h_seq,
                           pair_switch, temper, renew_iter_seq,
                           max_iter_seq, report_iter) {
  WM_k <- list()
  WM_k$W <- NULL
  WM_k$M <- NULL
  
  if (report_iter) {
    iter_list <- list()
  }
  
  if (all(p_a_seq==1)) {
    Wt <- NULL
    
    for (k in 1:length(n_cum)) {
      if (k == 1) {
        Wt_wait <- 1:n_cum[k]
      } else {
        Wt_wait <- (n_cum[k-1] + 1):n_cum[k]
      }
      Wt <- c(Wt, sample(Wt_wait, n1_seq[k]))
      
      if (report_iter) {
        iter_list[[k]] <- data.frame(k = k, n_jump = 1, n_iter = 1)
      }
    }
    
    WM_k <- list(W = Wt)
  } else {
    for (k in 1:length(n_cum)) {
      if (k == 1) {
        X_fix_k <- NULL
        X_wait_k <- X[1:n_cum[k],]
      } else {
        X_fix_k <- X[1:n_cum[k - 1],]
        X_wait_k <- X[(n_cum[k - 1] + 1):n_cum[k],]
      }
      if (p_a_seq[k] == 1) {
        Wt_fix <- WM_k$W
        if (k == 1) {
          Wt_wait <- 1:n_cum[k]
        } else {
          Wt_wait <- (n_cum[k-1] + 1):n_cum[k]
        }
        Wt <- c(Wt_fix, sample(Wt_wait, n1_seq[k]))
        Mt <- compute_M_set(Wt, X[1:n_cum[k],], prop_cov_inv_seq[[k]])
        WM_k <- list(W = Wt, M = Mt)
        if (report_iter) {
          iter_list[[k]] <- data.frame(k = k, n_jump = 1, n_iter = 1)
        }
      } else {
        WM_k <- rerand_part(X_fix = X_fix_k, X_wait = X_wait_k, n1_wait = n1_seq[k],
                            p_a = p_a_seq[k], Wt_fix = WM_k$W, M_fix = WM_k$M,
                            prop_cov_inv = prop_cov_inv_seq[[k]], 
                            H = H_seq[[k]], h = h_seq[[k]],
                            pair_switch = pair_switch, temper = temper, 
                            renew_iter = renew_iter_seq[k],
                            max_iter = max_iter_seq[k],
                            output_set = T, report_iter = report_iter)
        if (report_iter) {
          iter_list[[k]] <- data.frame(k = k, n_jump = WM_k$n_jump, n_iter = WM_k$n_iter)
        }
      }
    }
  }
  
  Wout <- rep(0, nrow(X))
  Wout[WM_k$W] <- 1
  if (report_iter) {
    list(iter_list = iter_list, W = Wout)
  } else {
    Wout
  }
}


seq_rerand <- function(X, n_seq, p_a_seq, n1_seq = round(n_seq / 2), n_fisher = 1, 
                       pair_switch = FALSE, temper = 1 / 10,  
                       renew_iter_seq = rep(Inf, length(n_seq)),
                       max_iter_seq = rep(Inf, length(n_seq)),
                       parallel = FALSE, n_cores = 1, report_iter = FALSE
) {
  n_cum <- cumsum(n_seq)
  n1_cum <- cumsum(n1_seq)
  
  if (all(p_a_seq==1)) {
    prop_cov_inv_seq <- NULL
    H_seq <- NULL
    h_seq <- NULL
  } else {
    pre_cum <- lapply(1:length(n_seq), function(k) {
      X_cum_k <- X[1:n_cum[k],]
      cov_inv <- solve(cov(X_cum_k))
      prop_cov_inv <- n1_cum[k] * (1 - n1_cum[k] / n_cum[k]) * cov_inv
      if (pair_switch) {
        H <- X_cum_k %*% cov_inv %*% t(X_cum_k) / (n1_cum[k] * (1 - n1_cum[k] / n_cum[k]))
        h <- 2 * (n1_cum[k] / n_cum[k]) * rowSums(H)
      } else {
        H <- h <- NULL
      }
      list(prop_cov_inv = prop_cov_inv, H = H, h = h)
    })
    prop_cov_inv_seq <- lapply(pre_cum, function(x) x$prop_cov_inv)
    H_seq <- lapply(pre_cum, function(x) x$H)
    h_seq <- lapply(pre_cum, function(x) x$h)
  }
  
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      seq_rerand_one(X, n_cum, n1_seq, p_a_seq,
                     prop_cov_inv_seq, H_seq, h_seq,
                     pair_switch, temper, renew_iter_seq, 
                     max_iter_seq, report_iter)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      seq_rerand_one(X, n_cum, n1_seq, p_a_seq,
                     prop_cov_inv_seq, H_seq, h_seq,
                     pair_switch, temper, renew_iter_seq,
                     max_iter_seq, report_iter)
    })
  }
  if (n_fisher == 1) {
    out[[1]]
  } else {
    out
  }
}


# # test
# # speed: pair-switching rerand > independent rerand
# n_seq <- c(20, 20, 20)
# n <- sum(n_seq)
# p <- 10
# X <- matrix(rnorm(n * p), n, p)
# p_a_seq <- 1 / c(6, 19, 75)
# system.time(seq_rerand(X, n_seq, p_a_seq = rep(1, length(n_seq)), n_fisher = 100)) # SeqCR
# system.time(seq_rerand(X, n_seq, p_a_seq, n_fisher = 100, pair_switch = F)) # SeqRR
# system.time(seq_rerand(X, n_seq, p_a_seq, n_fisher = 100, pair_switch = T)) # SeqPSRR
# 
# bench::mark(
#   seq_rerand(X, n_seq, p_a_seq, n_fisher = 10, pair_switch = F),
#   seq_rerand(X, n_seq, p_a_seq, n_fisher = 10, pair_switch = T),
#   check = F, max_iterations = 1
# )
# 
# bench::mark(
#   seq_rerand(X, n_seq, p_a_seq, n_fisher = 10, pair_switch = T),
#   seq_rerand(X, n_seq, p_a_seq, n_fisher = 10, pair_switch = T, report_iter = T),
#   check = F, max_iterations = 1
# )
# 
# seq_rerand(X, n_seq, p_a_seq = rep(1, length(n_seq)), n_fisher = 1, report_iter = T) # SeqCR
# seq_rerand(X, n_seq, p_a_seq, n_fisher = 1, pair_switch = F, report_iter = T) # SeqRR
# seq_rerand(X, n_seq, p_a_seq, n_fisher = 1, pair_switch = T, report_iter = T) # SeqPSRR
# 
# # rerand_part
# n_cum <- cumsum(n_seq)
# n1_seq <- round(n_seq / 2)
# n1_cum <- cumsum(n1_seq)
# pre_cum <- lapply(1:length(n_seq), function(k) {
#   X_cum_k <- X[1:n_cum[k],]
#   cov_inv <- solve(cov(X_cum_k))
#   prop_cov_inv <- n1_cum[k] * (1 - n1_cum[k] / n_cum[k]) * cov_inv
#   H <- X_cum_k %*% cov_inv %*% t(X_cum_k) / (n1_cum[k] * (1 - n1_cum[k] / n_cum[k]))
#   h <- 2 * (n1_cum[k] / n_cum[k]) * rowSums(H)
#   list(prop_cov_inv = prop_cov_inv, H = H, h = h)
# })
# prop_cov_inv_seq <- lapply(pre_cum, function(x) x$prop_cov_inv)
# H_seq <- lapply(pre_cum, function(x) x$H)
# h_seq <- lapply(pre_cum, function(x) x$h)
# pair_switch <- T
# temper <- 0.1
# max_iter <- Inf
# report_iter <- F
# WM_k <- list()
# WM_k$W <- NULL
# WM_k$M <- NULL
# iter_list <- list()
# for (k in 1:(length(n_seq) - 1)) {
#   if (k == 1) {
#     X_fix_k <- NULL
#     X_wait_k <- X[1:n_cum[k],]
#   } else {
#     X_fix_k <- X[1:n_cum[k - 1],]
#     X_wait_k <- X[(n_cum[k - 1] + 1):n_cum[k],]
#   }
#   WM_k <- rerand_part(X_fix = X_fix_k, X_wait = X_wait_k, n1_wait = n1_seq[k],
#                       p_a = p_a_seq[k], Wt_fix = WM_k$W, M_fix = WM_k$M,
#                       prop_cov_inv = prop_cov_inv_seq[[k]],
#                       H = H_seq[[k]], h = h_seq[[k]],
#                       pair_switch = pair_switch, temper = temper,
#                       max_iter = max_iter,
#                       output_set = T, report_iter = report_iter)
#   if (report_iter) {
#     iter_list[[k]] <- data.frame(k = k, n_jump = WM_k$n_jump, n_iter = WM_k$n_iter)
#   }
# }
# k <- length(n_seq)
# X_fix_k <- X[1:n_cum[k - 1],]
# X_wait_k <- X[(n_cum[k - 1] + 1):n_cum[k],]
# bench::mark(
#   rerand_part(X_fix = X_fix_k, X_wait = X_wait_k, n1_wait = n1_seq[k],
#               p_a = p_a_seq[k], Wt_fix = WM_k$W, M_fix = WM_k$M,
#               prop_cov_inv = prop_cov_inv_seq[[k]],
#               H = H_seq[[k]], h = h_seq[[k]],
#               pair_switch = F, temper = temper,
#               max_iter = max_iter,
#               output_set = T, report_iter = report_iter),
#   rerand_part(X_fix = X_fix_k, X_wait = X_wait_k, n1_wait = n1_seq[k],
#               p_a = p_a_seq[k], Wt_fix = WM_k$W, M_fix = WM_k$M,
#               prop_cov_inv = prop_cov_inv_seq[[k]],
#               H = H_seq[[k]], h = h_seq[[k]],
#               pair_switch = T, temper = temper,
#               max_iter = max_iter,
#               output_set = T, report_iter = report_iter),
#   check = F
# )[, -1]


allocate_S <- function(S, K, p, n_seq = NULL, budget = "increase", s_min = 10) {
  if (budget == "increase") {
    if (is.null(n_seq)) {
      n_seq <- rep(1, K)
    } 
    Cp <- 2 * p / (p + 2) * (gamma(p / 2 + 1))^(2 / p)
    # upper bound of sK
    sK_u <- S
    # determine lower bound of sK
    s <- rep(sK_u, K)
    for (k in K:2) {
      s[k - 1] <- ((Cp * n_seq[k - 1]) / (p * n_seq[k]) * s[k])^(p / (p + 2))
    }
    s <- ifelse(s < s_min, s_min, round(s))
    sK_l <- S - sum(s[-K])
    # find sK such that sum(s) = S
    while(sK_u - sK_l > 1) {
      sK_m <- round((sK_u + sK_l) / 2)
      s <- rep(sK_m, K)
      for (k in K:2) {
        s[k - 1] <- ((Cp * n_seq[k - 1]) / (p * n_seq[k]) * s[k])^(p / (p + 2))
      }
      s <- ifelse(s < s_min, s_min, round(s))
      if (sum(s) - S > 0) {
        sK_u <- sK_m
      } else if (sum(s) - S < 0) {
        sK_l <- sK_m
      } else {
        break
      }
    }
    s[K] <- s[K] - sum(s) + S
  } else if (budget == "equal") {
    s <- rep(floor(S / K), K)
    S_rest <- S - sum(s)
    s <- s + c(rep(0, K - S_rest), rep(1, S_rest))
  }
  s
}

# allocate_S(1000, 15, 10)
# allocate_S(1000, 15, 10, budget = "equal")


# measurement of randomness -----------------------------------------------

compute_random <- function(W_ref) {
  n <- length(W_ref[[1]])
  n1 <- sum(W_ref[[1]])
  ind_mat <- combn(n, 2)
  pair_in_same_group <- sapply(W_ref, function(W) {
    apply(ind_mat, 2, function(ind) {
      abs(sum(W[ind]) - 1) # (0,0) or (1,1) -> 1; (1,0) or (0,1) -> 0
    })
  })
  ps <- rowMeans(pair_in_same_group)
  
  # En
  sn <- mean(ps)
  En_vec <- ifelse(
    ps == 0 | ps == 1, 
    0,
    ps * log(ps) + (1 - ps) * log(1 - ps)
  )
  En <- mean(En_vec) / (sn * log(sn) + (1 - sn) * log(1 - sn))
  
  # Dn
  ps_det <- pair_in_same_group[,1]
  Dn <- sd(ps) / sd(ps_det)
  
  # lambda_max
  W_ref_alt <- sapply(W_ref, function(W) {
    2 * W - 1
  })
  Sigma_W <- var(t(W_ref_alt))
  Ln <- eigen(Sigma_W, symmetric = T, only.values = T)$values[1]
  
  # output
  data.frame(En = En, Dn = Dn, Ln = Ln)
}


# W_ref <- lapply(1:1000, function(x) {
#   sample(c(rep(1, 20), rep(0, 20)))
# })
# compute_random(W_ref)


