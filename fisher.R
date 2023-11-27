# Fisher randomization tests ----------------------------------------------

comp_stat <- function(W, Y) {
  mean(Y[W == 1]) - mean(Y[W == 0])
}

comp_stat_set <- function(Wt, Y) {
  mean(Y[Wt]) - mean(Y[-Wt])
}

fisher_p <- function(W, Y, W_ref = NULL, tau0 = 0, H1 = "!=",
                     Wt_ref = NULL, W_is_list = TRUE) {
  stat_obs <- comp_stat(W, Y)
  Y1 <- ifelse(W == 1, Y, Y + tau0)
  Y0 <- ifelse(W == 0, Y, Y - tau0)
  if (W_is_list) {
    if (is.null(Wt_ref)) {
      stat_star <- sapply(W_ref, function(Wstar) {
        Ystar <- Y1 * Wstar + Y0 * (1 - Wstar)
        comp_stat(Wstar, Ystar)
      })
    } else {
      stat_star <- sapply(Wt_ref, function(Wstar) {
        Ystar <- Y0
        Ystar[Wstar] <- Y1[Wstar]
        comp_stat_set(Wstar, Ystar)
      })
    }
  } else {
    if (is.null(Wt_ref)) {
      stat_star <- apply(W_ref, 1, function(Wstar) {
        Ystar <- Y1 * Wstar + Y0 * (1 - Wstar)
        comp_stat(Wstar, Ystar)
      })
    } else {
      stat_star <- apply(Wt_ref, 1, function(Wstar) {
        Ystar <- Y0
        Ystar[Wstar] <- Y1[Wstar]
        comp_stat_set(Wstar, Ystar)
      })
    }
  }
  if (H1 == ">") {
    sum(stat_star >= stat_obs) / length(stat_star)
  } else if (H1 == "<") {
    sum(stat_star <= stat_obs) / length(stat_star)
  } else if (H1 == "!=") {
    2 * min(
      sum(stat_star >= stat_obs) / length(stat_star),
      sum(stat_star <= stat_obs) / length(stat_star)
    )
  }
}


# # test
# n <- 100
# Y1 <- round(rnorm(n), 1)
# Y0 <- Y1
# W_ref <- lapply(1:1000, function(x) {
#   sample(c(rep(1, n / 2), rep(0, n - n / 2)))
# })
# Wt_ref <- lapply(W_ref, function(x) {
#   which(x == 1)
# })
# W_ref_mat <- t(sapply(W_ref, function(x) {
#   x
# }))
# Wt_ref_mat <- t(sapply(W_ref, function(x) {
#   which(x == 1)
# }))
# W <- W_ref[[1]]
# Y <- Y1 * W + Y0 * (1 - W)
# bench::mark(
#   fisher_p(W, Y, W_ref), # list + binary representation
#   fisher_p(W, Y, Wt_ref = Wt_ref), # list + set representation
#   fisher_p(W, Y, W_ref_mat, W_is_list = F), # matrix + binary representation
#   fisher_p(W, Y, Wt_ref = Wt_ref_mat, W_is_list = F) # matrix + set representation
# )
# # conclusion: set representation is faster


# Fisher interval estimation -----------------------------------------------------


# 1 exact approach (the proposed method) ----

solve_tau <- function(W, Wstar, Y) {
  if (identical(Wstar, W)) {
    NA
  } else {
    a <- Y[W == 1 & Wstar == 0]
    b <- Y[W == 0 & Wstar == 1]
    (sum(a) - sum(b)) / length(a)
  }
}

fisher_ci <- function(W, Y, W_ref, alpha = 0.05, H1 = "!=") {
  tau_solve <- sapply(W_ref, function(Wstar) {
    solve_tau(W, Wstar, Y)
  })
  if (H1 == ">") {
    ci_up <- Inf
    if (1 / length(W_ref) >= alpha) {
      ci_low <- -Inf
    } else {
      tau_solve_l <- ifelse(is.na(tau_solve), min(tau_solve, na.rm = T), tau_solve)
      ci_low <- sort(tau_solve_l)[floor(length(W_ref) * alpha) + 1]
    }
  } else if (H1 == "<") {
    ci_low <- -Inf
    if (1 / length(W_ref) >= alpha) {
      ci_up <- Inf
    } else {
      tau_solve_h <- ifelse(is.na(tau_solve), max(tau_solve, na.rm = T), tau_solve)
      ci_up <- sort(tau_solve_h, decreasing = T)[floor(length(W_ref) * alpha) + 1]
    }
  } else if (H1 == "!=") {
    if (length(W_ref) * alpha / 2 <= 1) {
      ci_low <- -Inf
      ci_up <- Inf
    } else {
      tau_solve_l <- ifelse(is.na(tau_solve), min(tau_solve, na.rm = T), tau_solve)
      ci_low <- sort(tau_solve_l)[floor(length(W_ref) * alpha / 2) + 1]
      tau_solve_h <- ifelse(is.na(tau_solve), max(tau_solve, na.rm = T), tau_solve)
      ci_up <- sort(tau_solve_h, decreasing = T)[floor(length(W_ref) * alpha / 2) + 1]
    }
  }
  data.frame(ci_low, ci_up)
}


# 2 bisection ----

fisher_ci_bi_one <- function(W, Y, W_ref, alpha = 0.05, H1 = ">", 
                             theta_l, theta_r, tol, max_iter) {
  Wt_ref <- lapply(W_ref, function(x) {
    which(x == 1)
  })
  if (1 / length(W_ref) >= alpha) {
    ci_up <- Inf
    ci_low <- -Inf
  } else {
    p_l <- fisher_p(W, Y, Wt_ref = Wt_ref, tau0 = theta_l, H1 = H1)
    p_r <- fisher_p(W, Y, Wt_ref = Wt_ref, tau0 = theta_r, H1 = H1)
    if ((alpha - p_l) * (p_r - alpha) < 0) {
      stop(paste0("p_l = ", p_l, ", p_r = ", p_r, ", reset theta_l and theta_r"))
    } else {
      if (H1 == ">") {
        for (i in 1:max_iter) {
          theta_m <- (theta_r + theta_l) / 2
          p_m <- fisher_p(W, Y, Wt_ref = Wt_ref, tau0 = theta_m, H1 = H1)
          if (p_m <= alpha) {
            theta_l <- theta_m
          } else {
            theta_r <- theta_m
          }
          if (theta_r - theta_l <= tol) break
        }
        ci_low <- theta_l
        ci_up <- Inf
      } else if (H1 == "<") {
        for (i in 1:max_iter) {
          theta_m <- (theta_r + theta_l) / 2
          p_m <- fisher_p(W, Y, Wt_ref = Wt_ref, tau0 = theta_m, H1 = H1)
          if (p_m > alpha) {
            theta_l <- theta_m
          } else {
            theta_r <- theta_m
          }
          if (theta_r - theta_l <= tol) break
        }
        ci_low <- -Inf
        ci_up <- theta_r
      }
    }
  }
  data.frame(ci_low, ci_up)
}

fisher_ci_bi <- function(W, Y, W_ref, alpha = 0.05, H1 = "!=", 
                         theta_l = -100, theta_r = 100, tol = 0.001,
                         max_iter = NULL) {
  if (is.null(max_iter)) {
    max_iter <- ceiling(log2((theta_r - theta_l) / tol))
  }
  if (H1 != "!=") {
    fisher_ci_bi_one(W, Y, W_ref, alpha, H1, theta_l, theta_r, tol, max_iter)
  } else if (H1 == "!=") {
    ci_low <- fisher_ci_bi_one(W, Y, W_ref, alpha / 2, ">", 
                               theta_l, theta_r, tol, max_iter)$ci_low
    ci_up <-  fisher_ci_bi_one(W, Y, W_ref, alpha / 2, "<", 
                               theta_l, theta_r, tol, max_iter)$ci_up
    data.frame(ci_low, ci_up)
  }
}

# 3 Robbins-Monro ----

solve_start <- function(W, Y, W_ref) {
  stat_obs <- comp_stat(W, Y)
  Y1 <- ifelse(W == 1, Y, Y + stat_obs)
  Y0 <- ifelse(W == 0, Y, Y - stat_obs)
  stat_star <- sort(sapply(W_ref, function(Wstar) {
    Ystar <- Y1 * Wstar + Y0 * (1 - Wstar)
    comp_stat(Wstar, Ystar)
  }))
  t1 <- stat_star[2]
  t2 <- rev(stat_star)[2]
  data.frame(ci_low = stat_obs - (t2 - t1) / 2, ci_up = stat_obs + (t2 - t1) / 2)
}

fisher_ci_rm_one <- function(W, Y, W_ref, alpha = 0.05, H1 = ">", ci_tmp, k) {
  z_a <- qnorm(alpha, lower.tail = F)
  g <- (z_a * (2 * pi)^(-0.5) * exp(- z_a^2 / 2))
  stat_obs <- comp_stat(W, Y)
  if (1 / length(W_ref) >= alpha) {
    ci_up <- Inf
    ci_low <- -Inf
  } else {
    if (H1 == ">") {
      ci_up <- Inf
      len_step <- k * abs(ci_tmp - stat_obs) / g
      for (i in 1:length(W_ref)) {
        Y1 <- ifelse(W == 1, Y, Y + ci_tmp)
        Y0 <- ifelse(W == 0, Y, Y - ci_tmp)
        Wstar <- W_ref[[i]]
        Ystar <- Y1 * Wstar + Y0 * (1 - Wstar)
        stat_star <- comp_stat(Wstar, Ystar)
        if (stat_star < stat_obs) {
          # indicate p(ci_tmp) < alpha
          ci_tmp <- ci_tmp + len_step * alpha / i
        } else {
          # indicate p(ci_tmp) > alpha
          ci_tmp <- ci_tmp - len_step * (1 - alpha) / i
        }
      }
      ci_low <- ci_tmp
    } else if (H1 == "<") {
      ci_low <- -Inf
      len_step <- k * abs(ci_tmp - stat_obs) / g
      for (i in 1:length(W_ref)) {
        Y1 <- ifelse(W == 1, Y, Y + ci_tmp)
        Y0 <- ifelse(W == 0, Y, Y - ci_tmp)
        Wstar <- W_ref[[i]]
        Ystar <- Y1 * Wstar + Y0 * (1 - Wstar)
        stat_star <- comp_stat(Wstar, Ystar)
        if (stat_star > stat_obs) {
          # indicate p(ci_tmp) < alpha
          ci_tmp <- ci_tmp - len_step * alpha / i
        } else {
          # indicate p(ci_tmp) > alpha
          ci_tmp <- ci_tmp + len_step * (1 - alpha) / i
        }
      }
      ci_up <- ci_tmp
    }
  }
  data.frame(ci_low, ci_up)
}

fisher_ci_rm <- function(W, Y, W_ref, alpha = 0.05, H1 = "!=", k = 2) {
  if (H1 == ">") {
    fisher_ci_rm_one(W, Y, W_ref, alpha, H1, 
                     ci_tmp = solve_start(W, Y, W_ref)$ci_low, k = k)
  } else if (H1 == "<") {
    fisher_ci_rm_one(W, Y, W_ref, alpha, H1, 
                     ci_tmp = solve_start(W, Y, W_ref)$ci_up, k = k)
  } else if (H1 == "!=") {
    ci_start <- solve_start(W, Y, W_ref)
    ci_low <- fisher_ci_rm_one(W, Y, W_ref, alpha / 2, H1 = ">", 
                               ci_tmp = ci_start$ci_low, k = k)$ci_low
    ci_up  <- fisher_ci_rm_one(W, Y, W_ref, alpha / 2, H1 = "<", 
                               ci_tmp = ci_start$ci_up, k = k)$ci_up
    data.frame(ci_low, ci_up)
  }
}



# # test
# # potential outcomes
# n <- 100
# Y1 <- round(rnorm(n), 1)
# Y0 <- Y1
# # given true tau, the distribution of p
# W_ref <- lapply(1:10, function(x) {
#   sample(c(rep(1, n / 2), rep(0, n - n / 2)))
# })
# p_two_sided <- sapply(W_ref, function(Wstar) {
#   Y <- Y1 * Wstar + Y0 * (1 - Wstar)
#   fisher_p_two(Wstar, Y, W_ref)
# })
# plot(ecdf(p_two_sided), verticals = T)
# abline(0, 1, col = "blue")
# 
# # given one assignment, p as a function of tau
# W_ref <- lapply(1:100, function(x) {
#   sample(c(rep(1, n / 2), rep(0, n - n / 2)))
# })
# W <- W_ref[[1]]
# Y <- Y1 * W + Y0 * (1 - W)
# tau_vec <- seq(-1, 1, length.out = 100)
# p_two_sided_vec <- unlist(
#   parallel::mclapply(tau_vec, function(tau0) {
#     fisher_p_two(W, Y, W_ref, tau0)
#   }, mc.cores = 4)
# )
# tau_est <- comp_stat(W, Y)
# tau_low <- max(tau_vec[tau_vec <= tau_est & p_two_sided_vec <= 0.05])
# tau_up <- min(tau_vec[tau_vec >= tau_est & p_two_sided_vec <= 0.05])
# plot(tau_vec, p_two_sided_vec)
# abline(v = tau_est, col = "red")
# abline(v = tau_low, col = "blue")
# abline(v = tau_up, col = "blue")
# 
# # solve tau0 such that comp_stat(W, Y) = comp_stat(Wstar, Y_imp)
# W <- W_ref[[1]]
# Wstar <- W_ref[[2]]
# tau0 <- solve_tau(W, Wstar, Y)
# Y1 <- ifelse(W == 1, Y, Y + tau0)
# Y0 <- ifelse(W == 0, Y, Y - tau0)
# Y_imp <- Y1 * Wstar + Y0 * (1 - Wstar)
# comp_stat(W, Y) - comp_stat(Wstar, Y_imp)
# 
# # empirical dist function of tau_solve = p value function
# tau_solve <- sapply(W_ref, function(Wstar) {
#   solve_tau(W, Wstar, Y)
# })
# tau_solve <- ifelse(is.na(tau_solve), min(tau_solve, na.rm = T), tau_solve)
# tau_vec <- tau_solve + 10e-10
# p_one_vec <- unlist(
#   parallel::mclapply(tau_vec, function(tau0) {
#     fisher_p_one(W, Y, W_ref, tau0)
#   }, mc.cores = 4)
# )
# plot(ecdf(tau_solve))
# segments(min(tau_solve) - 10^10, 1 / length(W_ref), min(tau_solve), 1 / length(W_ref))
# points(tau_vec, p_one_vec, col = "red")
# identical(ecdf(tau_solve)(tau_vec), p_one_vec)




# # check CI
# n <- 100
# Y1 <- round(rnorm(n), 1)
# Y0 <- Y1
# W_ref <- lapply(1:10000, function(x) {
#   sample(c(rep(1, n / 2), rep(0, n - n / 2)))
# })
# W <- W_ref[[1]]
# Y <- Y1 * W + Y0 * (1 - W)
# err <- 1e-14
# 
# library(tictoc)
# tic()
# ci <- fisher_ci(W, Y, W_ref)
# toc()
# tic()
# ci_bi <- fisher_ci_bi(W, Y, W_ref)
# toc()
# 
# ci_bi$ci_low < ci$ci_low
# ci_bi$ci_up > ci$ci_up
# 
# ci_low <- ci$ci_low
# ci_up <- ci$ci_up
# fisher_p(W, Y, W_ref, ci_low + err) # > 0.05
# fisher_p(W, Y, W_ref, ci_low - err) # <= 0.05
# fisher_p(W, Y, W_ref, ci_up - err) # > 0.05
# fisher_p(W, Y, W_ref, ci_up + err) # <= 0.05
# 
# ci_low <- ci_bi$ci_low
# ci_up <- ci_bi$ci_up
# fisher_p(W, Y, W_ref, ci_low + err) # <= 0.05
# fisher_p(W, Y, W_ref, ci_low - err) # <= 0.05
# fisher_p(W, Y, W_ref, ci_up - err) # <= 0.05
# fisher_p(W, Y, W_ref, ci_up + err) # <= 0.05
# 
# # fisher_ci_rm(W, Y, W_ref, k = 4)
# # fisher_ci_rm(W, Y, W_ref, k = 2)
# # fisher_ci_rm(W, Y, W_ref, k = 1)
# # fisher_ci_rm(W, Y, W_ref, k = 0.5)


