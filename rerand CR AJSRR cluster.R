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

get_id = function(W, clust_size, n1_clust) {
  rep(clust_size*W, each = clust_size) - rep((clust_size-1):0, n1_clust)
}

# rerandomization ---------------------------------------------------------

ClustAJSRR_one <- function(X, n_clust, n1_clust, prop_cov_inv, a, max_iter) {
  n = nrow(X)
  clust_size = n / n_clust

  # start point
  Wt <- sample(1:n_clust, n1_clust)
  if (a == Inf) {
    cont <- FALSE
    Mt <- NA
    Wbest <- Wt
    Mbest <- Mt
  } else {
    Mt <- compute_M_set(get_id(Wt,clust_size,n1_clust), X, prop_cov_inv)
    Wbest <- Wt
    Mbest <- Mt
    if (Mt <= a) {
      cont <- FALSE
    } else {
      cont <- TRUE
    }
  }
  
  # main ----
  
  if (max_iter == Inf) {
    # 1 rerand, not limit number of iterations
    # 1.1 independent rerand
    while (cont) {
      Wt <- sample(1:n_clust, n1_clust)
      Mt <- compute_M_set(get_id(Wt,clust_size,n1_clust), X, prop_cov_inv)
      if (Mt <= a) cont <- FALSE
    }
  } else {
    # 2 rerand, limit number of iterations 
    # 2.1 independent rerand
    if (cont) {
      for (iter in 2:max_iter) {
        Wt <- sample(1:n_clust, n1_clust)
        Mt <- compute_M_set(get_id(Wt,clust_size,n1_clust), X, prop_cov_inv)
        if (Mt < Mbest) {
          Wbest <- Wt
          Mbest <- Mt
        }
        if (Mt <= a) break
      }
    }
    
    # output the best rather than the last
    Wt <- Wbest
    Mt <- Mbest
  }
  
  # output binary representation of Wt
  Wout <- rep(0, n)
  Wout[get_id(Wt,clust_size,n1_clust)] <- 1
  return(Wout)
}

ClustAJSRR <- function(X, n_clust, n1_clust = round(n_clust / 2), n_fisher = 1, 
                       p_a = 0.001, max_iter = Inf, parallel = FALSE, n_cores = 1) {
  n <- nrow(X)
  clust_size = n / n_clust
  n1 <- n1_clust * clust_size
  
  if (p_a == 1) {
    cov_inv <- NULL
    prop_cov_inv <- NULL
  } else {
    cov_inv <- solve(cov(X))
    prop_cov_inv <- n1 * (1 - n1 / n) * cov_inv
  }
  
  a <- qchisq(p = p_a, df = ncol(X))
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      ClustAJSRR_one(X, n_clust, n1_clust, prop_cov_inv, a, max_iter)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      ClustAJSRR_one(X, n_clust, n1_clust, prop_cov_inv, a, max_iter)
    })
  }
  if (n_fisher == 1) {
    unlist(out)
  } else {
    out
  }
}


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
