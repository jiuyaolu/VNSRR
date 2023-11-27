SeqVNSRR_part <- function(X_fix, X_wait, n1_wait, p_a, Wt_fix, 
                          M_fix, prop_cov_inv, H, h,
                          max_iter, K_local, K_shaking) {
  p <- ncol(X_wait)
  n_wait <- nrow(X_wait)
  
  if (is.null(X_fix)) {
    X <- X_wait
    n_fix <- 0
    M_fix <- 0
  } else {
    X <- rbind(X_fix, X_wait)
    n_fix <- nrow(X_fix)
  }
  
  n <- nrow(X)
  
  ncp_chi <- n_fix / n_wait * M_fix
  a <- qchisq(p = p_a, df = p, ncp = ncp_chi) * n_wait / n

  # start point
  Wt_wait = (n_fix + 1):n
  Wt = sample(Wt_wait, n1_wait)
  nonWt = setdiff(Wt_wait, Wt)
  Mt = compute_M_set(c(Wt_fix, Wt), X, prop_cov_inv)

  cont = Mt > a
  
  Wbest = Wt
  Mbest = Mt
  
  iter = 0
  
  while (cont) {
    i <- sample(Wt, K_local)
    j <- sample(nonWt, K_local)
    
    swap = FALSE
    
    for (k in 1:K_local) {
      Wstar <- Wt
      Wstar[Wstar == i[k]] <- j[k]
      Mstar <- compute_Mstar_set(c(Wt_fix, Wt), c(Wt_fix, Wstar), i[k], j[k], Mt, H, h)
      
      if (Mstar < Mt) {
        Wt = Wstar
        Mt = Mstar
        nonWt[nonWt == j[k]] = i[k]
        
        swap = TRUE
        
        if (Mt <= a) {
          cont <- FALSE
          break
        }
      } 
    }
    
    if (Mt < Mbest) {
      Wbest = Wt
      Mbest = Mt
    }
    
    iter = iter + K_local
    if (iter >= max_iter) {
      cont <- FALSE
      break
    }
    
    # shaking
    if (!swap) {
      i <- sample(Wt, K_shaking)
      j <- sample(nonWt, K_shaking)
      
      for (k in 1:K_shaking) {
        Wt0 = Wt
        Wt[Wt == i[k]] <- j[k]
        nonWt[nonWt == j[k]] <- i[k]
        Mt <- compute_Mstar_set(c(Wt_fix, Wt0), c(Wt_fix, Wt), i[k], j[k], Mt, H, h)
      }
      
      iter = iter + K_shaking
      if (iter >= max_iter) {
        cont <- FALSE
        break
      }
    }
  }
  
  # output set representation of Wt
  return(list(W = c(Wt_fix, Wt), M = Mt))
}


SeqVNSRR_one <- function(X, n_cum, n1_seq, p_a_seq,
                         prop_cov_inv_seq, H_seq, h_seq,
                         max_iter_seq, K_local, K_shaking) {
  WM_k <- list()
  WM_k$W <- NULL
  WM_k$M <- NULL
  
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
    } else {
      WM_k <- SeqVNSRR_part(X_fix = X_fix_k, X_wait = X_wait_k, n1_wait = n1_seq[k],
                            p_a = p_a_seq[k], Wt_fix = WM_k$W, M_fix = WM_k$M,
                            prop_cov_inv = prop_cov_inv_seq[[k]], 
                            H = H_seq[[k]], h = h_seq[[k]],
                            max_iter = max_iter_seq[k],
                            K_local = K_local[k], K_shaking = K_shaking[k])
      
    }
  }
  
  Wout <- rep(0, nrow(X))
  Wout[WM_k$W] <- 1
  return(Wout)
}


SeqVNSRR <- function(X, n_seq, p_a_seq, n1_seq = round(n_seq / 2), n_fisher = 1, 
                     max_iter_seq = rep(Inf, length(n_seq)),
                     parallel = FALSE, n_cores = 1,
                     K_local = pmin(n1_seq, n_seq-n1_seq), K_shaking = rep(1, length(n_seq))
) {
  n_cum <- cumsum(n_seq)
  n1_cum <- cumsum(n1_seq)
  
  
  pre_cum <- lapply(1:length(n_seq), function(k) {
    X_cum_k <- X[1:n_cum[k],]
    cov_inv <- solve(cov(X_cum_k))
    prop_cov_inv <- n1_cum[k] * (1 - n1_cum[k] / n_cum[k]) * cov_inv
    
    H = X_cum_k %*% cov_inv %*% t(X_cum_k) / (n1_cum[k] * (1 - n1_cum[k] / n_cum[k]))
    h = 2 * (n1_cum[k] / n_cum[k]) * rowSums(H)

    list(prop_cov_inv = prop_cov_inv, H = H, h = h)
  })
  prop_cov_inv_seq <- lapply(pre_cum, function(x) x$prop_cov_inv)
  H_seq <- lapply(pre_cum, function(x) x$H)
  h_seq <- lapply(pre_cum, function(x) x$h)
  
  
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      SeqVNSRR_one(X, n_cum, n1_seq, p_a_seq,
                   prop_cov_inv_seq, H_seq, h_seq,
                   max_iter_seq, K_local, K_shaking)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      SeqVNSRR_one(X, n_cum, n1_seq, p_a_seq,
                   prop_cov_inv_seq, H_seq, h_seq,
                   max_iter_seq, K_local, K_shaking)
    })
  }
  if (n_fisher == 1) {
    out[[1]]
  } else {
    out
  }
}

