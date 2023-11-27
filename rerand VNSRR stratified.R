StratVNSRR_one <- function(X, n_seq, n1_seq, prop_cov_inv, H, h,
                      a, max_iter, K_local_seq, K_shaking_seq) {
  nstrata = length(n_seq)
  n_cum = c(0, cumsum(n_seq)[-nstrata])
  id_strata = lapply(1:nstrata, function(k) 1:n_seq[k]+n_cum[k])
  
  # start point
  Wt <- lapply(1:nstrata, function(k) sample(id_strata[[k]], n1_seq[k]))
  nonWt = lapply(1:nstrata, function(k) setdiff(id_strata[[k]], Wt[[k]]))
  Mt = compute_M_set(unlist(Wt), X, prop_cov_inv)
  
  cont = Mt > a
  
  Wbest = Wt
  Mbest = Mt
  
  iter = 0
  
  K_local_total = sum(K_local_seq)
  K_shaking_total = sum(K_shaking_seq)

  while (cont) {
    whichstrata = lapply(1:nstrata, function(k) rep(k,K_local_seq[k])) %>% unlist
    i <- lapply(1:nstrata, function(k) sample(Wt[[k]], K_local_seq[k])) %>% unlist
    j <- lapply(1:nstrata, function(k) sample(nonWt[[k]], K_local_seq[k])) %>% unlist
    
    perm = sample(1:K_local_total)
    whichstrata = whichstrata[perm]
    i = i[perm]
    j = j[perm]
    
    swap = FALSE
    
    for (k in 1:K_local_total) {
      Wstar <- Wt
      Wstar[[whichstrata[k]]][Wstar[[whichstrata[k]]] == i[k]] <- j[k]
      Mstar <- compute_Mstar_set(unlist(Wt), unlist(Wstar), i[k], j[k], Mt, H, h)
      
      if (Mstar < Mt) {
        Wt = Wstar
        Mt = Mstar
        nonWt[[whichstrata[k]]][nonWt[[whichstrata[k]]] == j[k]] = i[k]
        
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
    
    iter = iter + K_local_total
    if (iter >= max_iter) {
      cont <- FALSE
      break
    }
    
    # shaking
    if (!swap) {
      whichstrata = lapply(1:nstrata, function(k) rep(k,K_shaking_seq[k])) %>% unlist
      i <- lapply(1:nstrata, function(k) sample(Wt[[k]], K_shaking_seq[k])) %>% unlist
      j <- lapply(1:nstrata, function(k) sample(nonWt[[k]], K_shaking_seq[k])) %>% unlist
      
      for (k in 1:K_shaking_seq) {
        Wt0 = Wt
        Wt[[whichstrata[k]]][Wt[[whichstrata[k]]] == i[k]] <- j[k]
        nonWt[[whichstrata[k]]][nonWt[[whichstrata[k]]] == j[k]] <- i[k]
        Mt <- compute_Mstar_set(unlist(Wt0), unlist(Wt), i[k], j[k], Mt, H, h)
      }
      
      iter = iter + K_shaking_total
      if (iter >= max_iter) {
        cont <- FALSE
        break
      }
    }
  }
  
  # output binary representation of Wt
  Wout <- rep(0, n)
  Wout[unlist(Wbest)] <- 1
  return(Wout)
  
}

StratVNSRR <- function(X, n_seq, n1_seq = round(n_seq / 2), n_fisher = 1, 
                  p_a = 0.001, max_iter = Inf, parallel = FALSE, 
                  n_cores = 1, 
                  K_local_seq = pmin(n1_seq, n_seq - n1_seq), 
                  K_shaking_seq = rep(1, length(n_seq))) {
  n <- nrow(X)
  n1 = sum(n1_seq)
  cov_inv <- solve(cov(X))
  prop_cov_inv <- n1 * (1 - n1 / n) * cov_inv
  H <- X %*% cov_inv %*% t(X) / (n1 * (1 - n1 / n))
  h <- 2 * (n1 / n) * rowSums(H)
  a <- qchisq(p = p_a, df = ncol(X))
  
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      StratVNSRR_one(X, n_seq, n1_seq, prop_cov_inv, H, h,
                a, max_iter, K_local_seq, K_shaking_seq)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      StratVNSRR_one(X, n_seq, n1_seq, prop_cov_inv, H, h,
                a, max_iter, K_local_seq, K_shaking_seq)
    })
  }
  
  if (n_fisher == 1) {
    unlist(out)
  } else {
    out
  }
}
