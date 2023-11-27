compute_Mstar_set_clust <- function(Wt, i, j, Mt, H, h) {
  out = Mt
  
  Wtmp1 = Wtmp2 = Wt
  for (t in seq_along(i)) {
    Wtmp1 = Wtmp2
    Wtmp2[Wtmp2 == i[t]] = j[t]
    out = compute_Mstar_set(Wtmp1, Wtmp2, i[t], j[t], out, H, h)
  }
  
  return(out)
}

ClustVNSRR_one <- function(X, n_clust, n1_clust, prop_cov_inv, H, h,
                           a, max_iter, K_local, K_shaking) {
  n = nrow(X)
  clust_size = n / n_clust
  id_clust = lapply(1:n_clust, function(k) 1:clust_size + (k-1)*clust_size)
  
  # start point
  Wt = sample(1:n_clust, n1_clust)
  nonWt = setdiff(1:n_clust, Wt)
  Mt = compute_M_set(id_clust[Wt] %>% unlist, X, prop_cov_inv)
  
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
      Mstar <- compute_Mstar_set_clust(id_clust[Wt] %>% unlist, 
                                       id_clust[i[k]] %>% unlist, 
                                       id_clust[j[k]] %>% unlist, 
                                       Mt, H, h)
      
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
        Mt <- compute_Mstar_set_clust(id_clust[Wt0] %>% unlist, 
                                      id_clust[i[k]] %>% unlist, 
                                      id_clust[j[k]] %>% unlist, 
                                      Mt, H, h)
      }
      
      iter = iter + K_shaking
      if (iter >= max_iter) {
        cont <- FALSE
        break
      }
    }
  }
  
  # output binary representation of Wt
  Wout <- rep(0, n)
  Wout[id_clust[Wbest] %>% unlist] <- 1
  return(Wout)
}

ClustVNSRR <- function(X, n_clust, n1_clust = round(n_clust / 2), n_fisher = 1, 
                       p_a = 0.001, max_iter = Inf, parallel = FALSE, n_cores = 1, 
                       K_local = min(n1_clust,n_clust-n1_clust), K_shaking = 1) {
  n <- nrow(X)
  clust_size = n / n_clust
  n1 = n1_clust * clust_size
  
  cov_inv <- solve(cov(X))
  prop_cov_inv <- n1 * (1 - n1 / n) * cov_inv
  H <- X %*% cov_inv %*% t(X) / (n1 * (1 - n1 / n))
  h <- 2 * (n1 / n) * rowSums(H)
  
  a <- qchisq(p = p_a, df = ncol(X))
  
  if (parallel) {
    out <- mclapply(1:n_fisher, function(x) {
      ClustVNSRR_one(X, n_clust, n1_clust, prop_cov_inv, H, h,
                     a, max_iter, K_local, K_shaking)
    }, mc.cores = n_cores)
  } else {
    out <- lapply(1:n_fisher, function(x) {
      ClustVNSRR_one(X, n_clust, n1_clust, prop_cov_inv, H, h,
                     a, max_iter, K_local, K_shaking)
    })
  }
  
  if (n_fisher == 1) {
    unlist(out)
  } else {
    out
  }
}
