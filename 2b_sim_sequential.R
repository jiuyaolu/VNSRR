setwd("/Users/jiuyaolu/Desktop/rerandomization/VNSRR/")

library(tidyverse)
library(tictoc)
library(parallel)
source("rerand CR AJSRR PSRR.R")
source("rerand VNSRR sequential.R")
source("fisher.R")


test = FALSE
run_inference = FALSE

if (test) {
  n_cores = 5
  n_rep = 10
  n_fisher = 100
} else {
  n_cores = 5
  n_rep = 100
  n_fisher = 1000
}



# par
alpha = 0.1


for (ngroup in c(2)) {
  for (p in c(50, 250)) {
    for (S in c(1000)) {
      
      set.seed(2023)
      
      n = ngroup * p * 2
      n1 = round(n / 2)
      n_seq = rep(p*2, ngroup)
      load(sprintf("data/sim_sequential/sim_sequential_n%s_p%s.RData",n,p))
      prop_cov_inv <- lapply(1:n_rep, function(i) {
        n1 * (1 - n1 / n) * solve(cov(X[[i]]))
      })
      
      print(sprintf("n=%s ngroup=%s p=%s S=%s",n,ngroup,p,S))

# computation ----------------------------------------------------------------

# methods ----

if (test) {
  # method_list = c("SeqVNSRR")
  # method_list = c("SeqPSRR", "SeqVNSRR")
  method_list = c("SeqAJSRR","SeqPSRR","SeqVNSRR")
} else {
  method_list = c("SeqAJSRR","SeqPSRR","SeqVNSRR")
}
    
s_k_tab <- expand_grid(
  K = length(n_seq),
  S = S,
  budget = c("increase")
) %>%
  mutate(
    s_k = pmap(list(K, S, budget), function(K, S, budget) {
      allocate_S(S, K, p, n_seq, budget)
    })
  )
    
met_tab <- bind_rows(
  tibble(
    K      = length(n_seq),
    S      = 1,
    budget = "-",
    method = "SeqCR",
    s_k    = list(rep(1, K))
  ),
  expand_grid(
    method = method_list,
    s_k_tab
  )
) %>% 
  mutate(
    across(where(is.character), as_factor),
    n = n,
  ) %>% 
  arrange(S, budget, method) %>% 
  mutate(met_id = row_number(), .before = everything())

met_fun <- function(i_met, i_rep) {
  if (met_tab$method[i_met] == "SeqVNSRR") {
    time <- system.time(
      design <- SeqVNSRR(X[[i_rep]], n_seq, n_fisher = n_fisher, 
                         p_a_seq = 1 / met_tab$s_k[[i_met]],
                         max_iter_seq = 10 * met_tab$s_k[[i_met]])
    )[3] %>% unname()
  } else {
    pair_switch <- met_tab$method[i_met] == "SeqPSRR"
    time <- system.time(
      design <- seq_rerand(X[[i_rep]], n_seq, n_fisher = n_fisher, 
                           p_a_seq = 1 / met_tab$s_k[[i_met]],
                           pair_switch = pair_switch, 
                           max_iter_seq = 10 * met_tab$s_k[[i_met]])
    )[3] %>% unname()
  }

  lst(design, time)
}


# design & time ----
# list structure: met, rep, fisher
design_time_met_rep_fisher <- map(1:nrow(met_tab), function(i_met) {
  tic()
  out <- mclapply(1:n_rep, function(i_rep) {
    met_fun(i_met, i_rep)
  }, mc.cores = n_cores)
  toc()
  out
})

design_met_rep_fisher <- map(design_time_met_rep_fisher, ~ {
  map(.x, ~ {
    .x$design
  })
})

destime_met_rep_fisher <- map(design_time_met_rep_fisher, ~ {
  map_dbl(.x, ~ {
    .x$time
  })
})

if (!test) {
  save(design_met_rep_fisher, destime_met_rep_fisher,
       file = sprintf("output/sim_sequential/sim_sequential_n%s_p%s_nt%s_S%s_design_%s.Rdata",
                      n, p, n1, S, "all"))
}


# inference & time ----
# list structure: R2, met, rep
if (run_inference) {
infer_time_R2_met_rep <- map(design_met_rep_fisher, function(design_rep_fisher) {
  tic()
  out <- mclapply(seq_along(design_rep_fisher), function(i_design_fisher) {
    design_fisher = design_rep_fisher[[i_design_fisher]]
    W <- design_fisher[[1]]
    
    M <- compute_M(W, X[[i_design_fisher]], prop_cov_inv[[i_design_fisher]])
    rand_metric <- compute_random(design_fisher)
    # M <- NA
    # rand_metric <- tibble(En = NA, Dn = NA, Ln = NA)
    
    # H0
    Y <- Y0[[i_design_fisher]]
    tauhat <- mean(Y[W == 1]) - mean(Y[W == 0])
    p0 <- fisher_p(W, Y, W_ref = design_fisher)
    time_solve <- system.time(
      ci <- fisher_ci(W, Y, W_ref = design_fisher, alpha = alpha)
    )[3] %>% unname()
    # time_bi <- system.time(
    #   ci_bi <- fisher_ci_bi(W, Y, W_ref = design_fisher, alpha = alpha)
    # )[3] %>% unname()
    
    # H1
    Y <- Y1[[i_design_fisher]] * W + Y0[[i_design_fisher]] * (1 - W)
    p1 <- fisher_p(W, Y, W_ref = design_fisher)
    
    # output
    # ci_bi <- ci_bi %>% dplyr::rename(ci_low_bi = ci_low, ci_up_bi = ci_up)
    lst(
      # infer = tibble(M, rand_metric, tauhat, p0, p1, ci, ci_bi),
      infer = tibble(M, rand_metric, tauhat, p0, p1, ci),
      # time  = tibble(time_solve, time_bi)
      time  = tibble(time_solve)
    )
  }, mc.cores = n_cores)
  toc()
  out
})

infer_R2_met_rep <- map(infer_time_R2_met_rep, ~ {
  map(.x, ~ {
    .x$infer
  })
})

inftime_R2_met_rep <- map(infer_time_R2_met_rep, ~ {
  map(.x, ~ {
    .x$time
  })
})

if (!test) {
  save(infer_R2_met_rep, inftime_R2_met_rep,
       file = sprintf("output/sim_sequential/sim_sequential_n%s_p%s_nt%s_S%s_inference_%s.RData",
                      n, p, n1, S, "all"))
}
}

# summarize ---------------------------------------------------------------

# statistical performance ----

if (run_inference) {
infer_detail_tab <- infer_R2_met_rep %>%
  tibble() %>%
  mutate(met_tab %>% dplyr::select(-s_k), .before = everything()) %>%
  relocate(K, S) %>% 
  unnest(cols = c(.)) %>%
  unnest(cols = c(.))

infer_sum_tab <- infer_detail_tab %>%
  group_by(K, met_id, S, method) %>%
  summarise(
    across(c(M, En, Dn, Ln), mean),
    bias    = abs(mean(tauhat) - 0),
    sd      = sd(tauhat),
    moment1 = mean(tauhat),
    moment2 = mean(tauhat^2),
    size    = mean(p0 <= alpha),
    power   = mean(p1 <= alpha),
    cp      = mean(ci_low <= 0 & ci_up >= 0),
    length  = mean(ci_up - ci_low),
    # cp_bi      = mean(ci_low_bi <= 0 & ci_up_bi >= 0),
    # length_bi  = mean(ci_up_bi - ci_low_bi),
    .groups = "drop"
  )
}

# computational performance ----

# remove the time of the first {n_cores} replications
time_sum_tab <- met_tab %>% 
  dplyr::select(-s_k) %>%
  mutate(
    des_time = map_dbl(destime_met_rep_fisher, ~ {
      # mean(.x[-(1:n_cores)])
      mean(.x)
    }) / n_fisher * 1000,
    # ci_time_sol = map_dbl(inftime_R2_met_rep, ~ {
    #   # mean(map_dbl(.x, ~ {.x$time_solve})[-(1:n_cores)])
    #   mean(map_dbl(.x, ~ {.x$time_solve}))
    # }),
    # ci_time_bi = map_dbl(inftime_R2_met_rep, ~ {
    #   # mean(map_dbl(.x, ~ {.x$time_bi})[-(1:n_cores)])
    #   mean(map_dbl(.x, ~ {.x$time_bi}))
    # })
  ) %>% 
  relocate(S)

if (!test) {
  if (run_inference) {
    save(infer_sum_tab, time_sum_tab,
         file = sprintf("output/sim_sequential/sim_sequential_n%s_p%s_nt%s_S%s_summary_%s.RData",
                        n, p, n1, S, "all"))
  } else {
    save(time_sum_tab,
         file = sprintf("output/sim_sequential/sim_sequential_n%s_p%s_nt%s_S%s_summary_%s.RData",
                        n, p, n1, S, "all"))
  }
}


if (run_inference) print(infer_sum_tab %>% as.data.frame)
print(time_sum_tab %>% as.data.frame)




    }
  }
}
