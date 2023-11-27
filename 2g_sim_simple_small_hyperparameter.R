setwd("/Users/jiuyaolu/Desktop/rerandomization/VNSRR/")

library(tidyverse)
library(tictoc)
library(parallel)
source("rerand CR AJSRR PSRR.R")
source("rerand VNSRR simple.R")
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


for (n in c(30, 100, 500)) {
  for (p in c(2)) {
    for (r in c(0.5)) {
      for (S in c(1000)) {
        
        set.seed(2023)
        
        rp = p / n
        n1 = round(n * r)
        load(sprintf("data/sim_simple_small/sim_simple_small_n%s_p%s.RData",n,p))
        prop_cov_inv <- lapply(1:n_rep, function(i) {
          n1 * (1 - n1 / n) * solve(cov(X[[i]]))
        })
        
        print(sprintf("n=%s p=%s nt=%s S=%s",n,p,n1,S))
        
        # computation ----------------------------------------------------------------
        
        # methods ----
        
        K_list = c(1, 
                   round(min(n1,n-n1) * 0.25), 
                   round(min(n1,n-n1) * 0.5), 
                   round(min(n1,n-n1) * 0.75), 
                   min(n1,n-n1))
        met_tab = bind_rows(
          expand_grid(
            S = S,
            K_local = K_list,
            K_shaking = K_list,
          )
        ) %>% 
          mutate(
            across(where(is.character), as_factor),
            n = n,
            rp = rp,
            p = p,
            r = r,
            n1 = n1,
          ) %>%
          mutate(met_id = row_number(), .before = everything())
        
        met_fun = function(i_met, i_rep) {
            time = system.time(
              design <- VNSRR(X[[i_rep]], n1 = n1, n_fisher = n_fisher, p_a = 1 / met_tab$S[i_met],
                              K_local = met_tab$K_local[i_met], K_shaking = met_tab$K_shaking[i_met],
                              max_iter = 10 * met_tab$S[i_met])
            )[3] %>% unname()

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
               file = sprintf("output/sim_simple_small_hyperparameter/sim_simple_small_hyperparameter_n%s_p%s_nt%s_S%s_design_%s.Rdata",
                              n, p, n1, S, ifelse(n==500&S==10000,"all","all")))
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
                 file = sprintf("output/sim_simple_small_hyperparameter/sim_simple_small_hyperparameter_n%s_p%s_nt%s_S%s_inference_%s.RData",
                                n, p, n1, S, ifelse(n==500&S==10000,"all","all")))
          }
        }
        
        # summarize ---------------------------------------------------------------
        
        # statistical performance ----
        
        if (run_inference) {
          infer_detail_tab <- infer_R2_met_rep %>%
            tibble() %>%
            mutate(met_tab, .before = everything()) %>%
            relocate(n, rp, p, r, n1, S) %>% 
            unnest(cols = c(.)) %>%
            unnest(cols = c(.))
          
          infer_sum_tab <- infer_detail_tab %>%
            group_by(n, rp, p, r, n1, met_id, S, method) %>%
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
          relocate(n, rp, p, r, n1, S)
        
        if (!test) {
          if (run_inference) {
            save(infer_sum_tab, time_sum_tab,
                 file = sprintf("output/sim_simple_small_hyperparameter/sim_simple_small_hyperparameter_n%s_p%s_nt%s_S%s_summary_%s.RData",
                                n, p, n1, S, ifelse(n==500&S==10000,"all","all")))
          } else {
            save(time_sum_tab,
                 file = sprintf("output/sim_simple_small_hyperparameter/sim_simple_small_hyperparameter_n%s_p%s_nt%s_S%s_summary_%s.RData",
                                n, p, n1, S, ifelse(n==500&S==10000,"all","all")))
          }
        }
        
        
        if (run_inference) print(infer_sum_tab %>% as.data.frame)
        print(time_sum_tab %>% select(-rp) %>% as.data.frame)
        
        
        
      }
    }
  }
}

