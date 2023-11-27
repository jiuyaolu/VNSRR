setwd("/Users/jiuyaolu/Desktop/rerandomization/VNSRR/")

library(tidyverse)
library(tictoc)
library(parallel)
source("rerand CR AJSRR PSRR.R")
source("rerand VNSRR simple.R")
source("fisher.R")
set.seed(2023)

test <- FALSE

if (test) {
  n_cores <- 5
  n_rep <- 5
  n_fisher <- 1000
} else {
  n_cores <- 5
  n_rep <- 1000
  n_fisher <- 1000
}

# data
load("data/reserpine/reserpine.RData")
X <- reserpine %>% 
  select(where(
    ~ attr(.x, "label") == "predictor"
  )) %>% 
  as.matrix()
X = X[,c("BIRTHDT", "WEIGHT", "HR")]
Yz <- reserpine %>% 
  select(where(
    ~ str_detect(attr(.x, "label"), "imputed outcome")
  ))
Y_real <- reserpine %>% 
  select(where(
    ~ attr(.x, "label") == "outcome"
  )) %>% 
  pluck(1)
S_value <- 1000
data_name <- "reserpine"


n <- nrow(X)
n1 <- 20
prop_cov_inv <- n1 * (1 - n1 / n) * solve(cov(X))


# par
alpha <- 0.1


# computation ----------------------------------------------------------------

# methods ----

met_tab <- bind_rows(
  expand_grid(
    S      = 1,
    method = "CR"
  ),
  expand_grid(
    S      = S_value,
    method = c("RR","PSRR","VNSRR")
  )
) %>% 
  mutate(
    across(where(is.character), as_factor),
    n = n
  ) %>% 
  mutate(met_id = row_number(), .before = everything())

# if (test) {
#   met_tab = met_tab %>% filter(method != "RR")
# }

met_fun <- function(i) {
  if (met_tab$method[i] == "VNSRR") {
    time = system.time(
      design <- VNSRR(X, n1 = n1, n_fisher = n_fisher, p_a = 1 / met_tab$S[i])
    )[3] %>% unname()
  } else {
    pair_switch <- met_tab$method[i] == "PSRR"
    time <- system.time(
      design <- rerand(X, n1 = n1, n_fisher = n_fisher, p_a = 1 / met_tab$S[i], 
                       pair_switch = pair_switch)
    )[3] %>% unname()
  }
  
  lst(design, time)
}

# design & time ----
# list structure: met, rep, fisher
design_time_met_rep_fisher <- map(1:nrow(met_tab), function(i_met) {
  tic()
  out <- mclapply(1:n_rep, function(i_rep) {
    met_fun(i_met)
  }, mc.cores = n_cores)
  toc()
  
  if (!test) {
    save(out, file = str_glue("output/real_reserpine/rdns_des_m{i_met}.RData"))
  }
  
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



# inference & time ----
# list structure: dataset, met, rep

infer_time_met_rep <- map(design_met_rep_fisher, function(design_rep_fisher) {
  tic()
  out <- mclapply(design_rep_fisher, function(design_fisher) {
    W <- design_fisher[[1]]
    M <- compute_M(W, X, prop_cov_inv)
    rand_metric <- compute_random(design_fisher)
    
    # H0
    Y <- Y_real # impute missing outcomes with sharp null
    tauhat <- mean(Y[W == 1]) - mean(Y[W == 0])
    p0 <- fisher_p(W, Y, W_ref = design_fisher)
    time_solve <- system.time(
      ci <- fisher_ci(W, Y, W_ref = design_fisher, alpha = alpha)
    )[3] %>% unname()
    # time_bi <- system.time(
    #   ci_bi <- fisher_ci_bi(W, Y, W_ref = design_fisher, alpha = alpha)
    # )[3] %>% unname()
    
    # H1
    Y <- Yz$Y1 * W + Yz$Y0 * (1 - W)
    p1 <- fisher_p(W, Y, W_ref = design_fisher)
    
    # output
    # ci_bi <- ci_bi %>% rename(ci_low_bi = ci_low, ci_up_bi = ci_up)
    lst(
      # infer = tibble(M, rand_metric, tauhat, p0, p1, ci, ci_bi),
      infer = tibble(M, rand_metric, tauhat, p0, p1, ci),
      # time  = tibble(time_solve, time_bi),
      time  = tibble(time_solve)
    )
  }, mc.cores = n_cores)
  toc()
  out
})

infer_time_R2_met_rep <- list(infer_time_met_rep)

infer_R2_met_rep <- map(infer_time_R2_met_rep, ~ {
  map(.x, ~ {
    map(.x, ~ {
      .x$infer
    })
  })
})

inftime_R2_met_rep <- map(infer_time_R2_met_rep, ~ {
  map(.x, ~ {
    map(.x, ~ {
      .x$time
    })
  })
})

if (!test) {
  save(infer_R2_met_rep, inftime_R2_met_rep,
       file = str_glue("output/real_reserpine/rdns_inf.RData"))
}

# summarize ---------------------------------------------------------------

met_tab <- met_tab %>% 
  mutate(across(c(S), ~ {ifelse(is.na(.x), "-", .x)}))

# statistical performance ----

infer_detail_tab <- infer_R2_met_rep %>%
  tibble() %>%
  mutate(dataset = data_name) %>%
  unnest(cols = c(.)) %>%
  group_by(dataset) %>%
  mutate(met_tab, .before = everything()) %>%
  ungroup() %>% 
  relocate(n, dataset) %>% 
  unnest(cols = c(.)) %>%
  unnest(cols = c(.))

infer_sum_tab <- infer_detail_tab %>%
  group_by(n, dataset, met_id, S, method) %>%
  summarise(
    across(c(M, En, Dn, Ln), mean),
    bias    = abs(mean(tauhat) - 0),
    sd      = sd(tauhat),
    size    = mean(p0 <= alpha),
    power   = mean(p1 <= alpha),
    cp      = mean(ci_low <= 0 & ci_up >= 0),
    length  = mean(ci_up - ci_low),
    # cp_bi      = mean(ci_low_bi <= 0 & ci_up_bi >= 0),
    # length_bi  = mean(ci_up_bi - ci_low_bi),
    .groups = "drop"
  )


# computational performance ----

# remove the time of the first {n_cores} replications
time_sum_tab <- met_tab %>% 
  mutate(
    des_time = map_dbl(destime_met_rep_fisher, ~ {
      mean(.x)
    }) / n_fisher * 10000,
    # ci_time_sol = map_dbl(inftime_R2_met_rep[[1]], ~ {
    #   mean(map_dbl(.x, ~ {.x$time_solve}))
    # }),
    # ci_time_bi = map_dbl(inftime_R2_met_rep[[1]], ~ {
    #   mean(map_dbl(.x, ~ {.x$time_bi}))
    # })
  )

if (!test) {
  save(
    infer_detail_tab, infer_sum_tab,
    time_sum_tab, 
    file = str_glue("output/real_reserpine/rdns_sum.RData")
  )
}

time_sum_tab



infer_sum_tab %>% 
  select(-(n:S)) %>% 
  mutate(across(c(M), ~round(.x,1))) %>%
  mutate(across(c(En), ~round(100*.x))) %>%
  mutate(across(c(Dn), ~round(100*.x))) %>%
  mutate(across(c(Ln), ~round(100*.x))) %>%
  mutate(across(c(bias, size, cp), ~round(100*.x,1))) %>%
  mutate(across(c(sd, power, length), ~round(100*.x)))
