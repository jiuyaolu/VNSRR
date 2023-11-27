library(tidyverse)

setwd("~/Desktop/rerandomization/VNSRR/JHPCE/output/sim_simple/")


for (n in c(30, 100, 500)) {
  for (rp in c(0.5)) {
    for (r in c(0.5)) {
      for (S in c(1000, 10000)) {
        p = round(n * rp)
        n1 = round(n * r)
        
        df0 = NULL
        for (bbb in 1:200) {
          
          filename = sprintf("sim_simple_n%s_p%s_nt%s_S%s_newsummary_%s_%s.RData",
                             n, p, n1, S, ifelse(n==500&S==10000,"noAJS","all"), bbb)
          if (filename %in% list.files()) {
            load(filename)
            df0 = bind_rows(df0, infer_sum_tab)
            rm(infer_sum_tab)
          }
          
        }
        
        if (!is.null(df0)) {
          print(nrow(df0)/ifelse(n==500&S==10000,3,4)*5)
          df = df0 %>%
            group_by(n, S, method) %>%
            summarise(
              across(c(M, En, Dn, Ln, moment1, moment2, size, power, cp, length), mean),
              .groups = "drop"
            ) %>%
            mutate(bias = abs(moment1),
                   sd = sqrt(moment2 - moment1^2),
                   .keep = "unused",
                   .before = size) %>%
            mutate(across(c(M), ~round(.x,1))) %>%
            # mutate(across(c(En), ~round(.x,3))) %>%
            mutate(across(c(En), ~round(100*.x))) %>%
            # mutate(across(c(Dn), ~round(.x,2))) %>%
            mutate(across(c(Dn), ~round(100*.x))) %>%
            # mutate(across(c(Ln), ~round(.x,1))) %>%
            mutate(across(c(Ln), ~round(100*.x))) %>%
            mutate(across(c(bias, size, cp), ~round(100*.x,1))) %>%
            mutate(across(c(sd, power, length), ~round(100*.x))) %>%
            relocate(M, .after = length) %>%
            as.data.frame()
          print(df)
        }
      }
    }
  }
}
