library(tidyverse)

setwd("~/Desktop/rerandomization/VNSRR/JHPCE/output/save/sim_simple_small_0.2/")


for (n in c(30, 100, 500)) {
  for (p in c(2)) {
    for (r in c(0.5)) {
      for (S in c(1000, 10000)) {
        rp = p / n
        n1 = round(n * r)
        
        df0 = NULL
        for (bbb in 1:200) {
          
          filename = sprintf("sim_simple_small_n%s_p%s_nt%s_S%s_newsummary_%s_%s.RData",
                             n, p, n1, S, ifelse(n==500&S==10000,"all","all"), bbb)
          if (filename %in% list.files()) {
            load(filename)
            df0 = bind_rows(df0, infer_sum_tab)
            rm(infer_sum_tab)
          } else {
            # print(sprintf("%s missing!", filename))
          }
          
        }
        
        if (!is.null(df0)) {
          print(nrow(df0)/ifelse(n==500&S==10000,4,4)*5)
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
            mutate(across(c(M), ~round(.x,4))) %>%
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
