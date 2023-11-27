library(tidyverse)

setwd("~/Desktop/rerandomization/VNSRR/output/sim_simple_small/")

df = NULL

for (n in c(30, 100, 500)) {
  for (p in c(2)) {
    for (r in c(0.5)) {
      for (S in c(1000, 10000)) {
        rp = p / n
        n1 = round(n * r)
        
        filename = sprintf("sim_simple_small_n%s_p%s_nt%s_S%s_summary_%s.RData",
                           n, p, n1, S, ifelse(n==500&S==10000,"all","all"))
        if (filename %in% list.files()) {
          load(filename)
          df = bind_rows(df, time_sum_tab)
        }
      }
    }
  }
}

df %>%
  mutate(time = round(des_time,1)) %>%
  dplyr::select(n, S, method, time) %>%
  as.data.frame()
