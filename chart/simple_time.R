library(tidyverse)

setwd("~/Desktop/rerandomization/VNSRR/output/sim_simple/")

df = NULL

for (n in c(100, 500)) {
  for (rp in c(0.5)) {
    for (r in c(0.5)) {
      for (S in c(1000)) {
        p = round(n * rp)
        n1 = round(n * r)
        
        filename = sprintf("sim_simple_n%s_p%s_nt%s_S%s_summary_%s.RData",
                           n, p, n1, S, ifelse(n==500&S==10000,"noAJS","all"))
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
