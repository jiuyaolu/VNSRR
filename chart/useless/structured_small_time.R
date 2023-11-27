library(tidyverse)


for (scheme in c("sequential","stratified","cluster")) {
  df = NULL
  
  for (n in c(60, 1000)) {
    for (p in c(2)) {
      for (r in c(0.5)) {
        for (S in c(1000)) {
        
          rp = p / n
          n1 = round(n * r)
          
          setwd(sprintf("~/Desktop/rerandomization/VNSRR/output/sim_%s/",scheme))
          filename = sprintf("sim_%s_small_n%s_p%s_nt%s_S%s_summary_all.RData",
                             scheme, n, p, n1, S)
          if (filename %in% list.files()) {
            print(filename)
            load(filename)
            df = bind_rows(df, time_sum_tab %>% mutate(scheme=scheme) %>% relocate(scheme))
          }
        }
      }
    }
  }
  print(df %>% 
          mutate(time = round(des_time,1), .keep = "unused") %>% 
          as.data.frame)
}


