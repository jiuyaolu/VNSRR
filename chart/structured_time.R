library(tidyverse)


for (scheme in c("sequential","stratified","cluster")) {
  df = NULL
  
  for (n in c(100, 500)) {
    for (rp in c(0.5)) {
      for (r in c(0.5)) {
        for (S in c(1000)) {
        
          p = round(n * rp)
          n1 = round(n * r)
          
          setwd(sprintf("~/Desktop/rerandomization/VNSRR/output/sim_%s/",scheme))
          filename = sprintf("sim_%s_n%s_p%s_nt%s_S%s_summary_all.RData",
                             scheme, 2*n, p, 2*n1, S)
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


