library(tidyverse)



for (scheme in c("sequential","stratified","cluster")) {
  setwd(sprintf("~/Desktop/rerandomization/VNSRR/JHPCE/output/sim_%s/",scheme))
  
  for (n in c(60, 1000)) {
    for (S in c(1000)) {
      p = 2
      n1 = n / 2

      df0 = NULL
      for (bbb in 1:200) {
        
        filename = sprintf("sim_%s_small_n%s_p%s_nt%s_S%s_summary_%s_%s.RData",
                           scheme,n, p, n1, S, "all", bbb)
        if (filename %in% list.files()) {
          load(filename)
          df0 = bind_rows(df0, infer_sum_tab)
          rm(infer_sum_tab)
        } else {
          print(paste(filename,"missing!"))
        }
        
      }
      
      if (!is.null(df0)) {
        print(nrow(df0)/ifelse(scheme=="sequential",4,3)*5)
        df = df0 %>%
          group_by(method) %>%
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
          mutate(n = n) %>%
          relocate(n) %>%
          as.data.frame()
        print(df)
      }
    }
  }
}
