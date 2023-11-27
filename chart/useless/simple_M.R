library(tidyverse)

setwd("~/Desktop/rerandomization/VNSRR/JHPCE/output/sim_simple/")


df0 = NULL

for (n in c(30, 100, 500)) {
  for (rp in c(0.5)) {
    for (r in c(0.5)) {
      for (S in c(1000, 10000)) {
        p = round(n * rp)
        n1 = round(n * r)
        
        met_tab = bind_rows(
          expand_grid(
            S = 1,
            method = "CR",
          ),
          expand_grid(
            S = S,
            method = c("AJSRR","PSRR","VNSRR"),
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
        
        if (n==500 & S==10000) {
          met_tab = met_tab %>% filter(method != "AJSRR")
        }
        
        
        for (bbb in 1:200) {
          filename = sprintf("sim_simple_n%s_p%s_nt%s_S%s_inference_%s_%s.RData",
                             n, p, n1, S, ifelse(n==500&S==10000,"noAJS","all"), bbb)
          if (filename %in% list.files()) {
            load(filename)
            infer_detail_tab <- infer_R2_met_rep %>%
              tibble() %>%
              mutate(met_tab, .before = everything()) %>%
              relocate(n, rp, p, r, n1, S) %>% 
              unnest(cols = c(.)) %>%
              unnest(cols = c(.))
            infer_detail_tab <- infer_detail_tab %>% dplyr::select(n, S, method, M)
            df0 = bind_rows(df0, infer_detail_tab)
          }
          
        }
        
      }
    }
  }
}


df = df0 %>%
  filter(method != "CR") %>%
  mutate(S_label = paste0("S = ",S),
         n_label = paste0("n = ",n))
df$n_label = factor(df$n_label, levels = paste0("n = ",c(30,100,500)))
df %>%
  ggplot(aes(x = method,
             y = M)) +
  facet_grid(n_label ~ S_label, scales = "free_y") +
  geom_violin(trim = FALSE) +
  stat_summary(fun = mean) +
  geom_hline(aes(yintercept=qchisq(1/S, n/2)), 
             col="red",
             lty = 5)
