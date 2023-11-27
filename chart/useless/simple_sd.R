n = 500
S = 1000
rp = 0.5
r = 0.5
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

df0 = NULL
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
    
    df0 = bind_rows(df0, infer_detail_tab)
    
    rm(infer_detail_tab)
  }
  
}

for (m in unique(df0$method)) {
  tau = df0$tauhat[df0$method==m]
  cat(
    m,
    sqrt(sum(tau^2) / length(tau)),
    sqrt(sum(tau^2) / qchisq(0.025,length(tau))),
    sqrt(sum(tau^2) / qchisq(0.975,length(tau))),
    "\n"
  )
}
