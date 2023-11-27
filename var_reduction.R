library(tidyverse)

p_a <- c(0.001, 0.0001, 0.00001)
p <- 10
a <- qchisq(p = p_a, df = p)
R2 <- 0.5

(va_psrr <- a / p)
(va_rr <- pchisq(a, p + 2) / pchisq(a, p))

(var_red_psrr <- (1 - va_psrr) * R2)
(var_red_rr <- (1 - va_rr) * R2)

(sum_tab <- tibble(
  p_a = 10^c(-3:-6),
  a   = qchisq(p = p_a, df = p),
  va_psrr = a / p,
  va_rr   = pchisq(a, p + 2) / pchisq(a, p),
  var_red_bound = (1 - va_psrr) * R2,
  var_red_rr   = (1 - va_rr) * R2
))

sum_tab %>% 
  select(p_a, var_red_bound, var_red_rr) %>% 
  pivot_longer(cols = 2:3, names_to = "type", values_to = "var_red") %>% 
  ggplot(aes(x = p_a, y = var_red, linetype = type)) +
  geom_line() +
  scale_x_log10() +
  theme_bw()






