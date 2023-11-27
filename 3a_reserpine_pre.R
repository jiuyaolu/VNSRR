setwd("/Users/jiuyaolu/Desktop/rerandomization/VNSRR/")

library(tidyverse)
set.seed(2023)


# 1 wrangling -------------------------------------------------------------

# X
DEMOG <- read_csv("data/reserpine/ascii-crf-files_nida-cpu-0006/DEMOG.csv") %>% 
  dplyr::select(SUBJID, GENDER, BIRTHDT, HEIGHT, WEIGHT) %>% 
  mutate(GENDER = ifelse(GENDER == 2, 0, 1)) %>% 
  mutate(across(
      !SUBJID, 
      ~ structure(.x, "label" = "predictor")
  ))

DDRCPHYS_PRE <- read_csv("data/reserpine/ascii-crf-files_nida-cpu-0006/DDRCPHYS.csv") %>% 
  filter(VISID == "Day 1", RELTM == "PRE -15 min") %>% 
  dplyr::select(SUBJID, HR, SBP, DBP, RR) %>% 
  mutate(across(
    !SUBJID, 
    ~ structure(.x, "label" = "predictor")
  ))

# W
TRTCOND <- read_csv("data/reserpine/ascii-crf-files_nida-cpu-0006/TRTCOND.csv") %>% 
  dplyr::select(SUBJID, DRUGCD) %>% 
  mutate(DRUGCD = factor(case_when(
    DRUGCD == 1 ~ "reserpine 0.5mg",
    DRUGCD == 2 ~ "reserpine 1mg",
    DRUGCD == 3 ~ "placebo"
  ), levels = c("reserpine 1mg", "reserpine 0.5mg", "placebo"))) %>% 
  mutate(across(
    !SUBJID, 
    ~ structure(.x, "label" = "treatment")
  ))

# Y
DDRCPHYS <- read_csv("data/reserpine/ascii-crf-files_nida-cpu-0006/DDRCPHYS.csv") %>%
  filter(
    VISID == "Day 4", 
    RELTM != "PRE -15 min"
  ) %>% 
  group_by(SUBJID) %>% 
  summarise(HRPOST = round(mean(HR))) %>% 
  mutate(across(
    !SUBJID, 
    ~ structure(.x, "label" = "outcome")
  ))

dt_final <- list(DEMOG, DDRCPHYS_PRE, TRTCOND, DDRCPHYS) %>% 
  reduce(left_join, by = "SUBJID") %>% 
  dplyr::select(-SUBJID) %>% 
  mutate(DRUGCD = structure(
    ifelse(str_detect(DRUGCD, "reserpine"), "reserpine", "placebo"),
    "label" = "treatment"
  )) %>% 
  arrange(DRUGCD)

dt_final %>% group_by(DRUGCD) %>% summarise(across(everything(), mean))



# 2 matching & imputation -------------------------------------------------

X <- dt_final %>% 
  dplyr::select(where(
    ~ attr(.x, "label") == "predictor"
  )) %>% 
  as.matrix()
Y <- dt_final %>% 
  dplyr::select(where(
    ~ attr(.x, "label") == "outcome"
  )) %>% 
  mutate(across(everything(), ~ {attributes(.x) <- NULL; .x})) %>% 
  pull()
W <- dt_final %>% 
  dplyr::select(where(
    ~ attr(.x, "label") == "treatment"
  )) %>% 
  pull()


Y1 <- Y0 <- Y
Z <- ifelse(W == "reserpine", 1, 0)
impute_Y <- function(z) {
  idz <- which(Z == z)
  X_pred <- as_tibble(X)
  data_z <- tibble(Y, X_pred)[idz,]
  fit <- lm(Y ~ ., data = data_z)
  Yhat <- predict(fit, newdata = X_pred)
  Y_imp <- Yhat + rnorm(length(Yhat), mean = 0, sd = sqrt(var(Yhat[-idz]) / 10))
  Y_imp
}
Y1[Z == 0] <- round(impute_Y(1)[Z == 0])
Y0[Z == 1] <- round(impute_Y(0)[Z == 1])

reserpine <- dt_final %>% mutate(
  Y1 = structure(Y1, "label" = "imputed outcome under reserpine 1mg"),
  Y0 = structure(Y0, "label" = "imputed outcome under placebo")
  )

# 3 examination -------------------------------------------------


hist(reserpine$Y1)
hist(reserpine$Y0)

lm(Y1 ~ X) %>% summary
lm(Y0 ~ X) %>% summary

t.test(Y1, Y0)


# 4 save ------------------------------------------------------------------

save(reserpine, file = "data/reserpine/reserpine.RData")






