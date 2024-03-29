
df <- read.csv("./data/Imputes/0706imputed10.csv")
unique(df$Treatment_Assignment)
df$p1outcome
# Note: Compare these numbers to Fig 1 in RQL
table(df$P1Def5day, df$Treatment_Assignment)

library(dplyr)
trt_df <- df %>% 
  mutate(A1 = P1Def5day,
         A2 = 1*(Treatment_Assignment %in% c("TDM","Naltrexone")),
         response_ind = 1*(Treatment_Assignment %in% c("TDM","Usual Care"))
         ) %>% 
  select(col1, A2, response_ind, A1)
trt_df %>% select(!col1) %>%  table
# compare the "active" groups in figure 1 to the lower and upper branches

final_df <- df %>% 
  select(
    col1,
    duration1,
    gender,
    educyears,
    race,
    alcyears,
    alcintoxlife,
    married,
    ethnic,
    ocdsP0,
    pacsP0,
    phase1AvgApc,
    pdhdS1,
    phase1AvgPacs,
    mcsP1
  ) %>% 
  mutate(y = duration1 / (7*24)) %>% 
  inner_join(trt_df, by="col1") %>% 
  select(!c(col1, duration1))

baseline_vars <- c(
  "gender","educyears", "race", "alcyears", "alcintoxlife", 
  "married", "ethnic", "ocdsP0", "pacsP0"
)

stg1_vars <- c(
  "phase1AvgApc", "pdhdS1", "phase1AvgPacs", "mcsP1"
)

h1 <- c(
  "gender", "race", "alcyears", "alcintoxlife",
  "ocdsP0"
)

h2 <- c(baseline_vars, stg1_vars, "A1", "response_ind")

z <- c(
  "gender", "A1", "alcintoxlife", "ocdsP0",
  "phase1AvgPacs", "mcsP1"
)
z_form <- paste0(z, collapse=" + ")

h2_form <- paste0(
  "y2c ~ a2c:(response_ind + I(1-response_ind) + ",
  "response_ind:gender + response_ind:A1 +",
  "I(1-response_ind):(", z_form, "))"
)

mu_a2 <- glm(A2 ~ response_ind*gender, data=final_df)
mu_a2_hat <- predict(mu_a2, type="response")
summary(mu_a2_hat)

mu_a1 <- glm(A1 ~ gender, data=final_df)
mu_a1_hat <- predict(mu_a1, type="response")
summary(mu_a1_hat)

library(ranger)

mu_y2 <- ranger(paste0("y ~ ", paste0(h2, collapse = " + ")),
                data=final_df,
                num.trees = 2000,
                mtry=5)

mu_y2_hat <- mu_y2$predictions

y2c <- final_df$y - mu_y2_hat
a2c <- final_df$A2 - mu_a2_hat

rob_Q2_fit <- lm(h2_form, data=final_df)

N_boot <- 5000
omegas <- replicate(N_boot, rexp(nrow(final_df)))
Q2_mat <- model.matrix(as.formula(h2_form), data=final_df)

boot_lm_2 <- function(idx) {
  coef(lm(y2c ~ Q2_mat - 1, weights=omegas[,idx]))
}

boot_thetas <- sapply(1:N_boot, boot_lm_2)
cis_2 <- apply(boot_thetas, 1, quantile, probs=c(.025, .975))
which(apply(sign(cis_2), 2, function(x) x[1] * x[2] >= 0))

rob_Q2_fit %>% coef %>% as.matrix
coef_Q2 <- rob_Q2_fit %>% coef
ncol(Q2_mat)

# build pretty cis
output_stg2 <- sapply(1:ncol(Q2_mat), function(col_idx){
  sprintf("%.2f (%.2f, %.2f)", coef_Q2[col_idx], cis_2[1,col_idx], cis_2[2,col_idx])
})


boot_uposi_h0_2_fun <- function(idx){
  e <- omegas[,idx]
  max(abs(crossprod(y2c * (e-1), Q2_mat)))
}
boot_uposi_2_h0 <- sapply(1:N_boot, boot_uposi_h0_2_fun)
uposi_2_h0_stat <- max(abs(crossprod(y2c, Q2_mat)))
mean(uposi_2_h0_stat >= boot_uposi_2_h0)