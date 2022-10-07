library(future)
library(future.apply)
plan(multicore, workers=4L)
# plan(sequential)
RNGkind("L'Ecuyer-CMRG")
library(mlr3verse)
library(lgr)
library(dplyr)

#---------------------------------------------------
# Try a better model selection procedure
#---------------------------------------------------
load("./results/run.RData")
print("------------------------")
print("Beginning alternative model selection...")
print("Starting from Stage 2...")
print("------------------------")


set.seed(2022)
# noise <- replicate(100-ncol(Q2_mat), rnorm(nrow(final_df)))
# # scale_df$noise <- scale(noise, T, F)
# noise <- scale(noise, T, F)
# colnames(noise) <- paste0("N.", 1:ncol(noise))
# Q2_mat <- cbind(Q2_mat, noise * a2c)

### Begin Selection
use_glmnet <- F
use_mcp <- T
use_lar <- F
M2_SIZE <- 3
pen_fac <- rep(1, ncol(Q2_big))
pen_fac[which(colnames(Q2_big) %in% c("a2c", "R:a2c", "NR:a2c"))] <- 0
if(use_mcp){
  library(ncvreg)
  Q2_sel_fit <- cv.ncvreg(Q2_big, y2c, penalty.factor = pen_fac, dfmax=20)
  col_idx <- Q2_sel_fit$min
  coefs_pen <- coef(Q2_sel_fit$fit, which=col_idx)[-1]
  plot(Q2_sel_fit)
} else if(use_glmnet) {
  library(glmnet)
  Q2_sel_fit <- glmnet(Q2_big, y2c, penalty.factor = pen_fac)
  col_idx <- which.min(abs(Q2_sel_fit$df - M2_SIZE))
  coefs_pen <- coef(Q2_sel_fit, s=Q2_sel_fit$lambda[col_idx])
  plot(Q2_sel_fit)
} else if(use_lar){
  # larfit2 <- selectiveInference::lar(Q2_mat, y2c)
  # coefs_pen <- c(1, coef(larfit2, s=M2_SIZE, mode="step"))
  # no intercept
  # larfit2 <- selectiveInference::lar(Q2_mat, y2c, intercept=F)
  # coefs_pen <- c(coef(larfit2, s=M2_SIZE+1, mode="step"))
  
  # separate responder/non-responder coefficients
  # larfit2_r <- lar(scale(Q2_mat[r2c==1,], T, F), scale(y2c[r2c==1], T, F), intercept=T)
  larfit2_r <- selectiveInference::lar(Q2_mat[r2c==1,], 
                                       y2c[r2c==1], intercept=T)
  larfit2_nr <- selectiveInference::lar(Q2_mat[r2c!=1,], 
                                        y2c[r2c!=1], intercept=T)
  coefs_pen <- c(coef(larfit2_r, s=M2_SIZE+1, mode="step"),
                 coef(larfit2_nr, s=M2_SIZE+1, mode="step"))
  # colnames(Q2_mat)[which(abs(coefs_pen) > 0)]
}


M2_hat <- which(abs(coefs_pen) > 0)
# M2_hat_Q2 <- M2_hat[-1] - 1
# no intercept
M2_hat_Q2 <- M2_hat



# plot(Q2_sel_fit)
# which(Q2_sel_fit$df == 6)
# coefs <- coef(Q2_sel_fit, s=Q2_sel_fit$lambda[which(Q2_sel_fit$df == 6)])
# M2_hat <- which(abs(coefs) > 1e-5)
# rownames(coefs)[M2_hat]
# 
# M2_hat_Q2 <- M2_hat[-1] - 1

# Q2_sel_mat <- cbind(1, scale(Q2_mat[,M2_hat_Q2], F, T))
# colnames(Q2_sel_mat)
n <- nrow(Q2_mat)
# whitening <- solve(chol(crossprod(Q2_mat[,M2_hat_Q2])/n))
# Q2_mat[,M2_hat_Q2] <- Q2_mat[,M2_hat_Q2] %*% whitening
Q2_sel_mat <- Q2_big[,M2_hat_Q2]
# Q2_sel_mat <- cbind(1, Q2_mat[,M2_hat_Q2])

rob_Q2_fit <- lm(y2c ~ Q2_sel_mat - 1)
summary(rob_Q2_fit)

ci_sand_stg2 <- construct_cis_sandwich(rob_Q2_fit) %>% 
  add_length(colnames=c("Lower", "Upper"),
             rownames=colnames(Q2_sel_mat))
print(ci_sand_stg2)

set.seed(2022)
construct_cis_Q <- function(boot_mat, coefs) {
  stopifnot(nrow(boot_mat) == length(coefs))
  deltas <- apply(boot_mat, 1, quantile, probs=c(0.025, 0.975))
  deltas <- unname(deltas)
  lower <- deltas[1,] + coefs
  upper <- deltas[2,] + coefs
  cbind(lower=lower, upper=upper)
}
boot_iter_2 <- function() {
  # omega <- rexp(n, rate=1)
  omega <- rpois(n, lambda=1)
  omega <- omega / mean(omega)
  coef(lm(y2c ~ Q2_sel_mat - 1, weights = omega)) - coef(rob_Q2_fit)
}
boots_rob_Q2 <- future_replicate(10000, boot_iter_2())
print("Bootstrap CIs:")
print("------------------")
ci_stg2_boot <- construct_cis_Q(boots_rob_Q2, coef(rob_Q2_fit))
ci_stg2_boot <- add_length(ci_stg2_boot, colnames=c("Lower", "Upper"),
                           rownames=colnames(Q2_sel_mat))
print(ci_stg2_boot)

### Begin UPoSI

set.seed(2022)
source("./R/uposi.R")
alpha <- .05
uposi_stg2 <- Uposi$new(
  Q2_big, y2c, stage=2, M = M2_hat_Q2, 
  coef=coef(rob_Q2_fit), alpha = alpha, Nboot = 10000, seed = 2022
)
uposi_stg2$do()
Pval2 <- uposi_stg2$pval_H0 # mean(uposi_stg1$boots[1,] <= test_stat_H01)
test_stat_H02 <- uposi_stg2$test_stat_H0
C_alpha_H02 <- uposi_stg2$C_alpha_H0
message(sprintf(paste("Stage: %i",
                      "UPoSI omnibus test statistic: %.3f",
                      "Critical value at alpha=%.2f: %.3f",
                      "UPoSI omnibus P-value: %.3f", sep="\n"), 
                2,
                test_stat_H02, alpha, C_alpha_H02, Pval2))

print("UPOSI CIs:")
print("------------------")
ci_stg2 <- add_length(uposi_stg2$ci_2, colnames=c("Lower", "Upper"),
                      rownames=colnames(Q2_sel_mat))
print(ci_stg2)


if(use_lar){
  print("SI CIs:")
  print("------------------")
  sigmahat <- selectiveInference::estimateSigma(Q2_big, y2c, standardize=F, intercept=F)$sigmahat
  # larinf2 <- selectiveInference::larInf(larfit2, alpha=alpha, type = "all", k = M2_SIZE)
  larinf2_r <- selectiveInference::larInf(larfit2_r, alpha=alpha, type = "all", k = M2_SIZE, sigma=sigmahat)
  larinf2_nr <- selectiveInference::larInf(larfit2_nr, alpha=alpha, type = "all", k = M2_SIZE, sigma=sigmahat)
  # print(larinf2)
  ci_si_stg2 <- rbind(larinf2_r$ci[order(larinf2_r$vars),], larinf2_nr$ci[order(larinf2_nr$vars),])
  ci_si_stg2 <- add_length(ci_si_stg2, colnames=c("Lower", "Upper"),
                           rownames=colnames(Q2_sel_mat))
  print(ci_si_stg2)
}


##### Stage 1
set.seed(9272022)

Delta2_hat <- as.numeric(Q2_sel_mat %*% coef(rob_Q2_fit))
blip <- (1*(Delta2_hat > 0) - final_df$A2)*Delta2_hat
Y1 <- y2c + blip

sl_fit1 <- metalearn(outcome = Y1,
                     treat = final_df$A1, 
                     stage = 1)
y1c <- Y1 - sl_fit1$outcome
a1c <- final_df$A1 - sl_fit1$treat
a1c <- final_df$A1 - 0.5

# sl_fit1$outcome_coefs %>% colMeans


# h2_form <- as.formula(paste0("y2c ~ a2c + a2c:r2c + a2c:I(1-r2c) + (a2c):(", paste0(h2_tailor, collapse=" + "), ") - 1"))
# h1_form <- as.formula(paste0("y1c ~ a1c + a1c:(", paste0(h1_tailor, collapse=" + "), ")^3 - 1"))
# h1_form <- as.formula(paste0("y1c ~ a1c + a1c:(", paste0(h1_tailor, collapse=" + "), ") - 1"))
# scale_df <- final_df %>%
  mutate(across(dplyr::all_of(setdiff(h1_tailor, dummy_vars)), scale))
Q1_mat <- model.matrix(as.formula(h1_form), data=final_df)
# Q1_mat <- scale(Q1_mat, T, F)


### Begin Selection
use_glmnet <- F
use_mcp <- T
use_lar <- F
M1_SIZE <- 5
pen_fac <- rep(1, ncol(Q1_mat))
pen_fac[which(colnames(Q1_mat) %in% c("a1c"))] <- 0
if(use_mcp){
  library(ncvreg)
  Q1_sel_fit <- cv.ncvreg(Q1_mat, y1c, penalty.factor = pen_fac, dfmax=20)
  col_idx <- Q1_sel_fit$min
  coefs_pen <- coef(Q1_sel_fit, which=col_idx)[-1]
  plot(Q1_sel_fit)
} else if(use_glmnet) {
  library(glmnet)
  Q1_sel_fit <- glmnet(Q1_mat, y1c, penalty.factor = pen_fac)
  col_idx <- which.min(abs(Q1_sel_fit$df - M1_SIZE))
  coefs_pen <- coef(Q1_sel_fit, s=Q1_sel_fit$lambda[col_idx])
  plot(Q1_sel_fit)
} else if(use_lar){
  # larfit1 <- selectiveInference::lar(Q1_mat, y1c)
  # coefs_pen <- c(1, coef(larfit1, s=M1_SIZE, mode="step"))
  # no intercept
  larfit1 <- selectiveInference::lar(Q1_mat, y1c, intercept=F)
  coefs_pen <- c(coef(larfit1, s=M1_SIZE+1, mode="step"))
}

M1_hat <- which(abs(coefs_pen) > 0)
# M1_hat_Q1 <- M1_hat[-1] - 1
# no intercept
M1_hat_Q1 <- M1_hat

Q1_sel_mat <- cbind(Q1_mat[,M1_hat_Q1])
rob_Q1_fit <- lm(y1c ~ Q1_sel_mat - 1)
summary(rob_Q1_fit)

ci_sand_stg1 <- construct_cis_sandwich(rob_Q1_fit) %>% 
  add_length(colnames=c("Lower", "Upper"),
             rownames=colnames(Q1_sel_mat))


set.seed(2022)
blip = function(new_coef) {
  delta <- as.numeric(
    Q2_sel_mat %*% new_coef
  )
  return(pmax(delta, 0) - final_df$A2 * delta)
}
boot_iter_1 <- function() {
  # omega <- rexp(n, rate=1)
  omega <- rpois(n, lambda=1)
  omega <- omega / mean(omega)
  new_coef <- qr.coef(qr(sqrt(omega) * Q2_sel_mat),
                      sqrt(omega) * y2c)
  y1c_new <- y1c + blip(new_coef)
  coef(lm(y1c_new ~ Q1_sel_mat - 1, weights = omega)) - coef(rob_Q1_fit)
}
boots_rob_Q1 <- future_replicate(10000, boot_iter_1())
print("Bootstrap CIs:")
print("------------------")
ci_stg1_boot <- construct_cis_Q(boots_rob_Q1, coef(rob_Q1_fit))
ci_stg1_boot <- add_length(ci_stg1_boot, colnames=c("Lower", "Upper"),
                           rownames=colnames(Q1_sel_mat))
print(ci_stg1_boot)

alpha <- .05
uposi_stg1 <- Uposi$new(
  x=cbind(Q1_mat), y=y1c, stage=1,
  uposi_stage2 = uposi_stg2, a2=final_df$A2, 
  original_blip = blip(coef(rob_Q2_fit)),
  M = M1_hat_Q1, coef=coef(rob_Q1_fit), 
  alpha = alpha, Nboot = 10000, seed = 2022,
  target_type = "fixed"
)
uposi_stg1$do()
Pval1 <- uposi_stg1$pval_H0 # mean(uposi_stg1$boots[1,] <= test_stat_H01)
test_stat_H01 <- uposi_stg1$test_stat_H0
C_alpha_H01 <- uposi_stg1$C_alpha_H0
message(sprintf(paste("Stage: %i",
                      "UPoSI omnibus test statistic: %.3f",
                      "Critical value at alpha=%.2f: %.3f",
                      "UPoSI omnibus P-value: %.3f", sep="\n"), 
                1,
                test_stat_H01, alpha, C_alpha_H01, Pval1))

print("UPOSI CIs:")
print("------------------")
ci_stg1 <- add_length(uposi_stg1$ci_2, colnames=c("Lower", "Upper"),
                      rownames=colnames(Q1_sel_mat))
print(ci_stg1)

if(use_lar){
  print("SI CIs:")
  print("------------------")
  sigmahat <- selectiveInference::estimateSigma(Q1_mat, y1c, standardize=F, intercept=F)$sigmahat
  larinf1 <- selectiveInference::larInf(larfit1, alpha=alpha, type = "all", k = M1_SIZE, sigma=sigmahat)
  # print(larinf1)
  ci_si_stg1 <- larinf1$ci[order(larinf1$vars),]
  ci_si_stg1 <- add_length(ci_si_stg1, colnames=c("Lower", "Upper"),
                           rownames=colnames(Q1_sel_mat))
  print(ci_si_stg1)
}

save.image("./results/run2.RData")

