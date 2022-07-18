# rql_data.R
source("./R/data-read.R", echo=T)

library(SuperLearner)
library(earth)
library(mgcv)
library(kernlab)
library(randomForest)
library(nloptr)

.define_sl = function(){
  
  create_rbf <- create.Learner("SL.ksvm",
                               tune=list(epsilon=c(.1, .01)),
                               name_prefix = "SL.rbfsvm")
  
  create_tanhsvm <- create.Learner("SL.ksvm",
                                   params = list(kernel = "tanhdot",
                                                 cache = 100),
                                   tune=list(epsilon=c(.1, .01)),
                                   name_prefix = "SL.tanhsvm")
  
  screenrf <- create.Learner("screen.randomForest", 
                             params=list(ntree=30L, nVar=20L))
  
  cont_learners <- c(
    # "SL.randomForest", 
    "SL.lm",
    "SL.glmnet",
    "SL.earth",
    # create_mgcv$names,
    create_rbf$names
    # create_tanhsvm$names,
    # "SL.mean"
  )
  
  bin_learners <- c(
    # "SL.randomForest", 
    "SL.glm",
    "SL.glmnet",
    "SL.earth",
    # create_mgcv$names,
    create_rbf$names
    # create_tanhsvm$names,
    # "SL.mean"
  )
  
  sl_lib <- vector("list", 2L)
  sl_lib$cont_lib <- lapply(cont_learners, function(x) c(x, screenrf$names[1]))
  sl_lib$bin_lib <- lapply(bin_learners, function(x) c(x, screenrf$names[1]))
  
  return(sl_lib)
}
# sl_lib <- .define_sl()
create_rbf <- create.Learner("SL.ksvm",
                             tune=list(epsilon=c(.1, .01)),
                             name_prefix = "SL.rbfsvm")
sl_lib <- list(cont_lib=c(
  "SL.ranger",
  "SL.lm",
  "SL.glmnet",
  "SL.earth",
  # create_mgcv$names,
  create_rbf$names,
  # create_tanhsvm$names,
  "SL.mean"
))

metalearn = function(outcome, stage) {
  require(SuperLearner)
  stopifnot(!is.null(folds), 
            !is.null(inner_folds),
            !is.null(K),
            !is.null(final_df))
  
  switch (as.character(stage),
          "2" = {
            design_df <- final_df %>% 
              select(all_of(h2))
          },
          "1" = {
            design_df <- final_df %>% 
              select(all_of(h1))
          },
          stop("Error: Stage must be 1 or 2.")
  )
  
  SL.CV.control <- list(V=K, validRows=folds, shuffle=FALSE)
  SL.inner.control <- list(list(V=5, validRows=inner_folds, shuffle=FALSE))
  mu.y.hat <- CV.SuperLearner(outcome, design_df, 
                              SL.library=sl_lib$cont_lib,
                              cvControl = SL.CV.control,
                              innerCvControl = SL.inner.control
  )
  mu.y.h <- mu.y.hat$SL.predict
  # mu.a.hat <- CV.SuperLearner(treat, design_df, 
  #                             SL.library=sl_lib$bin_lib,
  #                             family="binomial",
  #                             method="method.CC_nloglik",
  #                             cvControl = SL.CV.control,
  #                             innerCvControl = SL.inner.control
  # )
  # mu.a.h <- mu.a.hat$SL.predict
  
  mu_a_hat <- vector("list", K)
  for( k in 1:K ){
    if(stage == 2) {
      mu.a.hat <- glm(A2 ~ response_ind*gender, data=final_df, subset=do.call(c, folds[-k]))
    } else {
      mu.a.hat <- glm(A1 ~ gender, data=final_df, subset=do.call(c, folds[-k]))
    }
    
    mu_a_hat[[k]] <- predict(mu.a.hat, newdata = final_df[folds[[k]],], type="response")
  }
  mu.a.h <- do.call(c, mu_a_hat)
  
  return(list(
    outcome = mu.y.h,
    treat = mu.a.h,
    outcome_coefs = coef(mu.y.hat),
    treat_coefs = coef(mu.a.hat)
  ))
}

# create folds
K <- 10
K2 <- 5
n2 <- nrow(final_df)
folds <- lapply(1:K, function(v) ((v-1)*n2/K)+(1:(n2/K)))
stopifnot(n2 %% (K*K2) == 0) # not dealing with all of these unequal folds
inner_folds <- lapply(1:K2, function(v) ((v-1)*n2*(K-1)/K/K2)+(1:(n2*(K-1)/K/K2)))

sl_fit2 <- metalearn(outcome = final_df$y, stage = 2)
y2c <- final_df$y - sl_fit2$outcome
a2c <- final_df$A2 - sl_fit2$treat


h2_form <- as.formula(paste0("y2c ~ (a2c*response_ind)*(", paste0(h2, collapse=" + "), ")^2 - 1"))
Q2_mat <- model.matrix(as.formula(h2_form), data=final_df)


### Begin Seletion
library(glmnet)
pen_fac <- rep(1, ncol(Q2_mat))
pen_fac[which(colnames(Q2_mat) %in% c("a2c", "response_ind", "a2c:response_ind"))] <- 0
Q2_sel_fit <- glmnet(Q2_mat, y2c, penalty.factor = pen_fac)
#TODO: switch to ncvreg, or just do a best-subsets
Q2_sel_fit <- ncvreg(Q2_mat, y2c, penalty.factor = pen_fac, dfmax=20)
plot(Q2_sel_fit)
coefs_pen <- coef(Q2_sel_fit, which=24)
M2_hat <- which(abs(coefs_pen) > 1e-5)
M2_hat_Q2 <- M2_hat[-1] - 1

colnames(Q2_mat)[]


plot(Q2_sel_fit)
which(Q2_sel_fit$df == 6)
coefs <- coef(Q2_sel_fit, s=Q2_sel_fit$lambda[which(Q2_sel_fit$df == 6)])
M2_hat <- which(abs(coefs) > 1e-5)
rownames(coefs)[M2_hat]

M2_hat_Q2 <- M2_hat[-1] - 1

# Q2_sel_mat <- cbind(1, scale(Q2_mat[,M2_hat_Q2], F, T))
Q2_sel_mat <- scale(cbind(1, Q2_mat[,M2_hat_Q2]), F, T)
rob_Q2_fit <- lm(y2c ~ Q2_sel_mat - 1)

### Begin UPoSI
N_boot <- 5000
omegas <- replicate(N_boot, rexp(nrow(Q2_mat)))
boot_lm_2 <- function(idx) {
  e <- omegas[,idx]
  max(abs(crossprod(y2c * (e-1), Q2_sel_mat)))
}
boot_uposi <- sapply(1:N_boot, boot_lm_2)
C_alpha_H0 <- quantile(boot_uposi, probs=c(.95))
test_stat_H0 <- max(abs(crossprod(y2c, Q2_sel_mat)))
mean(test_stat_H0 >= boot_uposi)


