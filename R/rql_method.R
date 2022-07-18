require(glmnet)
require(SuperLearner)
require(earth)
require(mgcv)
require(kernlab)
require(randomForest)
library(nloptr)
RQLMethod <- R6::R6Class(
  "RQLMethod",
  cloneable = FALSE,
  portable = FALSE,
  class = FALSE,
  lock_objects = FALSE,
  
  public = list(
    sim = NULL,
    gram = NULL,
    
    mu_method = NULL,
    mu_methods = c("true", "sl"),
    lambda_shrink = NULL,
    cont_lib = NULL,
    bin_lib = NULL,
    K = NULL,
    K2 = NULL,
    alpha=NULL,
    intercept=NULL,
    eps_seed=NULL,
    boot_num=NULL,
    boot_seed=NULL,
    mc_sim_num=NULL,
    mc_integrate=NULL,
    
    # calculated
    folds = NULL,
    inner_folds = NULL,
    Y1_Q = NULL,
    hats2 = NULL,
    hats1 = NULL,
    
    # the history matrices will include intercept
    Y2 = NULL,
    Y1 = NULL,
    H2 = NULL,
    H1 = NULL,
    a1 = NULL,
    a2 = NULL,
    a2c = NULL,
    a1c = NULL,
    y2c = NULL,
    y1c = NULL,
    y1fix = NULL,
    y1fix0 = NULL,
    
    H2cs = NULL,
    H1cs = NULL,
    
    Delta1 = NULL,
    mu_1y_a2_rule = NULL,
    xi_hat = NULL,
    xi_star = NULL,
    est_a2_rule = NULL,
    est_delta2_fun = NULL,
    
    
    # only responders in the second stage
    n2 = NULL,
    p2 = NULL,
    n1 = NULL,
    p1 = NULL,
    
    M2 = NULL,
    M1 = NULL,
    beta2_M = NULL,
    beta1_M = NULL,
    
    beta20_M = NULL,
    beta10_M = NULL,
    
    wgt_coef2 = NULL,
    beta2_lasso = NULL,
    wgt_coef1 = NULL,
    beta1_lasso = NULL,
    wgt_stg2 = NULL,
    wgt_stg1 = NULL,
    lambda1 = NULL,
    lambda2 = NULL,
    
    outputStg2 = NULL,
    outputPredStg2 = NULL,
    outputStg1 = NULL,
    outputPredStg1 = NULL,
    
    uposi_stg2 = NULL,
    uposi_stg2_wgt = NULL,
    uposi_stg2_fix = NULL,
    si_stg2 = NULL,
    uposi_stg1 = NULL,
    uposi_stg1_wgt = NULL,
    uposi_stg1_fix = NULL,
    si_stg1 = NULL,
    
    initialize = function(
      sim, 
      # gram,
      alpha, intercept, eps_seed, mc_sim_num,
      boot_seed,
      boot_num, K,
      mu_method, lambda_shrink=2
    ) {
      stopifnot(
        inherits(sim, "Simulation"),
        # is.null(gram) || inherits(gram, "Gram"),
        mu_method %in% self$mu_methods,
        is.numeric(alpha) && length(alpha) == 1 && 0 < alpha && alpha < 1,
        is.logical(intercept) || (
          is.numeric(intercept) && intercept %in% 0:1
        ) && length(intercept) == 1,
        is.numeric(eps_seed) && length(eps_seed) == 1,
        is.numeric(mc_sim_num) && length(mc_sim_num) == 1,
        is.numeric(boot_seed) && length(boot_seed) == 1
        # is.null(cont_lib) || (is.character(cont_lib) && length(cont_lib) > 0),
        # is.null(bin_lib) || (is.character(bin_lib) && length(bin_lib) > 0)
        # ,
        # is.logical(mc_integrate) || (
        #   is.numeric(mc_integrate) && mc_integrate %in% 0:1
        # ) && length(mc_integrate) == 1
      )
      self$sim = sim
      self$gram = gram
      self$mu_method = mu_method
      self$lambda_shrink= lambda_shrink
      self$alpha = alpha
      self$intercept = as.logical(intercept)
      self$eps_seed = as.integer(eps_seed)
      self$boot_seed = as.integer(boot_seed)
      self$boot_num = as.integer(boot_num)
      self$mc_sim_num = as.integer(mc_sim_num)
      self$K = K
      
      self$.define_sl()
      invisible(self)
    },
    
    apply_stg2 = function(inf_method=NULL) {
      if("posi-wgt" %in% inf_method & !("posi" %in% inf_method))
        stop("At this time, weighted UPoSI must be run with un-weighted.")
      
      nr <<- which(self$sim$sim_data$nr==1)
      n2 <<- length(nr)
      X1 <- unname(self$sim$sim_data$X1[nr,])
      X2 <- unname(self$sim$sim_data$X2[nr,])
      Y2 <<- self$sim$sim_data$y[nr]
      a2 <<- self$sim$sim_data$a2[nr]
      a1 <<- self$sim$sim_data$a1[nr]
      
      # H2.df <- data.frame(X1=X1, X2=X2, a1=a1)
      
      H2 <<- cbind(
        1,
        self$sim$sim_data$X1,
        self$sim$sim_data$X2,
        self$sim$sim_data$a1
      )[nr,]
      colnames(self$H2) <- c(
        "(Intercept)", paste0("X1_", 1:ncol(self$sim$sim_data$X1)),
        paste0("X2_", 1:ncol(self$sim$sim_data$X2)), "a1"
      )
      
      stopifnot(is.matrix(H2), nrow(H2) == self$n2)
      p2 <<- ncol(H2)
      
      # stage 1 init
      n1 <<- self$sim$n
      H1 <<- cbind(1, self$sim$sim_data$X1)
      stopifnot(is.matrix(H1), nrow(H1) == self$n1)
      p1 <<- ncol(H1)
      stopifnot(p1*2 == p2)
      colnames(H1) <<- colnames(H2)[1:p1]
      
      # create folds
      folds <<- lapply(1:K, function(v) ((v-1)*n2/K)+(1:(n2/K)))
      # TODO: Adjust this for non-response
      K2 <<- 5
      stopifnot(n2 %% (K*K2) == 0) # not dealing with all of these unequal folds
      inner_folds <<- lapply(1:K2, function(v) ((v-1)*n2*(K-1)/K/K2)+(1:(n2*(K-1)/K/K2)))
      
      hats <- self$center(outcome=Y2, treat=a2, stage=2)
      # hats <- self$metalearn(outcome=Y, treat=a2, stage=2)
      mu.y.h <- as.numeric(hats[["outcome"]])
      mu.a2.h <- as.numeric(hats[["treat"]])
      stopifnot(length(mu.y.h) == self$n2, length(mu.a2.h) == self$n2)
      self$hats2 <- hats
      hats <- NULL
      
      y2c <<- as.numeric(Y2 - mu.y.h)
      # y2c <<- y2c
      a2c <<- as.numeric(a2 - mu.a2.h)
      H2c <<- private$row_mult(H2[,-1], a2c)
      # H2c <<- scale(H2c, scale=F, center=T)
      H2c <<- cbind(a2c, H2c)
      colnames(H2c) <<- colnames(H2)
      
      H2c.scale <- apply(H2c, 2, sd)
      H2c.scale <- sapply(H2c.scale, function(s) ifelse(s==0, 1, s))
      stopifnot(length(H2c.scale) == p2)
      
      H2cs <<- scale(H2c, center=F, scale=H2c.scale)
      
      selection_out <- self$.use_selection_on_stage(2)
      beta2_lasso <<- selection_out$beta_lasso
      wgt_coef <- selection_out$wgt_coef
      coef_opts <- selection_out$coef_opts
      wgt_coef2 <<- wgt_coef
      # wgt_stg2 <<- pexp(abs(wgt_coef))
      wgt_stg2 <<- 1*(abs(wgt_coef) > 1e-5)
      lambda2 <<- coef_opts$s
      
      # selectiveInference is using some heuristic fudging here,
      # which we reproduce for comparability
      M2 <<- which(abs(beta2_lasso) > 1e-5 / sqrt(colSums(H2cs^2)))
      
      beta2_M <<- rep(0, p2)
      rob_Q2_fit <- lm(y2c ~ H2c[, M2, drop = F] - 1)
      beta2_M[M2] <<- unname(coef(rob_Q2_fit))
      
      
      ##### INFERENCE #####
      outputList <- outputPredictList <- vector("list", 10L)
      
      naive_ci = confint(rob_Q2_fit, level=1-alpha) # get rid of grand mean
      naive.lower = c(naive_ci[,1])
      naive.upper = c(naive_ci[,2])
      naive.f <- summary(rob_Q2_fit)$fstatistic
      naive.p <- 1-pf(naive.f[1], naive.f[2], naive.f[3])
      # browser()
      stage2df <- data.frame(
        sim = self$sim$sim_idx,
        stage = 2,
        alpha = self$alpha,
        var=M2,
        betahat = beta2_M[M2],
        mu_a_error = sqrt(mean((mu.a2.h - self$sim$sim_data$mu.a2)^2)),
        mu_a_bias = mean(mu.a2.h - self$sim$sim_data$mu.a2),
        mu_y_error = sqrt(mean((mu.y.h - self$sim$sim_data$mu.y2x)^2)),
        mu_y_bias = mean(mu.y.h - self$sim$sim_data$mu.y2x),
        pseudo.method = "NA"
      )
      rownames(stage2df) <- colnames(H2c)[M2]
      
      outputList[[1]] <- data.frame(
        stage2df,
        inf_method="naive",
        lower = naive.lower, upper = naive.upper,
        overall_p = naive.p,
        ci_len = naive.upper - naive.lower
      )
      out_pred_idx <- output_idx <- 1
      if("posi" %in% tolower(inf_method)){
        output_idx <- output_idx + 1
        uposi_stg2 <<- UposiRandom$new(
          H2c, y2c, M=M2, coef=beta2_M[M2],
          alpha = self$alpha, seed = self$boot_seed)
        box <- uposi_stg2$do(self$boot_num)$ci
        pval <- uposi_stg2$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                               inf_method="PoSI",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
        # now look at the alternative cis
        output_idx <- output_idx + 1
        box <- uposi_stg2$ci_2
        pval <- uposi_stg2$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                               inf_method="PoSI-coordinate",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
      }
      
      if("posi-wgt" %in% tolower(inf_method)){
        output_idx <- output_idx + 1
        uposi_stg2_wgt <<- UposiRandomWeight$new(
          H2c, y2c, M=M2, coef=beta2_M[M2],
          alpha = self$alpha, seed = self$boot_seed,
          weights=wgt_stg2)
        box <- uposi_stg2_wgt$do(self$boot_num)$ci
        pval <- uposi_stg2_wgt$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                               inf_method="PoSI-wgt",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
        # now look at the alternative cis
        output_idx <- output_idx + 1
        box <- uposi_stg2_wgt$ci_2
        pval <- uposi_stg2_wgt$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                               inf_method="PoSI-wgt-coordinate",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
      }
      
      if("posi-fix" %in% tolower(inf_method)){
        output_idx <- output_idx + 1
        uposi_stg2_fix <<- UposiFixed$new(H2c, y2c, M=M2, coef=beta2_M[M2],
                                      alpha = self$alpha, seed = self$boot_seed)
        uposi_stg2_fix$boots <- uposi_stg2$boots[1,] # copy over boots
        uposi_stg2_fix$Nboot <- uposi_stg2$Nboot
        
        box <- uposi_stg2_fix$extract_quantile()$confint()$ci
        pval <- uposi_stg2_fix$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"] 
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                               inf_method="PoSI-fix",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        output_idx <- output_idx + 1
        box <- uposi_stg2_fix$ci_2
        pval <- uposi_stg2_fix$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                               inf_method="PoSI-fix-coordinate",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
      }
      
      if("si" %in% tolower(inf_method)){
        stop("Not fully implemented.")
        output_idx <- output_idx + 1
        out <- selectiveInference$new(
          H2cs, y2c, beta2_lasso, lambda=coef_opts$s,
          scale_vars = H2c.scale, intercept=coef_opts$intercept,
          penalty.factor = coef_opts$penalty.factor, alpha = self$alpha,
          gridpts = 1e4
        )
        # box <- sweep(out$ci, 1, 1/H2c.scale[out$vars], "*")
        box <- out$ci[-1,]
        if(nrow(box) != length(M2))
          browser()
        si_stg2 <<- out
        lower = c(box[,1])
        upper = c(box[,2])
        
        outputList[[output_idx]] <- data.frame(stage2df,
                                      inf_method="SI",
                                      lower = lower,
                                      upper = upper,
                                      overall_p = NA,
                                      ci_len = upper - lower)
        
      }
      
      outputStg2 <<- do.call(rbind, outputList)
      outputPredStg2 <<- do.call(rbind, outputPredictList)
      
      invisible(self)
    },
    
    prepare_blips = function() {
      
      delta2_hat <- H2 %*% beta2_M
      # delta2_star <- H2 %*% beta20_M
      
      Y1_Q <<- as.numeric((1*(delta2_hat > 0) - hats2$treat) * delta2_hat) +
        as.numeric(hats2[["outcome"]])
      
      xi_fun <- function(delta, a2) (1*(delta > 0) - a2) * delta
      
      xi_star <<- xi_hat <<- rep(0, n1)
      xi_hat[nr] <<- xi_fun(delta2_hat, a2)
      # xi_star[nr] <<- xi_fun(delta2_star, a2)
      
      
      ##### Calculate conditional expectations of future trajectories #####
      beta2_rule <- as.numeric(beta2_M)
      if(p2 != length(beta2_rule))
        stop("Error: beta2 is not of the correct size.")
      est_a2_rule <<- function(X2, X1, A1) {
        as.integer(sum(c(1, X1, X2, A1) * beta2_rule) > 0)
      }
      est_delta2_fun <<- function(X2, X1, A1) {
        sum(c(1, X1) * beta2_rule[1:p1]) + 
          A1*beta2_rule[p2] +
          X2 %*% beta2_rule[p1 + (1:(p1 - 1))]
      }
      
      # hard-code use xi_hat
      Y1 <<- sim$sim_data$y + xi_hat
      
      invisible(self)
    },
    
    apply_stg1 = function(inf_method = NULL) {
      if("posi-wgt" %in% inf_method & !("posi" %in% inf_method))
        stop("At this time, weighted UPoSI must be run with un-weighted.")
      
      hats <- self$center(outcome=Y1, treat=a1, stage=1)
      mu_y1_h <- as.numeric(hats[["outcome"]])
      mu_a1_h <- as.numeric(hats[["treat"]])
      stopifnot(length(mu_y1_h) == self$n1, length(mu_a1_h) == self$n1)
      self$hats1 <- hats
      hats <- NULL
      
      # browser()
      
      y1c <<- as.numeric(Y1 - mu_y1_h)
      y1fix <<- sim$sim_data$y - mu_y1_h
      a1c <<- as.numeric(a1 - mu_a1_h)
      H1c <<- private$row_mult(H1, a1c)
      colnames(H1c) <<- colnames(H1)
      
      H1c_scale <- apply(H1c, 2, sd)
      H1c_scale <- sapply(H1c_scale, function(s) ifelse(s==0, 1, s))
      H1cs <<- scale(H1c, center=F, scale=H1c_scale)
      
      
      selection_out <- self$.use_selection_on_stage(1)
      beta1_lasso <<- selection_out$beta_lasso
      wgt_coef <- selection_out$wgt_coef
      coef_opts <- selection_out$coef_opts
      wgt_coef1 <<- wgt_coef
      # wgt_stg1 <<- pexp(abs(wgt_coef))
      wgt_stg1 <<- 1*(abs(wgt_coef) > 1e-5)
      lambda1 <<- coef_opts$s
      
      M1 <<- which(abs(beta1_lasso) > 1e-5 / sqrt(colSums(H1cs^2)))
      
      
      
      beta1_M <<- rep(0, p1)
      rob_Q1_fit <- lm(y1c ~ H1c[, M1, drop = F] - 1)
      beta1_M[M1] <<- unname(coef(rob_Q1_fit))
      
      ##### INFERENCE #####
      outputList <- outputPredictList <- vector("list", 10L)
      
      naive_ci = confint(rob_Q1_fit, level=1-alpha) # get rid of grand mean
      naive.lower = c(naive_ci[,1])
      naive.upper = c(naive_ci[,2])
      naive.f <- summary(rob_Q1_fit)$fstatistic
      naive.p <- 1-pf(naive.f[1], naive.f[2], naive.f[3])
      # browser()
      stage1df <- data.frame(
        sim = self$sim$sim_idx,
        stage = 1,
        alpha = self$alpha,
        var=M1,
        betahat = beta1_M[M1],
        mu_a_error = sqrt(mean((mu_a1_h - self$sim$sim_data$mu.a1)^2)),
        mu_a_bias = mean(mu_a1_h - self$sim$sim_data$mu.a1),
        mu_y_error = NA,
        mu_y_bias = NA,
        pseudo.method = "NA"
      )
      rownames(stage1df) <- colnames(H1c)[M1]
      
      outputList[[1]] <- data.frame(
        stage1df,
        inf_method="naive",
        lower = naive.lower, upper = naive.upper,
        overall_p = naive.p,
        ci_len = naive.upper - naive.lower
      )
      out_pred_idx <- output_idx <- 1
      fixed_implemented <- FALSE
      weights_implemented <- TRUE
      if("posi" %in% tolower(inf_method)){
        output_idx <- output_idx + 1
        uposi_stg1 <<- UposiStage1$new(
          self$uposi_stg2, self$a2, self$xi_hat,
          H1c, y1c, M=M1, coef=beta1_M[M1],
          alpha = self$alpha, seed = self$boot_seed)
        box <- uposi_stg1$do(self$boot_num)$ci
        pval <- uposi_stg1$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="PoSI",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
        # now look at the alternative cis
        output_idx <- output_idx + 1
        box <- uposi_stg1$ci_2
        pval <- uposi_stg1$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="PoSI-coordinate",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
      }
      
      if("posi-wgt" %in% tolower(inf_method) & weights_implemented){
        output_idx <- output_idx + 1
        uposi_stg1_wgt <<- UposiStage1Weight$new(
          self$uposi_stg2, self$a2, self$xi_hat,
          H1c, y1c, M=M1, coef=beta1_M[M1],
          alpha = self$alpha, seed = self$boot_seed,
          weights=wgt_stg1)
        box <- uposi_stg1_wgt$do(self$boot_num)$ci
        pval <- uposi_stg1_wgt$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="PoSI-wgt",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
        # now look at the alternative cis
        output_idx <- output_idx + 1
        box <- uposi_stg1_wgt$ci_2
        pval <- uposi_stg1_wgt$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="PoSI-wgt-coordinate",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        
      }
      
      if("posi-fix" %in% tolower(inf_method) & fixed_implemented){
        output_idx <- output_idx + 1
        uposi_stg1_fix <<- UposiFixed$new(H1c, y1c, M=M1, coef=beta1_M[M1],
                                          alpha = self$alpha, seed = self$boot_seed)
        uposi_stg1_fix$boots <- uposi_stg1$boots[1,] # copy over boots
        uposi_stg1_fix$Nboot <- uposi_stg1$Nboot
        
        box <- uposi_stg1_fix$extract_quantile()$confint()$ci
        pval <- uposi_stg1_fix$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="PoSI-fix",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
        output_idx <- output_idx + 1
        box <- uposi_stg1_fix$ci_2
        pval <- uposi_stg1_fix$pval_H0i
        lower = box[,"Lower"]
        upper = box[,"Upper"]
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="PoSI-fix-coordinate",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = pval,
                                               ci_len = upper - lower)
      }
      
      if("si" %in% tolower(inf_method)){
        stop("Not fully implemented.")
        output_idx <- output_idx + 1
        out <- selectiveInference$new(
          H1cs, y1c, beta2_lasso, lambda=coef_opts$s,
          scale_vars = H1c.scale, intercept=coef_opts$intercept,
          penalty.factor = coef_opts$penalty.factor, alpha = self$alpha,
          gridpts = 1e4
        )
        box <- out$ci[-1,]
        if(nrow(box) != length(M1))
          browser()
        si_stg2 <<- out
        lower = c(box[,1])
        upper = c(box[,2])
        
        outputList[[output_idx]] <- data.frame(stage1df,
                                               inf_method="SI",
                                               lower = lower,
                                               upper = upper,
                                               overall_p = NA,
                                               ci_len = upper - lower)
        
      }
      
      outputStg1 <<- do.call(rbind, outputList)
      outputPredStg1 <<- do.call(rbind, outputPredictList)
      
      
      invisible(self)
    },
    
    apply_to_data = function(inf_method = c("posi", "si")) {
      if("posi-wgt" %in% inf_method & !("posi" %in% inf_method))
        stop("At this time, weighted UPoSI must be run with un-weighted.")
      
      self$apply_stg2(inf_method = inf_method)
      self$prepare_blips()
      self$apply_stg1(inf_method = inf_method)
      invisible(self)
    },
    
    center = function(outcome, treat, stage=2) {
      if(is.matrix(outcome))
        outcome <- as.numeric(outcome)
      if(is.matrix(treat))
        treat <- as.numeric(treat)
      
      if(stage == 2) {
        stopifnot(
          length(outcome) == self$n2,
          length(treat) == self$n2
        )
        
        if(self$mu_method == "true") {
          ret <- list(
            outcome = self$sim$sim_data$mu.y2x,
            treat = self$sim$sim_data$mu.a2
          )
        } else {
          ret <- self$metalearn(outcome, treat, stage)
        }
        
      } else if(stage == 1) {
        stopifnot(
          length(outcome) == self$n1,
          length(treat) == self$n1
        )
        if(self$mu_method == "true") {
          ret <- list(
            outcome = self$mu_1y_a2_rule,
            treat = self$sim$sim_data$mu.a1
          )
        } else {
          ret <- self$metalearn(outcome, treat, stage)
        }
      } else {
        stop("Error: Stage must be 1 or 2.")
      }
      return(ret)
    },
    
    metalearn = function(outcome, treat, stage) {
      require(SuperLearner)
      stopifnot(!is.null(folds), 
                !is.null(inner_folds),
                !is.null(K),
                !is.null(H1),
                !is.null(H2))
      
      switch (as.character(stage),
        "2" = {
          design_df <- as.data.frame(H2[,-1])
        },
        "1" = {
          design_df <- as.data.frame(H1[,-1])
        },
        stop("Error: Stage must be 1 or 2.")
      )
      
      SL.CV.control <- list(V=K, validRows=folds, shuffle=FALSE)
      SL.inner.control <- list(list(V=5, validRows=inner_folds, shuffle=FALSE))
      mu.y.hat <- CV.SuperLearner(outcome, design_df, 
                                  SL.library=self$cont_lib,
                                  cvControl = SL.CV.control,
                                  innerCvControl = SL.inner.control
      )
      mu.y.h <- mu.y.hat$SL.predict
      mu.a.hat <- CV.SuperLearner(treat, design_df, 
                                  SL.library=self$bin_lib,
                                  family="binomial",
                                  method="method.CC_nloglik",
                                  cvControl = SL.CV.control,
                                  innerCvControl = SL.inner.control
      )
      mu.a.h <- mu.a.hat$SL.predict
      
      return(list(
        outcome = mu.y.h,
        treat = mu.a.h,
        outcome_coefs = coef(mu.y.hat),
        treat_coefs = coef(mu.a.hat)
      ))
    },
    
    .define_sl = function(){
      
      create_rbf <- create.Learner("SL.ksvm",
                                   tune=list(epsilon=c(.1, .01)),
                                   env = self,
                                   name_prefix = "SL.rbfsvm")
      
      create_tanhsvm <- create.Learner("SL.ksvm",
                                       params = list(kernel = "tanhdot",
                                                     cache = 100),
                                       tune=list(epsilon=c(.1, .01)),
                                       env = self,
                                       name_prefix = "SL.tanhsvm")
      
      screenrf <- create.Learner("screen.randomForest", 
                                 params=list(ntree=30L, nVar=20L),
                                 env = self)
      
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
      
      self$cont_lib <- lapply(cont_learners, function(x) c(x, screenrf$names[1]))
      self$bin_lib <- lapply(bin_learners, function(x) c(x, screenrf$names[1]))
      
      invisible(self)
    },
    
    .use_selection_on_stage = function(stage) {
      
      if(stage == 2) {
        p <- p2
        n <- n2
        design_cs <- H2cs
        outcome_c <- y2c
        
      } else if(stage == 1) {
        p <- p1
        n <- n1
        design_cs <- H1cs
        outcome_c <- y1c
        xi_fun <- function(delta, a2) (1*(delta > 0) - a2) * delta
      } else {
        stop("Invalid stage!")
      }
      
      pen <- c(0, rep(1,p-1)) # don't penalize ate
      opts <- list(
        x=design_cs, y=outcome_c, standardize=F, intercept = T, thresh=1e-14,
        penalty.factor = pen
      )
      
      glmnet_folds <- rep(0, n)
      for(k in 1:K) {
        glmnet_folds[folds[[k]]] <- k
      }
      opts$foldid <- glmnet_folds
      gfit <- do.call(cv.glmnet, opts)
      lambda <- gfit$lambda.min
      coef_opts <- c(list(object=gfit, s="lambda.min"), opts)
      
      # if(!is.null(eps_seed)){
      #   set.seed(eps_seed)
      # }
      # 
      # eps <- matrix(rnorm(n * boot_num, sd = 1),
      #               ncol = boot_num)
      # lambda <- lambda_shrink * mean(
      #   apply(abs(crossprod(design_cs, eps)), 2, max)
      # )/n
      # gfit <- do.call(glmnet, opts)
      # coef_opts <- c(list(object=gfit, s=lambda), opts)
      # beta_lasso <- as.numeric(do.call(coef, coef_opts))
      # # 
      # # # that's just the preliminary lambda; get the real one
      # # if(stage == 1) {
      # #   browser()
      # #   Mprelim <- which(abs(beta_lasso) > 1e-5)
      # # 
      # #   # eps <- matrix(0, ncol=boot_num, nrow=n)
      # #   # for(col_idx in 1:boot_num) {
      # #   #   e <- rexp(n)
      # #   # 
      # #   #   beta_prelim <- rep(0, p2)
      # #   #   beta_prelim[Mprelim] <- qr.coef(qr(private$row_mult(H2c[,Mprelim,drop=F], e)), outcome_c)
      # #   # 
      # #   #   blip_err <- xi_fun(as.numeric(H2 %*% beta_prelim), a2) - xi_hat
      # #   # 
      # #   #   eps[,col_idx] <- e * (outcome_c + blip_err - (design_cs %*% beta_lasso[-1] + beta_lasso[1]))
      # #   # }
      # #   # beta_prelim <- rep(0, p1+1)
      # #   # beta_prelim[Mprelim] <- qr.coef(qr(cbind(1,H1c)[,Mprelim,drop=F]), y1c)
      # #   # Delta1_prelim <- as.numeric(beta_prelim[1] + H1 %*% beta_prelim[-1])
      # #   Delta1_prelim <- as.numeric(beta_lasso[1] + H1 %*% beta_lasso[-1])
      # #   eps <- matrix(0, ncol=boot_num, nrow=n)
      # #   for(col_idx in 1:boot_num) {
      # #     e <- rnorm(n)
      # #     
      # #     # beta2_pert <- rep(0, p2)
      # #     # beta2_pert[Mprelim] <- qr.coef(
      # #     #   qr(private$row_mult(H2c[,M2,drop=F], sqrt(e))), 
      # #     #   sqrt(e) * y2c)
      # #     # blip_err <- xi_fun(as.numeric(H2 %*% beta2_pert), a2) - xi_hat
      # # 
      # #     # eps[,col_idx] <- e * (outcome_c + blip_err - (Delta1))
      # #     eps[,col_idx] <- e * (outcome_c - (Delta1_prelim))
      # #   }
      # # } else {
      #   eps <- matrix(rnorm(
      #     n * boot_num,
      #     sd = abs(as.numeric(outcome_c - (design_cs %*% beta_lasso[-1] + beta_lasso[1])))
      #   ), ncol = boot_num)
      # # }
      # lambda <- lambda_shrink * median(
      #   apply(abs(crossprod(design_cs, eps)), 2, max)
      # )/n

      coef_opts$s <- lambda
      coef_opts$exact <- TRUE
      beta_lasso <- as.numeric(do.call(coef, coef_opts))[-1]
      
      wgt_opts <- coef_opts
      wgt_opts$s <- lambda/(5)
      wgt_coef <<- as.numeric(do.call(coef, wgt_opts))[-1]
      
      return(list(beta_lasso = beta_lasso, wgt_coef = wgt_coef,
                  coef_opts = coef_opts))
    },
    
    integrate_functionals = function(input_a2_rule, input_delta2_fun, save=F) {
      
      val.functionals <- foreach(
        row_idx=seq.int(nrow(self$H1)), .combine=cbind
      ) %dopar% {
        self$sim$mu_y1_x1_integrate(self$H1[row_idx,-1],
                                    mc_num_sims=self$mc_sim_num,
                                    a2_rule=input_a2_rule, 
                                    delta2_fun=input_delta2_fun)
      }
      
      Delta1 <- val.functionals[2,]
      mu_1y_a2_rule <- val.functionals[1,]
      if(save){
        Delta1 <<- Delta1
        mu_1y_a2_rule <<- mu_1y_a2_rule
      }
      return(list(Delta1=Delta1, mu_1y_a2_rule=mu_1y_a2_rule))
    },
    
    calculate_coverages = function() {
      stopifnot(!is.null(beta10_M), !is.null(beta20_M))
      
      outputStg1 <- self$outputStg1
      stopifnot(nrow(outputStg1) %% length(beta10_M) == 0)
      beta10_rep <- rep(beta10_M, nrow(outputStg1) / length(beta10_M))
      stopifnot(length(beta10_rep) == nrow(outputStg1))
      
      outputStg1$beta0 <- beta10_rep
      outputStg1$cover <- 1*(outputStg1$lower <=beta10_rep & 
                               beta10_rep <= outputStg1$upper)
      self$outputStg1 <- outputStg1
      
      outputStg2 <- self$outputStg2
      stopifnot(nrow(outputStg2) %% length(beta20_M) == 0)
      beta20_rep <- rep(beta20_M, nrow(outputStg2) / length(beta20_M))
      stopifnot(length(beta20_rep) == nrow(outputStg2))
      
      outputStg2$beta0 <- beta20_rep
      outputStg2$cover <- 1*(outputStg2$lower <=beta20_rep & 
                               beta20_rep <= outputStg2$upper)
      self$outputStg2 <- outputStg2
      self$calculate_prediction_covs()
      invisible(self)
      
    },
    
    calculate_prediction_covs = function() {
      stopifnot(!is.null(beta10_M), !is.null(beta20_M))
      outputPredictList1 <- vector("list", 10L)
      outputPredictList2 <- vector("list", 10L)
      out_pred_idx <- 0
      Delta0 <- as.numeric(H1[,M1,drop=F] %*% beta10_M)
      if(!is.null(uposi_stg1)){
        out_pred_idx <- out_pred_idx + 1
        Delt_int <- uposi_stg1$pred_int(H1[,M1,drop=F])
        Delta_cover <- sapply(seq.int(n1), function(i) {
          1*(Delt_int[1, i] <= Delta0[i] & Delta0[i] <= Delt_int[2, i])
        })
        Delta_cover_all <- prod(Delta_cover)
        Delta_cover_rate <- mean(Delta_cover)
        outputPredictList1[[out_pred_idx]] <- data.frame(
          stage=1,
          alpha = self$alpha,
          inf_method="PoSI",
          cover_rate=Delta_cover_rate,
          cover_all=Delta_cover_all
        )
      }
      if(!is.null(uposi_stg1_fix)){
        out_pred_idx <- out_pred_idx + 1
        Delt_int <- uposi_stg1_fix$pred_int(H1[,M1,drop=F])
        Delta_cover <- sapply(seq.int(n1), function(i) {
          1*(Delt_int[1, i] <= Delta0[i] & Delta0[i] <= Delt_int[2, i])
        })
        Delta_cover_all <- prod(Delta_cover)
        Delta_cover_rate <- mean(Delta_cover)
        outputPredictList1[[out_pred_idx]] <- data.frame(
          stage=1,
          alpha = self$alpha,
          inf_method="PoSI-fix",
          cover_rate=Delta_cover_rate,
          cover_all=Delta_cover_all
        )
      }
      if(!is.null(uposi_stg1_wgt)){
        out_pred_idx <- out_pred_idx + 1
        Delt_int <- uposi_stg1_wgt$pred_int(H1[,M1,drop=F])
        Delta_cover <- sapply(seq.int(n1), function(i) {
          1*(Delt_int[1, i] <= Delta0[i] & Delta0[i] <= Delt_int[2, i])
        })
        Delta_cover_all <- prod(Delta_cover)
        Delta_cover_rate <- mean(Delta_cover)
        outputPredictList1[[out_pred_idx]] <- data.frame(
          stage=1,
          alpha = self$alpha,
          inf_method="PoSI-wgt",
          cover_rate=Delta_cover_rate,
          cover_all=Delta_cover_all
        )
      }
      
      # Stage 2
      out_pred_idx <- 0
      Delta0 <- as.numeric(H2[,M2,drop=F] %*% beta20_M)
      if(!is.null(uposi_stg2)){
        out_pred_idx <- out_pred_idx + 1
        Delt_int <- uposi_stg2$pred_int(H2[,M2,drop=F])
        Delta_cover <- sapply(seq.int(n1), function(i) {
          1*(Delt_int[1, i] <= Delta0[i] & Delta0[i] <= Delt_int[2, i])
        })
        Delta_cover_all <- prod(Delta_cover)
        Delta_cover_rate <- mean(Delta_cover)
        outputPredictList2[[out_pred_idx]] <- data.frame(
          stage=2,
          alpha = self$alpha,
          inf_method="PoSI",
          cover_rate=Delta_cover_rate,
          cover_all=Delta_cover_all
        )
      }
      if(!is.null(uposi_stg2_fix)){
        out_pred_idx <- out_pred_idx + 1
        Delt_int <- uposi_stg2_fix$pred_int(H2[,M2,drop=F])
        Delta_cover <- sapply(seq.int(n1), function(i) {
          1*(Delt_int[1, i] <= Delta0[i] & Delta0[i] <= Delt_int[2, i])
        })
        Delta_cover_all <- prod(Delta_cover)
        Delta_cover_rate <- mean(Delta_cover)
        outputPredictList2[[out_pred_idx]] <- data.frame(
          stage=2,
          alpha = self$alpha,
          inf_method="PoSI-fix",
          cover_rate=Delta_cover_rate,
          cover_all=Delta_cover_all
        )
      }
      if(!is.null(uposi_stg2_wgt)){
        out_pred_idx <- out_pred_idx + 1
        Delt_int <- uposi_stg2_wgt$pred_int(H2[,M2,drop=F])
        Delta_cover <- sapply(seq.int(n1), function(i) {
          1*(Delt_int[1, i] <= Delta0[i] & Delta0[i] <= Delt_int[2, i])
        })
        Delta_cover_all <- prod(Delta_cover)
        Delta_cover_rate <- mean(Delta_cover)
        outputPredictList2[[out_pred_idx]] <- data.frame(
          stage=2,
          alpha = self$alpha,
          inf_method="PoSI-wgt",
          cover_rate=Delta_cover_rate,
          cover_all=Delta_cover_all
        )
      }
      
      outputPredStg1 <<- do.call(rbind, outputPredictList1)
      outputPredStg2 <<- do.call(rbind, outputPredictList2)
      invisible(self)
    }
  ),
  
  
  private = list(
    
    #' Helper: Row multiplication
    #'
    #' Multiply the rows of a matrix x by a vector or column matrix y.
    #' This is equivalent to (but faster than) diag(y) %*% x.
    #' @param x a matrix
    #' @param y a vector
    #'
    #' @return diag(y) %*% x, but faster
    #'
    row_mult = function(x, y) {
      # return(sweep(x, 1, y, "*"))
      return(x*y[row(x)])
    }
    
  )
)

