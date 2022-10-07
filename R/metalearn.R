metalearn <- function(outcome, treat, stage) {
  require(mlr3verse)
  require(mlr3hyperband)
  stopifnot(!is.null(K),
            !is.null(K2),
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
  design_df$outcome <- outcome
  design_df$treat <- factor(treat)
  
  # create Q outcome regression matrix for TMLE
  
  inner_resampling = rsmp("cv", folds = K2)
  search_space = ps(
    nrounds = p_int(lower=50, upper=5000, tags="budget"),
    eta = p_dbl(lower=0.01, upper=0.5),
    lambda = p_dbl(lower=0.1, upper=10, logscale=T),
    alpha = p_dbl(lower=0.01, upper=1, logscale=T),
    # colsample_by_tree = p_dbl(0.5, 1),
    subsample = p_dbl(.8, 1)
  )
  now <- Sys.time()
  stoptime <- now + 4*60*60 # 4 hours
  terminator = trm("combo", terminators=list(
    # trm("evals", n_evals = 5e3),
    trm("clock_time", stop_time=stoptime),
    # trm("stagnation", threshold=0.05, iters=20)
    trm("stagnation_batch", threshold=0.001, n=2)
  ))
  
  # tuner = tnr("grid_search", resolution = 50)
  tuner = tnr("hyperband", eta=2)
  
  Y_task <- mlr3::as_task_regr(design_df, target = "outcome")
  Y_task$set_col_roles("treat", remove_from="feature")
  A_task <- mlr3::as_task_classif(design_df, target = "treat")
  A_task$set_col_roles("outcome", remove_from="feature")
  
  resampling_Y <- mlr3::rsmp("custom_cv")
  resampling_A <- mlr3::rsmp("custom_cv")
  resampling_Y$instantiate(Y_task, f = as.factor(folds_mlr))
  resampling_A$instantiate(A_task, f = as.factor(folds_mlr))
  
  regr_lrn <- mlr3::lrn("regr.xgboost")
  at_regr = AutoTuner$new(regr_lrn, inner_resampling, terminator=terminator, 
                          tuner=tuner, search_space=search_space)
  # encode_po <- po("encode", method="one-hot")
  # graph_Y <- encode_po %>>% at_regr
  
  start <- Sys.time()
  rrY <- mlr3::resample(Y_task, at_regr, resampling_Y, store_models = T,
                        allow_hotstart = T)
  mu.y.h <- rrY$prediction() %>% as.data.table %>% arrange(row_ids) %>% 
    pull(response)
  stop <- Sys.time()
  message(sprintf("Finished learning the outcome in %.2f minutes.", 
                  difftime(stop, start, units="mins")))
  
  # cf_msr <- msr("classif.logloss")
  # cf_lrn <- mlr3::lrn("classif.xgboost", predict_type="prob")
  # cf_lrn <- mlr3::lrn("classif.log_reg", predict_type="prob")
  # at_cf = AutoTuner$new(cf_lrn, inner_resampling, measure=cf_msr, 
  #                       terminator=terminator, 
  #                       tuner=tuner, search_space=search_space)
  # graph_A <- encode_po %>>% cf_lrn
  
  # start <- Sys.time()
  # rrA <- mlr3::resample(A_task, cf_lrn, resampling_A, store_models = T,
  #                       allow_hotstart = T)
  # mu.a.h <- rrA$prediction() %>% as.data.table %>% arrange(row_ids) %>% 
  #   pull(prob.1)
  # stop <- Sys.time()
  # message(sprintf("Finished learning the treatment in %.2f minutes.", 
  #                 difftime(stop, start, units="mins")))
  
  mu_a_hat <- vector("list", K)
  for( k in 1:K ){
    if(stage == 2) {
      mu.a.hat <- glm(A2 ~ response_ind, data=final_df, subset=do.call(c, folds[-k]))
    } else {
      mu.a.hat <- glm(A1 ~ 1, data=final_df, subset=do.call(c, folds[-k]))
    }
    
    mu_a_hat[[k]] <- predict(mu.a.hat, newdata = final_df[folds[[k]],], type="response")
  }
  mu.a.h <- do.call(c, mu_a_hat)
  
  out <- list(
    outcome = mu.y.h,
    treat = mu.a.h,
    outcome_rr = rrY,
    treat_rr = NULL
  )
  class(out) <- "metalearner"
  
  return(out)
}


folds_mlr <- rep(0, n2)
for(k in 1:K){
  folds_mlr[folds[[k]]] <- k
}