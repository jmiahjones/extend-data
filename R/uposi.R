library(R6)
library(checkmate)
library(future.apply)

Uposi <- R6::R6Class(
  "Uposi",
  public = list(
    
    x = NULL,
    y = NULL,
    n = NULL,
    p = NULL,
    
    stage = NULL,
    
    alpha = NULL,
    Nboot = NULL,
    seed = NULL,
    
    M = NULL,
    coef = NULL,
    target_type = NULL,
    
    standardize = NULL,
    intercept = NULL,
    center_outcome = NULL,
    
    boots = NULL,
    C_alpha = NULL,
    C_alpha_H0 = NULL,
    
    Sigma = NULL,
    Sigmainv = NULL,
    
    ci = NULL,
    ci_2 = NULL,
    pred_i = NULL,
    
    test_stat_H0 = NULL,
    pval_H0 = NULL,
    
    uposi_stage2 = NULL,
    a2 = NULL,
    original_blip = NULL,
    
    initialize = function(
      x, y, stage, M, coef=NULL, 
      uposi_stage2 = NULL, a2 = NULL, original_blip = NULL,
      alpha = 0.05,
      Nboot = 1000, seed = NULL,
      target_type = "fixed",
      standardize = T,
      intercept = F,
      center_outcome = F
    ) {
      # TODO: How to handle intercept? Ensure X doesn't have the constant
      # column? Translate M to intercept numbering?
      
      
      assert_matrix(x, any.missing = F)
      self$x = x
      self$n = nrow(x)
      self$p = ncol(x)
      if(test_numeric(y, any.missing = F, len = self$n)) {
        y <- matrix(y, ncol=1)
      }
      assert_matrix(y, any.missing = F, nrows = self$n, ncols = 1)
      self$y = y
      
      stage <- assert_count(stage, coerce=T, positive = T)
      assert_integer(stage, upper=2L)
      self$stage = stage
      
      M <- assert_integerish(M, lower=1, upper=self$p, 
                             sorted=T, unique=T, coerce=T)
      self$M = M
      assert(
        check_numeric(coef, len=length(M)),
        check_matrix(coef, nrows=length(M), ncols=1),
        check_null(coef)
      )
      if(is.null(coef)) {
        coef = qr.coef(qr(self$x[,M]), self$y)
      }
      self$coef = matrix(coef, ncol=1)
      
      self$Sigma = crossprod(self$x[,self$M])/self$n
      
      assert_double(alpha, lower = 0, upper = 1, any.missing = F, len = 1)
      self$alpha = alpha
      
      Nboot = assert_count(Nboot, positive=T, coerce=T)
      self$Nboot = Nboot
      
      seed = assert_integerish(seed, null.ok = T, coerce = T)
      self$seed = seed
      
      assert_choice(target_type, c("fixed", "random"))
      self$target_type = target_type
      
      assert_flag(standardize)
      self$standardize = standardize
      assert_flag(intercept)
      self$intercept = intercept
      assert_flag(center_outcome)
      self$center_outcome = center_outcome
      
      if(intercept){
        #TODO: Stub
      }
      
      if(center_outcome){
        self$y = self$y - mean(self$y)
      }
      
      if(is.null(attr(self$x, "uposi:center"))) {
        # center columns of x
        means <- colMeans(self$x)
        self$x <- sweep(self$x, 2, means, "-")
        attr(self$x, "uposi:center") <- means
      }
      
      if(standardize) {
        # make standardization idempotent
        if(is.null(attr(self$x, "uposi:scale"))) {
          # scale columns of x
          l2 <- sqrt(colMeans(self$x^2))
          self$x <- sweep(self$x, 2, l2, "/")
          attr(self$x, "uposi:scale") <- l2
        }
      }
      
      if(self$stage == 1) {
        assert_r6(uposi_stage2, "Uposi", public=c("x","y","M"))
        assert_integerish(a2, lower=0, upper=1, any.missing = F, len=self$n)
        if(test_matrix(original_blip, nrows=self$n, ncols=1))
          original_blip = as.numeric(original_blip)
        assert_numeric(original_blip, any.missing=F, len=self$n)
      }
      
      self$uposi_stage2 = uposi_stage2
      self$a2 = a2
      self$original_blip = original_blip
      
      invisible(self)
      
    },
    
    do = function(Nboot = NULL) {
      self$
        bootstrap(Nboot)$
        extract_quantile()$
        confint()
    },
    
    bootstrap = function(Nboot = NULL) {
      if(!is.null(self$seed))
        set.seed(self$seed)
      
      if(!is.null(Nboot)) {
        Nboot = assert_count(Nboot, positive=T, coerce=T)
        self$Nboot = Nboot
      }
      
      self$boots = if(self$stage == 1) {
        future_replicate(self$Nboot, self$boot_iter_stage1())
      } else {
        future_replicate(self$Nboot, self$boot_iter_stage2())
      }
      
      invisible(self)
    },
    
    boot_iter_stage2 = function() {
      e <- private$gen_mult_norm()
      ret <- switch(
        self$target_type,
        fixed = private$xy_iter_stage2(e),
        random = c(private$xy_iter_stage2(e), private$xx_iter_stage2(e)),
        stop("Invalid target_type!")
      )
      return(ret)
    },
    
    boot_iter_stage1 = function() {
      e <- private$gen_mult_exp()
      
      new_coef <- qr.coef(
        qr(private$row_mult(
          self$uposi_stage2$x[,self$uposi_stage2$M, drop=F], sqrt(e)
        )),
        sqrt(e) * self$uposi_stage2$y
      )
      
      new_blip <- private$blip(new_coef)
      ret <- switch(
        self$target_type,
        fixed = private$xy_iter_stage1(e, new_blip),
        random = c(private$xy_iter_stage1(e, new_blip), 
                   private$xx_iter_stage1(e)),
        stop("Invalid target_type!")
      )
      # ret <- if(self$target_type == "fixed") {
      #   private$xy_iter_stage1(e, new_blip)
      # } else {
      #   c(private$xy_iter_stage1(e, new_blip), private$xx_iter_stage1(e))
      # }
      return(ret)
    },
    
    extract_quantile = function() {
      switch(
        self$target_type,
        fixed = assert(
          check_matrix(self$boots, nrows=1, ncols=self$Nboot, mode="numeric",
                       all.missing=F),
          check_numeric(self$boots, len=self$Nboot, all.missing = F)
        ),
        random = assert_matrix(self$boots, nrows=2, ncols=self$Nboot, 
                               mode="numeric", all.missing=F),
        stop("Invalid target_type in extract_quantile!")
      )
      # if(self$target_type == "fixed"){
      #   assert(
      #     check_matrix(self$boots, nrows=1, ncols=self$Nboot, mode="numeric",
      #                  all.missing=F),
      #     check_numeric(self$boots, len=self$Nboot, all.missing = F)
      #   )
      # } else {
      #   assert_matrix(self$boots, nrows=2, ncols=self$Nboot, mode="numeric",
      #                 all.missing=F)
      # }
      
      if(self$target_type == "random") {
        coef_l1 <- sum(abs(self$coef))
        RHS <- self$boots[1,] + self$boots[2,] * coef_l1
      } else {
        RHS <- self$boots
      }
      self$C_alpha = unname(quantile(RHS, probs=1-self$alpha))
      
      self$test_stat_H0 <- max(abs(crossprod(self$x, self$y)))/self$n
      if(self$target_type == "random") RHS <- self$boots[1,]
      self$pval_H0 <- mean(RHS >= self$test_stat_H0)
      self$C_alpha_H0 <- unname(quantile(RHS, probs=1-self$alpha))
      
      invisible(self)
      
    },
    
    confint = function() {
      assert_numeric(self$C_alpha, len=1, all.missing=F, finite=T)
      
      if(self$standardize){
        unscale <- attr(self$x, "uposi:scale")[self$M]
      } else {
        unscale <- 1
      }
      if(self$intercept)
        stop("Testing the intercept is unsupported at this time.")
      self$Sigmainv <- solve(self$Sigma)
      ci_half_len <- unscale * rowSums(abs(self$Sigmainv)) * self$C_alpha
      ci <- cbind(self$coef - ci_half_len, self$coef + ci_half_len)
      colnames(ci) <- c("Lower", "Upper")
      rownames(ci) <- colnames(self$x[,self$M])
      self$ci <- ci
      
      ci_half_len <- unscale * diag(self$Sigmainv) * self$C_alpha
      ci_2 <- cbind(self$coef - ci_half_len, self$coef + ci_half_len)
      colnames(ci_2) <- c("Lower", "Upper")
      rownames(ci_2) <- colnames(self$x[,self$M])
      self$ci_2 <- ci_2
      invisible(self)
    },
    
    pred_int = function(X) {
      stopifnot(ncol(X) == nrow(self$Sigmainv))
      xSinv_l1 <- rowSums(abs(X %*% self$Sigmainv))
      pred_half_len <- self$C_alpha * xSinv_l1
      pred <- as.numeric(X %*% self$coef)
      self$pred_i <- sapply(seq.int(self$n), function(i){
        pred[i] + c(-1,1)*pred_half_len[i]
      })
      return(self$pred_i)
    }
  ),
  
  private = list(
    # extract_quantile_fixed = function() {
    #   assert(
    #     check_matrix(self$boots, nrows=1, ncols=self$Nboot, mode="numeric",
    #                  all.missing=F),
    #     check_numeric(self$boots, len=self$Nboot, all.missing = F)
    #   )
    #   
    #   self$C_alpha = unname(quantile(self$boots, probs=1-self$alpha))
    #   
    #   test_stat_H0 <- max(abs(crossprod(self$x, self$y)))
    #   
    #   RHS <- self$boots
    #   self$pval_H0 <- mean(RHS >= test_stat_H0)
    #   self$C_alpha_H0 <- unname(quantile(RHS, probs=1-self$alpha))
    #   
    #   invisible(self)
    # },
    # 
    # extract_quantile_random = function() {
    #   
    #   assert_matrix(self$boots, nrows=2, ncols=self$Nboot, mode="numeric",
    #                 all.missing=F)
    #   
    #   test_stat_H0 <- max(abs(self$Sigma %*% self$coef))
    #   
    #   coef_l1 <- sum(abs(self$coef))
    #   
    #   RHS <- self$boots[1,] + self$boots[2,] * coef_l1
    #   self$C_alpha <- unname(quantile(RHS, probs=1-self$alpha))
    #   
    #   RHS <- self$boots[1,]
    #   self$pval_H0 <- mean(RHS >= test_stat_H0)
    #   self$C_alpha_H0 <- unname(quantile(RHS, probs=1-self$alpha))
    #   
    #   invisible(self)
    # },
    
    gen_mult = function() {
      if(stage == 2) {
        e <- private$gen_mult_norm()
      } else {
        e <- private$gen_mult_exp()
      }
      return(e)
    },
    
    gen_mult_norm = function() {
      e <- rnorm(self$n)
      return(e - mean(e))
    },
    
    gen_mult_exp = function() {
      e <- rexp(self$n, 1)
      return(e/mean(e))
    },
    
    xy_iter_stage2 = function(e) {
      return(
        max(abs(crossprod(self$x, self$y * e)/self$n))
      )
    },
    
    xx_iter_stage2 = function(e) {
      tmp <- abs(crossprod(self$x, private$row_mult(self$x, e))/self$n)
      return(
        max(tmp)
      )
    },
    
    xy_iter_stage1 = function(e, blip) {
      assert(
        check_numeric(blip, len=self$n),
        check_matrix(blip, mode="numeric", nrows=self$n, ncols=1)
      )
      return(
        max(abs(
          crossprod(
            self$x,
            e * (self$y + blip - self$original_blip) - (self$y)
          )/self$n
        ))
      )
    },
    
    xx_iter_stage1 = function(e) {
      return(max(abs(crossprod(self$x, private$row_mult(self$x, e-1))/self$n)))
    },
    
    blip = function(coefs) {
      delta <- as.numeric(self$uposi_stage2$x[,self$uposi_stage2$M, drop=F] %*% coefs)
      return(pmax(delta, 0) - self$a2 * delta)
    },
    
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
      return(x*y[row(x)])
    },
    
    #' Helper: Column multiplication
    #'
    #' Multiply the rows of a matrix x by a vector or column matrix y.
    #' This is equivalent to (but faster than) x %*% diag(y).
    #' @param x a matrix
    #' @param y a vector
    #'
    #' @return x %*% diag(y), but faster
    #'
    col_mult = function(x, y) {
      return(x*y[col(x)])
    }
    
  )
)

UposiFixed <- R6::R6Class(
  "UposiFixed",
  inherit = Uposi,
  public = list(
    
    
    boot_iter = function() {
      e <- super$gen_mult()
      return(super$xy_iter(e))
    },
    
    extract_quantile = function() {
      stopifnot(
        !is.null(self$boots),
        is.numeric(self$boots),
        length(self$boots) == self$Nboot
      )
      self$Sigma <- crossprod(self$x[,self$M])/self$n
      
      self$C_alpha = unname(quantile(self$boots, probs=1-self$alpha))
      
      test_stat_H0 <- max(abs(self$Sigma %*% self$coef))
      intercept <- sum(self$x[,self$M[1]]*self$y) / sum(self$x[,self$M[1]]^2)
      test_stat_H0i <- max(abs(
        self$Sigma %*% c(self$coef[1,] - intercept, self$coef[-1,])
      ))
      
      RHS <- self$boots
      self$pval_H0 <- mean(RHS >= test_stat_H0)
      self$C_alpha_H0 <- unname(quantile(RHS, probs=1-self$alpha))
      
      self$pval_H0i <- mean(RHS >= test_stat_H0i)
      self$C_alpha_H0i <- unname(quantile(RHS, probs=1-self$alpha))
      
      invisible(self)
    }
  )
)


UposiRandom <- R6::R6Class(
  "UposiRandom",
  inherit = Uposi,
  public = list(
    
    boot_iter = function() {
      e <- super$gen_mult()
      return(c(super$xy_iter(e), super$xx_iter(e)))
    },
    
    extract_quantile = function() {
      stopifnot(
        !is.null(self$boots),
        is.matrix(self$boots),
        ncol(self$boots) == self$Nboot,
        nrow(self$boots) == 2,
        !is.null(self$coef),
        !is.null(self$alpha)
      )
      self$Sigma <- crossprod(self$x[,self$M])/self$n
      test_stat_H0 <- max(abs(self$Sigma %*% self$coef))
      intercept <- sum(self$x[,self$M[1]]*self$y) / sum(self$x[,self$M[1]]^2)
      test_stat_H0i <- max(abs(
        self$Sigma %*% c(self$coef[1,] - intercept, self$coef[-1,])
      ))
      
      coef_l1 <- sum(abs(self$coef))
      # if(inherits(self, "UposiStage1"))
      #   browser()
      RHS <- self$boots[1,] + self$boots[2,] * coef_l1
      self$C_alpha <- unname(quantile(RHS, probs=1-self$alpha))
      
      RHS <- self$boots[1,]
      self$pval_H0 <- mean(RHS >= test_stat_H0)
      self$C_alpha_H0 <- unname(quantile(RHS, probs=1-self$alpha))
      
      RHS <- self$boots[1,] + self$boots[2,] * as.numeric(abs(self$coef[1,]))
      self$pval_H0i <- mean(RHS >= test_stat_H0i)
      self$C_alpha_H0i <- unname(quantile(RHS, probs=1-self$alpha))
      
      invisible(self)
    }
    
  )
)

UposiStage1 <- R6::R6Class(
  "UposiStage1",
  inherit = UposiRandom,
  public = list(
    
    uposi_stage2 = NULL,
    a2 = NULL,
    original_blip = NULL,
    
    initialize = function(
      uposi_stage2, a2, original_blip,
      ...
    ) {
      if(is.matrix(original_blip))
        original_blip <- as.numeric(original_blip)
      
      stopifnot(
        inherits(uposi_stage2, "UposiRandom"),
        is.numeric(a2), length(a2) == self$n,
        is.numeric(original_blip), length(original_blip) == self$n
      )
      
      super$initialize(...)
      self$uposi_stage2 = uposi_stage2
      self$a2 = a2
      self$original_blip = original_blip
      invisible(self)
    },
    
    gen_mult = function() { # overrode the multiplier to facilitate lm's
      stopifnot(!is.null(self$n))
      e <- rexp(self$n) # mean 1, variance 1
      e <- e / mean(e)
      return(e)
    },
    
    boot_iter = function() {
      e <- self$gen_mult() # self because we overrode the multiplier
      
      new_coef <- qr.coef(
        qr(super$row_mult(
          self$uposi_stage2$x[,self$uposi_stage2$M, drop=F], sqrt(e)
        )),
        sqrt(e) * self$uposi_stage2$y
      )
      
      new_blip <- private$blip(new_coef)
      return(c(private$xy_iter(e, new_blip), private$xx_iter(e)))
    }
    
  ),
  
  private = list(
    xy_iter = function(e, blip) {
      stopifnot(
        !is.null(self$x),
        !is.null(self$y),
        !is.null(self$n),
        is.numeric(e), length(e) == self$n,
        (
          is.numeric(blip) && length(blip) == self$n
        ) || (
          is.matrix(blip) && nrow(blip) == self$n && ncol(blip) == 1
        )
      )
      return(
        max(abs(
          crossprod(
            self$x,
            e * (self$y + blip - self$original_blip) - (self$y)
          )/self$n
        ))
      )
    },
    
    xx_iter = function(e) {
      stopifnot(
        !is.null(self$x),
        !is.null(self$y),
        !is.null(self$n),
        is.numeric(e), length(e) == self$n
      )
      return(max(abs(crossprod(self$x, super$row_mult(self$x, e-1))/self$n)))
    },
    
    blip = function(coefs) {
      stopifnot(
        is.numeric(coefs), length(coefs) == ncol(self$uposi_stage2$x[,self$uposi_stage2$M, drop=F]),
        is.numeric(self$a2), length(self$a2) == self$n
      )
      
      delta <- as.numeric(self$uposi_stage2$x[,self$uposi_stage2$M, drop=F] %*% coefs)
      return(pmax(delta, 0) - self$a2 * delta)
    }
  )
)
