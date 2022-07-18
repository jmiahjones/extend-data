library(R6)

Uposi <- R6::R6Class(
  "Uposi",
  public = list(
    
    x = NULL,
    Sigma = NULL,
    Sigmainv = NULL,
    y = NULL,
    n = NULL,
    p = NULL,
    alpha = NULL,
    Nboot = NULL,
    seed = NULL,
    
    M = NULL,
    coef = NULL,
    boots = NULL,
    C_alpha = NULL,
    C_alpha_H0 = NULL,
    C_alpha_H0i = NULL,
    ci = NULL,
    ci_2 = NULL,
    pred_i = NULL,
    pval = NULL,
    pval_H0 = NULL,
    pval_H0i = NULL,
    
    initialize = function(
      x, y, M = NULL, coef=NULL, alpha = 0.05, Nboot = 1000, seed = NULL
    ) {
      # TODO: How to handle intercept? Ensure X doesn't have the constant
      # column? Translate M to intercept numbering?
      self$x = as.matrix(x)
      self$y = private$reshape_colmat(y)
      self$n = nrow(x)
      self$p = ncol(x)
      
      stopifnot(
        nrow(y) == self$n, ncol(y) == 1,
        length(alpha) == 1, is.numeric(alpha), 0 < alpha, alpha < 1,
        length(Nboot) == 1, is.numeric(Nboot),
        is.null(seed) || (length(seed) == 1 && is.numeric(seed)),
        is.null(M) || (is.numeric(M)),
        is.null(coef) || is.matrix(coef) || is.numeric(coef)
      )
      
      if(!is.null(M)){
        M <- as.integer(M)
        stopifnot(
          min(M) >= 1, max(M) <= self$p
        )
        self$M = M
        
        if(is.null(coef)) {
          self$coef = qr.coef(qr(self$x[,M]), self$y)
        }
        coef <- private$reshape_colmat(coef)
        stopifnot(nrow(coef) == length(M))
        self$coef <- coef
      }
      
      Nboot = as.integer(Nboot)
      if(!is.null(seed))
        seed = as.integer(seed)
      
      self$alpha = alpha
      self$Nboot = Nboot
      self$seed = seed
      
    },
    
    do = function(Nboot = NULL) {
      self$
        bootstrap(Nboot)$
        extract_quantile()$
        confint()
    },
    
    bootstrap = function(Nboot) {
      if(!is.null(self$seed))
        set.seed(self$seed)
      
      stopifnot(is.null(Nboot) || (length(Nboot) == 1 && is.numeric(Nboot)))
      if(!is.null(Nboot))
        self$Nboot = as.integer(Nboot)
      
      self$boots = replicate(self$Nboot, self$boot_iter())
      invisible(self)
    },
    
    extract_quantile = function() {
      stop("Not implemented! Call the sub-classes.")
    },
    
    confint = function() {
      stopifnot(
        !is.null(self$x),
        !is.null(self$n),
        !is.null(self$M),
        !is.null(self$C_alpha),
        !is.null(self$coef)
      )
      # self$Sigma <- crossprod(self$x[,self$M])/self$n
      self$Sigmainv <- solve(self$Sigma)
      ci_half_len <- rowSums(abs(self$Sigmainv)) * self$C_alpha
      ci <- cbind(self$coef - ci_half_len, self$coef + ci_half_len)
      colnames(ci) <- c("Lower", "Upper")
      rownames(ci) <- colnames(self$x[,self$M])
      self$ci <- ci
      
      ci_half_len <- diag(self$Sigmainv) * self$C_alpha
      ci_2 <- cbind(self$coef - ci_half_len, self$coef + ci_half_len)
      colnames(ci_2) <- c("Lower", "Upper")
      rownames(ci_2) <- colnames(self$x[,self$M])
      self$ci_2 <- ci_2
      invisible(self)
    },
    
    pred_int = function(X) {
      stopifnot(ncol(X) == nrow(self$Sigmainv))
      # xSinvx <- rowSums((X %*% self$Sigmainv) * X)
      # x_sel_max <- apply(abs(X), 1, max)
      # stopifnot(is.numeric(x_sel_max))
      # pred_half_len <- self$C_alpha * (xSinvx / x_sel_max)
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
    reshape_colmat = function(x) {
      stopifnot(is.matrix(x) || is.numeric(x))
      if(is.matrix(x))
        stopifnot(ncol(x) == 1)
      
      if(is.numeric(x)) {
        x <- matrix(x, ncol = 1)
      }
      
      return(x)
    },
    
    gen_mult = function() {
      stopifnot(
        !is.null(self$n)
      )
      e <- rnorm(self$n)
      return(e - mean(e))
    },
    
    xy_iter = function(e) {
      stopifnot(
        !is.null(self$x),
        !is.null(self$y),
        !is.null(self$n),
        is.numeric(e), length(e) == self$n
      )
      return(
        max(abs(crossprod(self$x, self$y * e)/self$n))
      )
    },
    
    xx_iter = function(e) {
      stopifnot(
        !is.null(self$x),
        !is.null(self$y),
        !is.null(self$n),
        is.numeric(e), length(e) == self$n
      )
      
      tmp <- abs(crossprod(self$x, private$row_mult(self$x, e))/self$n)
      return(
        max(tmp)
      )
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


UposiRandomWeight <- R6::R6Class(
  "UposiRandomWeight",
  inherit = UposiRandom,
  public = list(
    
    weights = NULL,
    
    initialize = function(..., weights) {
      super$initialize(...)
      stopifnot(length(weights) == self$p)
      self$weights <- weights
    },
    
    boot_iter = function() {
      e <- super$gen_mult()
      return(c(private$xy_iter(e), private$xx_iter(e)))
    }
  ),
  
  private = list(
    
    xy_iter = function(e) {
      stopifnot(
        !is.null(self$x),
        !is.null(self$y),
        !is.null(self$n),
        is.numeric(e), length(e) == self$n
      )
      return(
        max(self$weights * abs(crossprod(self$x, self$y * e)/self$n))
      )
    },
    
    xx_iter = function(e) {
      stopifnot(
        !is.null(self$x),
        !is.null(self$y),
        !is.null(self$n),
        is.numeric(e), length(e) == self$n
      )
      
      tmp <- abs(crossprod(self$x, super$row_mult(self$x, e))/self$n)
      tmp <- super$row_mult(super$col_mult(tmp, self$weights), self$weights)
      return(
        max(tmp)
      )
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
            e * (self$y + blip) - (self$y + self$original_blip)
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

## Stg 1 Weight
UposiStage1Weight <- R6::R6Class(
  "UposiStage1Weight",
  inherit = UposiStage1,
  public = list(
    
    weights = NULL,
    
    initialize = function(..., weights) {
      super$initialize(...)
      stopifnot(length(weights) == self$p)
      self$weights <- weights
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
          self$weights * crossprod(
            self$x, 
            e * (self$y + blip) - (self$y + self$original_blip)
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
      
      tmp <- abs(crossprod(self$x, super$row_mult(self$x, e-1))/self$n)
      tmp <- super$row_mult(super$col_mult(tmp, self$weights), self$weights)
      return(max(tmp))
    }
  )
)