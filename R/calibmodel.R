#' @export
calibmodel <- function(X, y, engine=stats::lm, lambda=NULL, seed=123, ...) {
  set.seed(seed)  
  n_train <- floor(0.5 * nrow(X))
  y_mean <- mean(y)  
  X_mean <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_scaled <- scale(X, center=X_mean, scale=X_sd)
  y <- y - y_mean
  p <- ncol(X)
  if (!is.null(lambda)) {
    if (length(lambda) > 1){
      lambda <- select_ridge_lambda(X_scaled, y, lambda)$lambda.min
    }
    sqrt_lambda <- sqrt(lambda)    
  }
  # Obtain calibration residuals
  # Split data into calibration and training sets
  train_idx <- sample(nrow(X), size=n_train)
  X_train <- X_scaled[train_idx,]
  y_train <- y[train_idx]
  X_calib <- X_scaled[-train_idx,]
  y_calib <- y[-train_idx]                
  # Fit model on training data with augmented matrices
  if (!is.null(lambda)) {
  diag_sqrt_lambda <- diag(sqrt_lambda, p, p)
  Z_train <- rbind(X_train, diag_sqrt_lambda)
  Z_calib <- rbind(X_calib, diag_sqrt_lambda)
  } else {
    Z_train <- X_train
    Z_calib <- X_calib
  }
  # Compare functions using identical()
  if (!is.null(lambda)) {
    rep_0_p <- rep(0, p)
    y_train_aug <- c(y_train, rep_0_p)
    y_calib_aug <- c(y_calib, rep_0_p)      
  } else {
    y_train_aug <- y_train
    y_calib_aug <- y_calib
  }  
  df_train <- data.frame(y=y_train_aug, as.data.frame(Z_train))
  df_calib <- data.frame(y=y_calib_aug, as.data.frame(Z_calib))
  model <- engine(y ~ . -1, data=df_train, ...)        
  # Get calibration residuals
  calib_pred <- predict(model, newdata=as.data.frame(X_calib))
  calib_resid <- y_calib - calib_pred        
  # Store calibration residuals
  model <- engine(y ~ . -1, data=df_calib, ...)
  model$residuals <- drop(calib_resid)
  model$y_mean <- y_mean
  model$X_mean <- X_mean
  model$X_sd <- X_sd
  model$model <- model
  class(model) <- c("calibmodel", class(model$model))
  return(model)
}

#' @export
predict.calibmodel <- function(object, newdata, 
                             method=c("none", "gaussian", "surrogate", "bootstrap", "tsbootstrap"), 
                             level=95, nsim=100, seed=123, ...) {
  stopifnot(inherits(object, "calibmodel"))
  set.seed(seed)
  method <- match.arg(method)
    
    # Convert newdata to matrix if it isn't already
    if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
    
    # Check dimensions
    if (length(object$X_mean) != ncol(newdata)) {
      stop("Number of variables in newdata (", ncol(newdata), 
           ") must match the training data (", length(object$X_mean), ")")
    }
    
    # Scale new data using training scaling parameters
    newdata_scaled <- scale(newdata, 
                          center = object$X_mean,
                          scale = object$X_sd)
    
    # Convert to data frame and add column names if missing
    newdata_df <- as.data.frame(newdata_scaled)
    if (is.null(colnames(newdata_df))) {
      colnames(newdata_df) <- paste0("V", seq_len(ncol(newdata_df)))
    }
    
    # Use predict on the underlying model
    pred <- try(predict(object$model, newdata = newdata_df, ...), silent=TRUE)
    if (inherits(pred, "try-error")) {
      pred <- predict(object$model, as.matrix(newdata_df), ...)
    }
    # Get point predictions
    if (method == "none") {  
      return(drop(pred) + object$y_mean)
    } else {
      # Generate simulations
      sims <- simulate.calibmodel(object, newdata = newdata, nsim = nsim, 
                                  method = method, seed = seed)            
      # Calculate prediction intervals using empirical quantiles
      alpha <- (1 - level/100) / 2
      lower <- apply(sims, 1, quantile, probs = alpha)
      upper <- apply(sims, 1, quantile, probs = 1 - alpha)            
      # Return as matrix with columns fit, lwr, upr
      return(cbind(fit = drop(pred) + object$y_mean, 
                  lwr = lower,
                  upr = upper))
    }
  }


#' @export
simulate.calibmodel <- function(object, newdata, nsim = 100, 
                              method = c("gaussian", "surrogate", "bootstrap", "tsbootstrap"),
                              seed = 123, ...) {
  
  set.seed(seed)  
  method <- match.arg(method)
  
  # Scale newdata
  newdata_scaled <- scale(as.matrix(newdata), 
                         center = object$X_mean,
                         scale = object$X_sd)
  
  # Try prediction with data frame first
  pred <- try(drop(predict(object$model, newdata = as.data.frame(newdata_scaled))), silent = TRUE) + object$y_mean
  
  # If that fails, try with matrix
  if (inherits(pred, "try-error")) {
    pred <- drop(predict(object$model, newdata = newdata_scaled)) + object$y_mean
  }
  
  if (method == "gaussian") {    
    # Get the residual standard error from the model
    sigma <- sqrt(sum(object$residuals^2) / object$df.residual)
    # Generate random normal errors
    errors <- matrix(rnorm(length(pred) * nsim, 
                           mean = 0, 
                           sd = sigma), 
                     nrow = length(pred), 
                     ncol = nsim)
  }
  
  if (method == "surrogate") {
    errors <- matrix(0, nrow = length(pred), ncol = nsim)
    for (i in 1:nsim) {
      set.seed(seed + i*10)
      errors[,i] <- as.numeric(tseries::surrogate(sample(object$residuals, size=length(pred), replace=TRUE)))
    }
  }
  
  if (method == "bootstrap") {
    errors <- matrix(0, nrow = length(pred), ncol = nsim)
    for (i in 1:nsim) {
      set.seed(seed + i*100)
      errors[,i] <- sample(object$residuals, size=length(pred), replace=TRUE)
    }
  }
  
  if (method == "tsbootstrap") {
    errors <- matrix(0, nrow = length(pred), ncol = nsim)
    for (i in 1:nsim) {
      set.seed(seed + i*1000)     
      errors[,i] <- as.numeric(tseries::tsbootstrap(sample(object$residuals, size=length(pred), replace=TRUE), type="block"))
    }
  }
  
  # Add errors to fitted values to create simulations
  #misc::debug_print(errors)
  #misc::debug_print(pred)
  #misc::debug_print(dim(errors))
  #misc::debug_print(dim(pred))
  #misc::debug_print(length(pred))
  result <- errors + drop(pred)
  # Convert to data.frame to match simulate.lm output format
  result <- as.data.frame(result)
  names(result) <- paste0("sim_", 1:nsim)
  return(result)
}

#' @export
calibmodel.formula <- function(formula, data, engine=stats::lm, 
                               lambda=0.1, seed=123, ...) {
  # Extract X matrix and y vector from formula and data
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)[,-1, drop=FALSE]  # Remove intercept column  
  # Call the original caliblm function
  result <- calibmodel(X, y, engine=engine, 
                             lambda=lambda, seed=seed, ...)  
  # Add formula-related attributes
  result$call <- match.call()
  result$terms <- terms(formula, data=data)
  result$model <- mf  
  result$model$model <- result
  return(result)
}


#' @export
select_ridge_lambda <- function(X, y, lambda_seq = 10**seq(-10, 10, length.out=100), seed = 123) {
  # Ensure X is a matrix and y is a vector
  X <- as.matrix(X)
  y <- as.vector(y)
  
  # Center and scale X and y
  y_mean <- mean(y)
  X_mean <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_scaled <- scale(X, center = X_mean, scale = X_sd)
  y <- y - y_mean
  n <- nrow(X)
  p <- ncol(X)
  
  # Calculate SVD of X once
  svd_X <- svd(X_scaled)
  d <- svd_X$d
  U <- svd_X$u
  V <- svd_X$v
  
  # Calculate GCV for each lambda
  gcv_scores <- sapply(lambda_seq, function(lambda) {
    # Calculate effective degrees of freedom
    df <- sum(d^2 / (d^2 + lambda))
    
    # Calculate ridge coefficients
    ridge_coef <- V %*% (diag(d/(d^2 + lambda)) %*% (t(U) %*% y))
    
    # Calculate fitted values
    fitted <- X_scaled %*% ridge_coef
    
    # Calculate residual sum of squares
    rss <- sum((y - fitted)^2)
    
    # Calculate GCV
    gcv <- (rss/n) / (1 - df/n)^2
    return(gcv)
  })
  
  # Find lambda with minimum GCV
  best_idx <- which.min(gcv_scores)
  best_lambda <- lambda_seq[best_idx]
  
  # Return results
  result <- list(
    lambda = lambda_seq,
    GCV = gcv_scores,
    lambda.min = best_lambda
  )
  class(result) <- "ridge_lambda"
  return(result)
}