#' @export
ridgemodel <- function(X, y, workhorse=stats::lm, lambda=10**seq(-10, 10, length.out=100), 
                       type_residuals=c("gaussian", "calibrated"), seed=123, ...) {
  set.seed(seed)  
  n_train <- floor(0.5 * nrow(X))
  y_mean <- mean(y)
  X_mean <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_scaled <- scale(X, center=X_mean, scale=X_sd)
  y <- y - y_mean
  p <- ncol(X)
  if (length(lambda) > 1){
    lambda <- select_ridge_lambda(X_scaled, y, lambda)$lambda.min
  }
  sqrt_lambda <- sqrt(lambda)
  Z <- rbind(X_scaled, diag(sqrt_lambda, p, p))
  y_aug <- c(y, rep(0, p))
  # Create data frame with augmented data
  df <- data.frame(y = y_aug, as.data.frame(Z))
  type_residuals <- match.arg(type_residuals)        
  
  if (type_residuals == "gaussian") {        
    # Original code for other workhorse functions 
    model <- try(workhorse(y ~ . -1, data=df, ...), silent=TRUE)
    if (inherits(model, "try-error")) {
      model <- workhorse(X = as.matrix(X), y = y, ...)
    }
    model$type_residuals <- "gaussian"
  }
  
  if (type_residuals == "calibrated") {
    # Obtain calibration residuals
    # Split data into calibration and training sets
    train_idx <- sample(nrow(X), size=n_train)
    X_train <- X_scaled[train_idx,]
    y_train <- y[train_idx]
    X_calib <- X_scaled[-train_idx,]
    y_calib <- y[-train_idx]                
    # Fit model on training data with augmented matrices
    Z_train <- rbind(X_train, diag(sqrt_lambda, p, p))
    Z_calib <- rbind(X_calib, diag(sqrt_lambda, p, p))
    y_train_aug <- c(y_train, rep(0, p))
    y_calib_aug <- c(y_calib, rep(0, p))        
    df_train <- data.frame(y=y_train_aug, as.data.frame(Z_train))
    df_calib <- data.frame(y=y_calib_aug, as.data.frame(Z_calib))
    model <- workhorse(y ~ . -1, data=df_train, ...)        
    # Get calibration residuals
    calib_pred <- predict(model, newdata=as.data.frame(X_calib))
    calib_resid <- y_calib - calib_pred        
    # Store calibration residuals
    model <- workhorse(y ~ . -1, data=df_calib, ...)
    model$residuals <- calib_resid
    model$type_residuals <- "calibrated"
  }
  
  model$y_mean <- y_mean
  model$X_mean <- X_mean
  model$X_sd <- X_sd
  class(model) <- c("ridgemodel", "lm")
  return(model)
}

#' @export
predict.ridgemodel <- function(object, newdata, 
                               method=c("gaussian", "surrogate", "bootstrap", "tsbootstrap"), 
                               interval=c("none", "prediction"), seed=123, ...) {
  stopifnot(inherits(object, "ridgemodel"))
  set.seed(seed)
  method <- match.arg(method)
  interval <- match.arg(interval)
  if (object$type_residuals == "gaussian") {
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
    # Make predictions
    pred <- stats::predict.lm(object, as.data.frame(newdata_scaled), interval=interval, ...) + object$y_mean
    return(drop(pred))
  }
  if (object$type_residuals == "calibrated") {        
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
    # Get point predictions
    pred <- stats::predict.lm(object, as.data.frame(newdata_scaled), interval="none") + object$y_mean        
    # Check if interval prediction is requested
    if (interval == "prediction") {
      # Get level from ... or use default 0.95
      level <- if (!is.null(list(...)$level)) list(...)$level else 0.95            
      # Generate simulations
      sims <- simulate.ridgemodel(object, newdata = newdata, nsim = 1000, 
                                  method=method, seed=seed)            
      # Calculate prediction intervals using empirical quantiles
      alpha <- (1 - level) / 2
      lower <- apply(sims, 1, quantile, probs = alpha)
      upper <- apply(sims, 1, quantile, probs = 1 - alpha)            
      # Return as matrix with columns fit, lwr, upr
      return(cbind(fit = drop(pred), 
                   lwr = lower,
                   upr = upper))
    }        
    return(drop(pred))
  }
}

#' @export
simulate.ridgemodel <- function(object, newdata, nsim = 100L, seed = 123, 
                                method=c("gaussian", "surrogate", "bootstrap", "tsbootstrap"), ...) {
  stopifnot(inherits(object, "ridgemodel"))
  set.seed(seed)  
  method <- match.arg(method)
  # Get predictions for new data
  fitted_values <- predict.ridgemodel(object, newdata = newdata, method="gaussian")
  
  if (method == "gaussian") {    
    # Get the residual standard error from the model
    sigma <- sqrt(sum(object$residuals^2) / object$df.residual)
    # Generate random normal errors
    errors <- matrix(rnorm(length(fitted_values) * nsim, 
                           mean = 0, 
                           sd = sigma), 
                     nrow = length(fitted_values), 
                     ncol = nsim)
  }
  
  if (method == "surrogate") {
    errors <- matrix(0, nrow = length(fitted_values), ncol = nsim)
    for (i in 1:nsim) {
      set.seed(seed + i*10)
      errors[,i] <- as.numeric(tseries::surrogate(sample(object$residuals, size=length(fitted_values), replace=TRUE)))
    }
  }
  
  if (method == "bootstrap") {
    errors <- matrix(0, nrow = length(fitted_values), ncol = nsim)
    for (i in 1:nsim) {
      set.seed(seed + i*100)
      errors[,i] <- sample(object$residuals, size=length(fitted_values), replace=TRUE)
    }
  }
  
  if (method == "tsbootstrap") {
    errors <- matrix(0, nrow = length(fitted_values), ncol = nsim)
    for (i in 1:nsim) {
      set.seed(seed + i*1000)     
      errors[,i] <- as.numeric(tseries::tsbootstrap(sample(object$residuals, size=length(fitted_values), replace=TRUE), type="block"))
    }
  }
  
  # Add errors to fitted values to create simulations
  result <- errors + fitted_values
  # Convert to data.frame to match simulate.lm output format
  result <- as.data.frame(result)
  names(result) <- paste0("sim_", 1:nsim)
  return(result)
}

#' @export
ridgemodel.formula <- function(formula, data, workhorse=stats::lm, 
                               lambda=0.1, seed=123, ...) {
  # Extract X matrix and y vector from formula and data
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)[,-1, drop=FALSE]  # Remove intercept column  
  # Call the original caliblm function
  result <- misc::ridgemodel(X, y, workhorse=workhorse, 
                             lambda=lambda, seed=seed, ...)  
  # Add formula-related attributes
  result$call <- match.call()
  result$terms <- terms(formula, data=data)
  result$model <- mf  
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