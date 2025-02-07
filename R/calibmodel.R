#' @export
calibmodel <- function(X, y, engine=stats::lm, lambda=NULL, positive_response=FALSE, seed=123, ...) {
  set.seed(seed) 
  if (!is.null(lambda)) {
    if (length(lambda) > 1){
      lambda <- select_ridge_lambda(X, y, lambda)$lambda.min
    }
    sqrt_lambda <- sqrt(lambda)    
  }   
  if (positive_response && any(y < 0)) {
    stop("Positive response is set to TRUE but some values of y are negative")
  }
  y_mean <- 0 # keep this that way 
  X_mean <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_scaled <- scale(X, center=X_mean, scale=X_sd)
  if (positive_response == FALSE) {
    y_mean <- mean(y)
    response <- y - y_mean
  } else {
    response <- y
  }
  p <- ncol(X)
  # Obtain calibration residuals
  # Split data into calibration and training sets
  train_idx <- as.integer(create_data_partition(response, p=0.5)[[1]])
  X_train <- X_scaled[train_idx, ]
  y_train <- response[train_idx]
  X_calib <- X_scaled[-train_idx, ]
  y_calib <- response[-train_idx]                
  # Fit model on training data with augmented matrices
  if (!is.null(lambda) || positive_response == FALSE) {
    diag_sqrt_lambda <- diag(sqrt_lambda, p, p)
    Z_train <- rbind(X_train, diag_sqrt_lambda)
    rep_0_p <- rep(0, p)
    y_train_aug <- c(y_train, rep_0_p)
    df_train <- data.frame(y=y_train_aug, as.data.frame(Z_train))
    df_calib <- data.frame(y=y_calib, as.data.frame(X_calib))
    object <- engine(y ~ . -1, data=df_train, ...)     
    # Get calibration residuals
    calib_pred <- try(predict(object, as.data.frame(X_calib)), silent=TRUE)
    if (inherits(calib_pred, "try-error")) {
      calib_pred <- predict(object, as.matrix(X_calib))
    }
    calib_resid <- y_calib - calib_pred           
  } else { # positive response and no lambda
    df_train <- data.frame(y=y_train, as.data.frame(X_train))
    df_calib <- data.frame(y=y_calib, as.data.frame(X_calib))
    object <- engine(y ~ ., data=df_train, ...)
    calib_pred <- try(predict(object, as.data.frame(X_calib)), silent=TRUE)
    if (inherits(calib_pred, "try-error")) {
      calib_pred <- predict(object, as.matrix(X_calib))
    }
    calib_resid <- y_calib/calib_pred           
  }    
  # Store calibration residuals
  fit_obj <- vector("list", 6)
  fit_obj$positive_response <- positive_response
  fit_obj$residuals <- drop(calib_resid)
  fit_obj$y_mean <- y_mean
  fit_obj$X_mean <- X_mean
  fit_obj$X_sd <- X_sd
  fit_obj$model <- object
  class(fit_obj) <- c("calibmodel", class(object))
  return(fit_obj)
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
    pred <- try(predict(object$model, newdata = newdata_df, ...), silent=TRUE) + object$y_mean
    if (inherits(pred, "try-error")) {
      pred <- predict(object$model, as.matrix(newdata_df), ...) + object$y_mean
    }
    # Get point predictions
    if (method == "none") {  
      return(drop(pred))
    } else {
      # Generate simulations
      sims <- simulate.calibmodel(object, newdata = newdata, nsim = nsim, 
                                  method = method, seed = seed)            
      # Calculate prediction intervals using empirical quantiles
      alpha <- (1 - level/100) / 2
      pred <- apply(sims, 1, median)
      lower <- apply(sims, 1, quantile, probs = alpha)
      upper <- apply(sims, 1, quantile, probs = 1 - alpha)            
      # Return as matrix with columns fit, lwr, upr
      return(cbind(fit = drop(pred), 
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
  pred <- try(drop(predict(object$model, as.data.frame(newdata_scaled))), silent = TRUE) + object$y_mean
  if (inherits(pred, "try-error")) {
    pred <- drop(predict(object$model, as.matrix(newdata_scaled))) + object$y_mean
  }
  
  if (method == "gaussian") { 
    stopifnot(object$positive_response == FALSE)   
    # Get the residual standard error from the model
    sigma <- try(sqrt(sum(object$residuals^2) / object$model$df.residual), silent=TRUE)
    if (inherits(sigma, "try-error")) {
      sigma <- try(sqrt(sum(object$residuals^2) / (length(object$residuals) - ncol(newdata_scaled) - 1)), silent=TRUE)
    }
    # Generate random normal errors
    errors <- matrix(rnorm(length(pred) * nsim, 
                           mean = 0, 
                           sd = sigma), 
                     nrow = length(pred), 
                     ncol = nsim)
  }
  
  if (method == "surrogate") {
    #errors <- matrix(0, nrow = length(pred), ncol = nsim)
    #for (i in 1:nsim) {
    #  set.seed(seed + i*10)
    #  errors[, i] <- as.numeric(tseries::surrogate(sample(object$residuals, size=length(pred), replace=TRUE)))
    #}
    errors <- tseries::surrogate(sample(object$residuals, size=length(pred), replace=TRUE), ns=nsim)
  }
  
  if (method == "bootstrap") {
    errors <- matrix(sample(object$residuals, size=length(pred)*nsim, replace=TRUE), nrow=length(pred), ncol=nsim)    
  }
  
  if (method == "tsbootstrap") {
    #errors <- matrix(0, nrow = length(pred), ncol = nsim)
    #for (i in 1:nsim) {
    #  set.seed(seed + i*1000)     
    #  errors[, i] <- as.numeric(tseries::tsbootstrap(sample(object$residuals, size=length(pred), replace=TRUE), type="block"))
    #}
    errors <- tseries::tsbootstrap(sample(object$residuals, size=length(pred), replace=TRUE), type="block", nb=nsim)
  }
  # Add errors to fitted values to create simulations
  if (object$positive_response == FALSE) {
    result <- errors + drop(pred)
  } else {
    result <- errors * drop(pred)
  }
  # Convert to data.frame to match simulate.lm output format
  result <- as.data.frame(result)
  colnames(result) <- paste0("sim_", 1:nsim)
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
  
  # Calculate GCV for all lambda values at once
  div <- d^2 + rep(lambda_seq, rep(length(d), length(lambda_seq)))
  dim(div) <- c(length(d), length(lambda_seq))
  
  # Calculate coefficients for all lambda values
  a <- drop(d * (t(U) %*% y)) / div
  ridge_coef <- V %*% a
  
  # Calculate fitted values for all lambda values
  fitted <- X_scaled %*% ridge_coef
  
  # Calculate effective degrees of freedom for all lambda values
  df <- colSums(matrix(d^2/div, length(d)))
  
  # Calculate GCV scores
  resid <- y - fitted
  rss <- colSums(resid^2)
  gcv_scores <- (rss/n) / (1 - df/n)^2
  
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