#' @export
rvfl <- function(X, y, lambda, n_hidden_features=5L, 
                            activation_function=c("relu", "sigmoid", "tanh"),
                            coeff=NULL, workhorse=NULL) {
  activation_function <- match.arg(activation_function)
  scaled_X <- scale(X)
  y_mean <- mean(y)
  centered_y <- y - y_mean    
  # Create augmented features if needed
  if (n_hidden_features > 0) {
    if (activation_function == "relu") {
      hidden_features <- pmax(scaled_X, 0)
    }
    else if (activation_function == "sigmoid") {
      hidden_features <- ifelse(scaled_X >= 0,
                                1 / (1 + exp(-scaled_X)),
                                exp(scaled_X) / (1 + exp(scaled_X)))
    }
    else if (activation_function == "tanh") {
      hidden_features <- tanh(scaled_X)            
    }
    scaled_X <- cbind(scaled_X, hidden_features)
  }    
  # If workhorse is lm or glm, use ridge regression formula
  if (is.null(workhorse) || inherits(workhorse, "lm") || inherits(workhorse, "glm")) {
    # Ridge regression solution: (X'X + λI)^(-1)X'y
    XtX <- crossprod(scaled_X)
    p <- ncol(XtX)
    ridge_term <- lambda * diag(p)
    Xty <- crossprod(scaled_X, centered_y)
    
    # Solve ridge regression
    coeff <- solve(XtX + ridge_term, Xty)
    
    # Create a list with necessary components
    result <- list(
      coefficients = coeff,
      x = scaled_X,
      y = centered_y,
      fitted.values = scaled_X %*% coeff,
      lambda = lambda,
      y_mean = y_mean,
      scaling = list(
        center = attr(scaled_X, "scaled:center"),
        scale = attr(scaled_X, "scaled:scale")
      )
    )
    class(result) <- c("nlgaussianridge", "lm")
    return(result)
  } else {
    # Use provided workhorse function
    df_covariates <- data.frame(X=scaled_X, y=centered_y)
    return(workhorse(X=scaled_X, y=centered_y))
  }
}

#' @export
predict.rvfl <- function(object, newdata, ...) {
  # Scale new data using training scaling parameters
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  scaled_newdata <- scale(newdata, 
                          center = object$scaling$center,
                          scale = object$scaling$scale)
  
  # Make predictions
  pred <- scaled_newdata %*% object$coefficients + object$y_mean
  return(drop(pred))
}

