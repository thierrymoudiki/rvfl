#'
#' Simple and efficient implementation using the kernel trick on concatenated features
#'
#' @param x Matrix of predictors (n x p)
#' @param y Response vector (n x 1) or matrix (n x m)
#' @param lambda Regularization parameter
#' @param nb_hidden Number of hidden units
#' @param activation Activation function
#' @param sigma Scale parameter for weights
#' @param seed Random seed
#' @param ... Additional arguments
#'
#' @return krvfl object
#' @export
krvfl <- function(x, y, lambda = 0.1, nb_hidden = 100L,
                  activation = c("relu", "sigmoid", "tanh", "linear"),
                  sigma = 1.0, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  # Input validation
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(y)) y <- as.matrix(y)
  if (nrow(x) != nrow(y)) stop("x and y must have same number of rows")
  activation <- match.arg(activation)
  n <- nrow(x)
  p <- ncol(x)
  # Center and scale predictors
  x_mean <- colMeans(x)
  x_sd <- apply(x, 2, sd)
  x_sd[x_sd == 0] <- 1
  x_scaled <- scale(x, center = x_mean, scale = x_sd)
  # Center response
  y_mean <- mean(y)
  y_centered <- y - y_mean
  # Generate random features
  W <- matrix(rnorm(nb_hidden * p, 0, sigma), nrow = p, ncol = nb_hidden)
  # Compute hidden features
  Z <- x_scaled %*% W
  Phi <- pmax(Z, 0)
  # Concatenate original and transformed features
  X_combined <- cbind(x_scaled, Phi)
  # Apply kernel trick: K = X_combined %*% t(X_combined)
  # But we'll use the explicit feature mapping for efficiency
  K <- tcrossprod(X_combined)
  # Solve: (K + lambda * I) alpha = y_centered
  K_reg <- K + lambda * diag(n)
  # Use Cholesky decomposition for stability
  R <- chol(K_reg)
  alpha <- backsolve(R, forwardsolve(t(R), y_centered))
  # Store the dual coefficients and feature matrix for prediction
  fitted_values <- K %*% alpha + y_mean
  result <- list(
    x = x,
    x_scaled = x_scaled,
    y = y,
    alpha = alpha,
    W = W,
    activation = activation,
    lambda = lambda,
    nb_hidden = nb_hidden,
    x_mean = x_mean,
    x_sd = x_sd,
    y_mean = y_mean,
    fitted.values = fitted_values,
    residuals = y - fitted_values,
    X_train_combined = X_combined,
    Phi_train = Phi,
    call = match.call()
  )
  class(result) <- "krvfl"
  return(result)
}

#' Predict method for krvfl objects
#'
#' @param object A krvfl object from krvfl()
#' @param newdata Matrix of new predictors (m x p)
#' @param ... Additional arguments (unused)
#'
#' @return Predicted values (m x 1) or (m x q) matrix
#' @export
predict.krvfl <- function(object, newdata = NULL, ...) {
  # If no new data, return fitted values
  if (is.null(newdata)) {
    return(as.vector(object$fitted.values))
  }
  # Input validation
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  if (ncol(newdata) != ncol(object$x)) {
    stop("newdata must have same number of columns as training data")
  }
  # Scale new data using training statistics
  newdata_scaled <- scale(newdata, 
                          center = object$x_mean, 
                          scale = object$x_sd)
  # Apply the same transformation as in training
  # Generate hidden features using stored weights
  Z_new <- newdata_scaled %*% object$W
  # Apply activation function
  Phi_new <- switch(object$activation,
                    "relu" = pmax(Z_new, 0),
                    "sigmoid" = 1 / (1 + exp(-Z_new)),
                    "tanh" = tanh(Z_new),
                    "linear" = Z_new
  )
  # Concatenate original and transformed features (same as training)
  X_new_combined <- cbind(newdata_scaled, Phi_new)
  # Compute kernel matrix between new data and training data
  # K_new = X_new_combined %*% t(X_training_combined)
  # We need to reconstruct the training combined features
  Z_train <- object$x_scaled %*% object$W
  Phi_train <- switch(object$activation,
                      "relu" = pmax(Z_train, 0),
                      "sigmoid" = 1 / (1 + exp(-Z_train)),
                      "tanh" = tanh(Z_train),
                      "linear" = Z_train
  )
  X_train_combined <- cbind(object$x_scaled, Phi_train)
  # Compute cross-kernel matrix
  K_new <- tcrossprod(X_new_combined, X_train_combined)
  # Make predictions: K_new %*% alpha + y_mean
  return(as.vector(K_new %*% object$alpha + object$y_mean))
}

#' Print method for krvfl objects
#'
#' @param x A krvfl object
#' @param ... Additional arguments (unused)
#' @export
print.krvfl <- function(x, ...) {
  cat("Kernel Random Vector Functional Link (KRVFL) Model\n")
  cat("=================================================\n")
  cat("Call: ", deparse(x$call), "\n\n")
  cat("Model details:\n")
  cat("- Training samples:", nrow(x$x), "\n")
  cat("- Features:", ncol(x$x), "\n")
  cat("- Hidden units:", x$nb_hidden, "\n")
  cat("- Activation:", x$activation, "\n")
  cat("- Regularization (lambda):", x$lambda, "\n")
  cat("- Training RMSE:", sqrt(mean(x$residuals^2)), "\n")
  cat("\nUse predict() to make predictions on new data.\n")
}

#' Summary method for krvfl objects
#'
#' @param object A krvfl object
#' @param ... Additional arguments (unused)
#' @export
summary.krvfl <- function(object, ...) {
  cat("Kernel Random Vector Functional Link (KRVFL) Summary\n")
  cat("===================================================\n")
  cat("Call: ", deparse(object$call), "\n\n")
  
  # Model specifications
  cat("Model Specifications:\n")
  cat("- Training samples:", nrow(object$x), "\n")
  cat("- Original features:", ncol(object$x), "\n")
  cat("- Hidden units:", object$nb_hidden, "\n")
  cat("- Total features (original + hidden):", ncol(object$x) + object$nb_hidden, "\n")
  cat("- Activation function:", object$activation, "\n")
  cat("- Regularization parameter:", object$lambda, "\n")
  
  # Performance metrics
  residuals <- object$residuals
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  
  cat("\nPerformance Metrics:\n")
  cat("- RMSE:", rmse, "\n")
  cat("- MAE:", mae, "\n")
  
  # Residual statistics
  cat("\nResidual Statistics:\n")
  print(summary(as.vector(residuals)))
  
  # Feature scaling info
  cat("\nFeature Scaling:\n")
  cat("- Features are centered and scaled\n")
  cat("- Response is centered (mean:", object$y_mean, ")\n")
}

# # Example usage and testing
# if (TRUE) {
#   # Test the prediction function with your example
#   set.seed(123)
#   X <- as.matrix(MASS::Boston[, -ncol(MASS::Boston)])
#   n <- nrow(X)
#   p <- ncol(X) 
#   x <- X
#   y <- as.numeric(MASS::Boston[, ncol(MASS::Boston)])
#   
#   set.seed(123)
#   (idx_train <- sample(seq_len(n), size=floor(0.8*n)))
#   X_train <- x[idx_train, ]
#   X_test <- x[-idx_train, ]
#   y_train <- y[idx_train]
#   y_test <- y[-idx_train]
#   # Fit model
#   model_kernel <- krvfl(X_train, y_train, lambda = 0.1, nb_hidden = 100)
#   
#   # Make predictions
#   pred <- predict(model_kernel, X_test)
#   
#   # Print model summary
#   #print(model_kernel)
#   #summary(model_kernel)
#   
#   # Calculate test RMSE
#   test_rmse <- sqrt(mean((y_test - pred)^2))
#   cat("Test RMSE:", test_rmse, "\n")
# }


