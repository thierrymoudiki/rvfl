#' Dropout Layer
#'
#' Applies dropout to the input data for regularization during training.
#' The dropout rate determines the probability of retaining a unit.
#'
#' @param x A numeric matrix of input data.
#' @param dropout A numeric value between 0 and 1 representing the dropout rate. Default is 0.
#' @param seed An integer seed for random number generation. Default is 42.
#'
#' @return A matrix of the same shape as `x` with dropout applied.
#' @examples
#' # Apply dropout to a matrix
#' x <- matrix(runif(20), nrow = 5)
#' dropout_layer(x, dropout = 0.5)
dropout_layer <- function(x, dropout = 0, seed = 42) {
  if (dropout == 0) return(x)
  set.seed(seed)
  mask <- matrix(runif(length(x)) > dropout, nrow = nrow(x), ncol = ncol(x))
  return(x * mask / (1 - dropout))  # Scale to maintain expected value
}

#' Custom Regression Function (RVFL2)
#'
#' Fits a regression model using a Random Vector Functional Link (RVFL) network with dropout and regularization.
#' The model includes both the original features and transformed features from a hidden layer.
#'
#' @param x A numeric matrix of input features (training data).
#' @param y A numeric vector of target values (training labels).
#' @param nb_hidden An integer specifying the number of hidden units in the model. Default is 50.
#' @param sigma A numeric value controlling the standard deviation of the initial random weights. Default is 1.
#' @param dropout A numeric value between 0 and 1 representing the dropout rate. Default is 0.
#' @param seed An integer seed for random number generation. Default is 42.
#' @param lambda_1 A numeric value controlling the L2 regularization on the input features. Default is 0.
#' @param lambda_2 A numeric value controlling the L2 regularization on the hidden layer features. Default is 0.
#'
#' @return A list of class `rvfl2` containing the fitted model with the coefficients, hidden layer weights, and training parameters.
rvfl2 <- function(x, y, nb_hidden = 50, sigma = 1, dropout = 0, seed = 42, 
                  lambda_1 = 0, lambda_2 = 0) {
  # Set seed for reproducibility
  set.seed(seed)
  # Store scaling parameters for prediction
  x_center <- apply(x, 2, mean)
  x_scale <- apply(x, 2, sd)
  y_center <- mean(y)
  # Scale and center the input data
  x_scaled <- scale(x)
  y_centered <- scale(y, center = TRUE, scale = FALSE)
  # Ensure nb_hidden is not greater than the number of columns in x_scaled
  ncol_x <- ncol(x_scaled)
  nb_hidden <- min(nb_hidden, ncol_x)
  # Generate random weights for the hidden layer
  W <- matrix(rnorm(nb_hidden * ncol_x, 0, sigma), nrow = ncol_x, ncol = nb_hidden)
  # Compute hidden features Z and apply ReLU activation
  Z <- x_scaled %*% W
  Phi_X <- pmax(Z, 0)
  # Apply dropout to the transformed predictors
  Phi_X <- dropout_layer(Phi_X, dropout = dropout, seed = seed)
  # Create design matrix with original and transformed features
  design_matrix <- cbind(x_scaled, Phi_X)
  # Add regularization to the diagonal
  p <- ncol(design_matrix)
  reg_matrix_1 <- lambda_1 * diag(ncol(x_scaled))
  reg_matrix_2 <- lambda_2 * diag(ncol(Phi_X))
  # Create block diagonal regularization matrix
  reg_matrix <- matrix(0, nrow = p, ncol = p)
  reg_matrix[1:ncol(x_scaled), 1:ncol(x_scaled)] <- reg_matrix_1
  reg_matrix[(ncol(x_scaled)+1):p, (ncol(x_scaled)+1):p] <- reg_matrix_2
  # Solve the regularized least squares problem
  A <- crossprod(design_matrix) + reg_matrix
  b <- crossprod(design_matrix, y_centered)
  # Use solve with error handling
  coef <- try(solve(A, b), silent = TRUE)
  if (inherits(coef, "try-error")) {
    coef <- MASS::ginv(A) %*% b  # Use generalized inverse if solve fails
  }
  
  # Create the model object
  model <- list(
    coef = coef,
    W = W,
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    dropout = dropout,
    seed = seed,
    nb_hidden = nb_hidden,
    sigma = sigma
  )
  
  class(model) <- "rvfl2"
  return(model)
}

#' Make Predictions with RVFL2 Model
#'
#' Predicts the target values using the fitted RVFL2 model on new input data.
#'
#' @param object An object of class `rvfl2` (the fitted model).
#' @param new_x A numeric matrix of input features (test data).
#' @param dropout A numeric value between 0 and 1 representing the dropout rate. Default is the value used during training.
#' @param seed An integer seed for random number generation. Default is the value used during training.
#'
#' @return A numeric vector of predicted target values.
predict.rvfl2 <- function(object, new_x, dropout = object$dropout, seed = object$seed) {
  
  # Scale the new input data using training parameters
  new_x_scaled <- scale(new_x, center = object$x_center, scale = object$x_scale)
  
  # Compute hidden features for new data
  Z_new <- new_x_scaled %*% object$W
  Phi_X_new <- pmax(Z_new, 0)
  
  # Apply dropout
  Phi_X_new <- dropout_layer(Phi_X_new, dropout = dropout, seed = seed)
  
  # Create design matrix for new data
  design_matrix_new <- cbind(new_x_scaled, Phi_X_new)
  
  # Get predictions
  predictions <- design_matrix_new %*% object$coef + object$y_center
  
  return(predictions)
}

#' Summary of RVFL2 Model
#'
#' Displays a summary of the RVFL2 model, including the key hyperparameters.
#'
#' @param object An object of class `rvfl2` (the fitted model).
#'
#' @return None (prints summary information to the console).
summary.rvfl2 <- function(object) {
  cat("RVFL2 Model Summary\n")
  cat("====================\n")
  cat("Lambda_1: ", object$lambda_1, "\n")
  cat("Lambda_2: ", object$lambda_2, "\n")
  cat("Dropout Rate: ", object$dropout, "\n")
  cat("Number of hidden units: ", object$nb_hidden, "\n")
  cat("Sigma: ", object$sigma, "\n")
}

#' Compute Residuals for RVFL2 Model
#'
#' Computes the residuals (observed - predicted) for the RVFL2 model on the given data.
#'
#' @param object An object of class `rvfl2` (the fitted model).
#' @param x A numeric matrix of input features (test data).
#' @param y A numeric vector of observed target values (test labels).
#'
#' @return A numeric vector of residuals.
residuals.rvfl2 <- function(object, x, y) {
  predictions <- predict(object, x)
  residuals <- y - predictions
  return(residuals)
}

#' Compute Fitted Values for RVFL2 Model
#'
#' Returns the fitted (predicted) values for the RVFL2 model on the given data.
#'
#' @param object An object of class `rvfl2` (the fitted model).
#' @param x A numeric matrix of input features (test data).
#'
#' @return A numeric vector of fitted values.
fitted.rvfl2 <- function(object, x) {
  return(predict(object, x))
}
