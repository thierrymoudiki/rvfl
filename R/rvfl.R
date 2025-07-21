# Define the generic function
#' @title rvfl
#' @description Generic function for RVFL models
#' @param x Object input
#' @param ... Additional arguments passed to methods
#' @export
rvfl <- function(x, ...) {
  UseMethod("rvfl")
}

#' @export
rvfl.matrix <- function(x, y, 
                 engine=stats::lm,
                 lambda=10**seq(-10, 10, length.out=100), 
                 n_hidden_features=5L, 
                 nodes_sim = c("sobol", "halton", "unif"),
                 activ = c("relu", "sigmoid", "tanh", 
                 "leakyrelu", "elu", "linear"),
                 positive_response=FALSE,
                 seed=123, ...) {                
  nodes_sim <- match.arg(nodes_sim)                  
  activ <- match.arg(activ) 
  set.seed(seed)
  # Create hidden features
  new_predictors <- create_new_predictors(x, 
                                          nb_hidden = n_hidden_features,
                                          nodes_sim = nodes_sim,
                                          activ = activ)
  ##misc::debug_print(new_predictors)                                          
  # Fit model using calibmodel
  object <- rvfl::calibmodel(x=new_predictors$predictors, y=y, 
                             engine=engine, lambda=lambda, 
                             positive_response=positive_response,
                             seed=seed, ...)  
  # Store parameters needed for prediction
  object$new_predictors <- new_predictors
  object$n_hidden_features <- n_hidden_features
  object$activ <- activ  
  object$nodes_sim <- nodes_sim
  class(object) <- c("rvfl", "calibmodel")
  return(object)
}

#' @export
rvfl.formula <- function(formula, data, ...) {
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(formula, data)[,-1, drop=FALSE]  # Remove intercept column  
  ##misc::debug_print(y)
  ##misc::debug_print(x)
  rvfl.matrix(x, y, ...)
}

#' @export
predict.rvfl <- function(object, newdata, 
                               method=c("none", "gaussian", "surrogate", 
                               "bootstrap", "tsbootstrap"), 
                               level=95, nsim=100, seed=123, ...) {
  method <- match.arg(method)  
  # Input validation
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    newdata <- as.matrix(newdata)
  }  
  misc::debug_print(newdata)
  # Create hidden features with error checking
  #tryCatch({
  new_predictors <- create_new_predictors(newdata, 
                                          nb_hidden = object$n_hidden_features,
                                          nodes_sim = object$nodes_sim,
                                          activ = object$activ)    
  misc::debug_print(new_predictors)                                      
  # Validate the created features
  if (any(is.na(new_predictors$predictors))) {
    stop("Hidden features contain NA/NaN values")
  }    
  # Make prediction
  pred <- predict.calibmodel(object, new_predictors$predictors,
                              method, level, nsim, seed, ...)
  if (length(pred) == 0 || any(is.na(pred))) {
    stop("Invalid prediction result")
  }    
  return(pred)
  #}, error = function(e) {
  #  warning("Error in prediction: ", e$message)
  #  return(rep(NA, nrow(newdata)))
  #})
}

#' @export
coef.rvfl <- function(object) {
  res <- try(coef.calibmodel(object), silent=TRUE)
  if (inherits(res, "try-error")) {
    return(NULL)
  }
  return(res)
}

#' @export
simulate.rvfl <- function(object, newdata, nsim = 250L, 
                                method = c("gaussian", "surrogate", 
                                "bootstrap", "tsbootstrap", "splitconformal"),
                                seed = 123, ...) {
  method <- match.arg(method)
  new_predictors <- create_new_predictors(newdata, 
                                          nb_hidden = object$n_hidden_features,
                                          nodes_sim = object$nodes_sim,
                                          activ = object$activ)
  return(simulate.calibmodel(object, new_predictors$predictors, nsim, method, seed, ...))
}


#' @export
summary.rvfl <- function(object, newdata) {
  if (is.null(newdata)) {
    stop("newdata must be provided to calculate derivatives.")
  }
  
  # Validate newdata
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("newdata must be a numeric matrix or data frame.")
  }
  if (any(is.na(newdata))) {
    stop("newdata contains NA/NaN values, which are not allowed.")
  }
  
  # Initialize summary list
  rvfl_summary <- list()
  
  # Calculate numerical derivatives
  zero <- 1e-4
  eps_factor <- zero ** (1 / 3)
  n_vars <- ncol(newdata)
  derivatives <- matrix(0, nrow = nrow(newdata), ncol = n_vars)
  colnames(derivatives) <- colnames(newdata)
  
  for (j in 1:n_vars) {
    newdata_plus <- newdata_minus <- newdata
    value_x <- newdata[, j]
    cond <- abs(value_x) > zero
    h <- ifelse(cond, eps_factor * abs(value_x), zero)
    newdata_plus[, j] <- newdata_plus[, j] + h
    newdata_minus[, j] <- newdata_minus[, j] - h
    
    # Get predictions for perturbed inputs
    pred_plus <- predict(object, newdata_plus)
    pred_minus <- predict(object, newdata_minus)
    
    # Calculate centered difference
    derivatives[, j] <- (pred_plus - pred_minus) / (2 * h)
  }
  
  # Derivative statistics with confidence intervals and significance codes
  derivative_stats <- data.frame(
    #Variable = colnames(derivatives),
    Mean = apply(derivatives, 2, mean),
    StdDev = apply(derivatives, 2, sd),
    #Min = apply(derivatives, 2, min),
    #Max = apply(derivatives, 2, max),
    CI_Lower = NA,
    CI_Upper = NA,
    P_Value = NA,
    Significance = ""
  )
  
  for (j in 1:n_vars) {
    # Use tryCatch to handle constant data or other errors in t.test
    t_test_result <- tryCatch(
      t.test(derivatives[, j]),
      error = function(e) NULL
    )
    
    if (is.null(t_test_result) || is.na(t_test_result$p.value)) {
      # Assign NA for constant data or errors
      derivative_stats$CI_Lower[j] <- NA
      derivative_stats$CI_Upper[j] <- NA
      derivative_stats$P_Value[j] <- NA
      derivative_stats$Significance[j] <- " "
    } else {
      # Extract t-test results for non-constant data
      derivative_stats$CI_Lower[j] <- t_test_result$conf.int[1]
      derivative_stats$CI_Upper[j] <- t_test_result$conf.int[2]
      derivative_stats$P_Value[j] <- t_test_result$p.value
      
      # Assign significance codes
      if (!is.na(t_test_result$p.value)) {
        if (t_test_result$p.value < 0.001) {
          derivative_stats$Significance[j] <- "***"
        } else if (t_test_result$p.value < 0.01) {
          derivative_stats$Significance[j] <- "**"
        } else if (t_test_result$p.value < 0.05) {
          derivative_stats$Significance[j] <- "*"
        } else if (t_test_result$p.value < 0.1) {
          derivative_stats$Significance[j] <- "."
        } else {
          derivative_stats$Significance[j] <- " "
        }
      }
    }
  }
  
  # Add derivative statistics to the summary
  rvfl_summary$derivative_stats <- derivative_stats
  
  # Add model coefficients
  rvfl_summary$coefficients <- coef(object)
  
  # Add additional model information
  rvfl_summary$n_hidden_features <- object$n_hidden_features
  rvfl_summary$activation_function <- object$activ
  rvfl_summary$nodes_simulation <- object$nodes_sim
  
  # Add class for printing
  class(rvfl_summary) <- "summary.rvfl"
  
  return(rvfl_summary)
}

#' @export
print.summary.rvfl <- function(x, ...) {
  cat("Random Vector Functional Link (RVFL) Model Summary\n")
  cat("--------------------------------------------------\n")
  cat("Number of Hidden Features:", x$n_hidden_features, "\n")
  cat("Activation Function:", x$activation_function, "\n")
  cat("Node Simulation Method:", x$nodes_simulation, "\n\n")    
  cat("\nDerivative Statistics:\n")
  print(x$derivative_stats)
}

#' @export
plot.rvfl <- function(object) {
  return(plot.calibmodel(object))
}

#' @export
residuals.rvfl <- function(object) {
  return(residuals.calibmodel(object))
}

#' @export
fitted.rvfl <- function(object) {
  return(fitted.calibmodel(object))
}