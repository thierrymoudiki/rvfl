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
  #misc::debug_print(new_predictors)                                          
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
  #misc::debug_print(y)
  #misc::debug_print(x)
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
  #misc::debug_print(newdata)
  # Create hidden features with error checking
  #tryCatch({
  new_predictors <- create_new_predictors(newdata, 
                                        nb_hidden = object$n_hidden_features,
                                        nodes_sim = object$nodes_sim,
                                        activ = object$activ)    
  #misc::debug_print(new_predictors)                                      
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
  orig_summary <- list()
  #orig_summary$model_summary <- summary(object$model)
  #misc::debug_print(newdata)
  # If newdata is provided, calculate numerical derivatives
  zero <- 1e-4
  eps_factor <- zero ** (1 / 3)  
  # Small perturbation for centered differences
  #h <- 1e-6
  n_vars <- ncol(newdata)
  derivatives <- matrix(0, nrow=nrow(newdata), ncol=n_vars)
  colnames(derivatives) <- colnames(newdata)
  # Calculate centered differences for each variable
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
    # Add derivatives to the summary
    #orig_summary$derivatives <- derivatives
  return(summary(derivatives))  
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