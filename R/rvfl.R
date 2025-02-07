#' @export
rvfl <- function(X, y, 
                 engine=stats::lm,
                 lambda=10**seq(-10, 10, length.out=100), 
                 n_hidden_features=5L, 
                 nodes_sim = c("sobol", "halton", "unif"),
                 activ = c("relu", "sigmoid", "tanh", "leakyrelu", "elu", "linear"),
                 positive_response=FALSE,
                 seed=123, ...) {                
  nodes_sim <- match.arg(nodes_sim)                  
  activ <- match.arg(activ) 
  set.seed(seed)
  
  # Create hidden features
  #misc::debug_print("Creating hidden features:")
  #misc::debug_print(paste("Activation function:", activ))
  
  new_predictors <- create_new_predictors(X, 
                                        nb_hidden = n_hidden_features,
                                        nodes_sim = nodes_sim,
                                        activ = activ)
  
  #misc::debug_print("Hidden layer values:")
  #misc::debug_print(head(new_predictors$predictors[, (ncol(X) + 1):ncol(new_predictors$predictors)]))
  
  # Fit model using calibmodel
  object <- calibmodel(new_predictors$predictors, y, 
                      engine=engine, lambda=lambda, 
                      positive_response=positive_response,
                      seed=seed, ...)  
  
  # Store parameters needed for prediction
  object$new_predictors <- new_predictors
  object$n_hidden_features <- n_hidden_features
  object$activ <- activ  
  object$nodes_sim <- nodes_sim
  
  class(object) <- "rvfl"
  return(object)
}

#' @export
rvfl.formula <- function(formula, data, ...) {
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)[,-1, drop=FALSE]  # Remove intercept column  
  rvfl(X, y, ...)
}

#' @export
predict.rvfl <- function(object, newdata, coeff=NULL, ...) {
  newX <- create_new_predictors(newdata, 
                               nb_hidden = object$n_hidden_features,
                               nodes_sim = object$nodes_sim,
                               activ = object$activ, 
                               nn_xm = object$new_predictors$nn_xm, 
                               nn_scales = object$new_predictors$nn_scales)
   
  # Use calibmodel predict method
  #class(object) <- "calibmodel"
  pred <- predict.calibmodel(object, newX$predictors, coeff=coeff, ...)
  #class(object) <- c("rvfl", "calibmodel")
  return(pred)
}

#' @export
simulate.rvfl <- function(object, newdata, nsim = 100, 
                          method = c("surrogate", "bootstrap", "tsbootstrap"),
                          seed = 123, ...) {
  
  set.seed(seed)  
  method <- match.arg(method)
  #misc::debug_print(object)
  #misc::debug_print(class(object))
  # Try prediction with data frame first
  pred <- try(drop(predict.rvfl(object, as.data.frame(newdata), coeff=coeff, ...)), silent = TRUE)
  if (inherits(pred, "try-error")) {
    pred <- drop(predict.rvfl(object, as.matrix(newdata), coeff=coeff, ...))
  }
  #misc::debug_print(pred)
  if (method == "gaussian") { 
    stopifnot(object$positive_response == FALSE)   
    # Get the residual standard error from the model
    sigma <- try(sqrt(sum(object$residuals^2) / object$df.residual), silent=TRUE)
    if (inherits(sigma, "try-error")) {
      # Use number of predictors from the original model instead of newdata
      sigma <- try(sqrt(sum(object$residuals^2) / (length(object$residuals) - ncol(newdata) - 1)), silent=TRUE)
    }
    # Generate random normal errors
    errors <- matrix(rnorm(length(pred) * nsim, 
                           mean = 0, 
                           sd = sigma), 
                     nrow = length(pred), 
                     ncol = nsim)
  }
  
  if (method == "surrogate") {
    errors <- tseries::surrogate(sample(object$residuals, size=length(pred), replace=TRUE), ns=nsim)
  }
  
  if (method == "bootstrap") {
    errors <- matrix(sample(object$residuals, size=length(pred)*nsim, replace=TRUE), nrow=length(pred), ncol=nsim)    
  }
  
  if (method == "tsbootstrap") {
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


