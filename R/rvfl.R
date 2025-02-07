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
  new_predictors <- create_new_predictors(X, 
                                        nb_hidden = n_hidden_features,
                                        nodes_sim = nodes_sim,
                                        activ = activ)
  
  # Fit model using calibmodel
  object <- calibmodel(new_predictors$predictors, y, 
                      engine=engine, lambda=lambda, 
                      positive_response=positive_response,
                      seed=seed, ...)  
  
  # Store parameters needed for prediction
  object$new_predictors <- new_predictors
  object$n_hidden_features <- n_hidden_features
  object$activation_function <- activ  
  object$nodes_sim <- nodes_sim
  
  class(object) <- c("rvfl", "calibmodel")
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
  stopifnot(inherits(object, "rvfl"))
  
  # Create new predictors
  newX <- create_new_predictors(newdata, 
                               nb_hidden = object$n_hidden_features,
                               nodes_sim = object$nodes_sim,
                               activ = object$activation_function, 
                               nn_xm = object$new_predictors$nn_xm, 
                               nn_scales = object$new_predictors$nn_scales)
  
  # Use calibmodel predict method
  class(object) <- "calibmodel"
  pred <- predict(object, newX$predictors, coeff=coeff, ...)
  class(object) <- c("rvfl", "calibmodel")
  
  return(pred)
}


