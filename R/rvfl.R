#' @export
rvfl <- function(X, y, 
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
  new_predictors <- create_new_predictors(X, 
                                          nb_hidden = n_hidden_features,
                                          nodes_sim = nodes_sim,
                                          activ = activ)
  # Fit model using calibmodel
  object <- rvfl::calibmodel(new_predictors$predictors, y, 
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
  X <- model.matrix(formula, data)[,-1, drop=FALSE]  # Remove intercept column  
  rvfl::rvfl(X, y, ...)
}

#' @export
predict.rvfl <- function(object, newdata, 
                               method=c("none", "gaussian", "surrogate", 
                               "bootstrap", "tsbootstrap"), 
                               level=95, nsim=100, seed=123, ...) {
  return(predict.calibmodel(object, newdata, method, level, nsim, seed, ...))
}

#' @export
coef.rvfl <- function(object) {
  return(coef.calibmodel(object))
}

#' @export
simulate.rvfl <- function(object, newdata, nsim = 100, 
                                method = c("gaussian", "surrogate", "bootstrap", "tsbootstrap"),
                                seed = 123, ...) {
  return(simulate.calibmodel(object, newdata, nsim, method, seed, ...))
}


#' @export
summary.rvfl <- function(object, newdata=NULL) {
    # Get original summary
  orig_summary <- summary(object$model)  
  # Calculate numerical derivatives
  if (is.null(newdata)) {
    X <- model.matrix(object$model)
  } else {
    X <- model.matrix(object$model, newdata)
  }
  y <- model.response(model.frame(object$model))
  n <- nrow(X)
  p <- ncol(X) - object$n_hidden_features
  h <- 1e-5  # Small step size for numerical derivatives  
  # Initialize matrices for derivatives
  derivatives <- matrix(0, nrow=n, ncol=p)  
  colnames(derivatives) <- colnames(X)[1:p]
  # Calculate derivatives for each predictor
  for(j in 1:p) {
    X_plus <- X
    X_minus <- X
    X_plus[,j] <- X_plus[,j] + h
    X_minus[,j] <- X_minus[,j] - h    
    pred_plus <- predict(object, X_plus)
    pred_minus <- predict(object, X_minus)    
    derivatives[,j] <- (pred_plus - pred_minus)/(2*h)
  }  
  # Calculate mean effects and standard errors
  mean_effects <- colMeans(derivatives)
  se_effects <- apply(derivatives, 2, sd)/sqrt(n)
  
  # Perform t-tests for each predictor
  t_test_results <- try(apply(derivatives, 2, function(x) {
    test <- t.test(x)
    c(
      Effect = test$estimate,
      Std.Error = test$stderr,
      p.value = test$p.value,
      CI.lower = test$conf.int[1],
      CI.upper = test$conf.int[2]
    )
  }), silent=TRUE)
  if (inherits(t_test_results, "try-error")) {
    t_test_results <- t(apply(derivatives, 2, summary))
  } else {
    t_test_results <- t(t_test_results)
  }
  # Create numerical derivatives summary table
  deriv_summary <- as.data.frame(t_test_results)      
  # Return both summaries
  list(
    model_summary = orig_summary,
    derivatives_summary = deriv_summary
  )
}

#' @export
plot.rvfl <- function(object) {
  return(plot.calibmodel(object))
}