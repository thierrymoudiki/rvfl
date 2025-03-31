#' Fit and forecast a time series using RVFL
#'
#' This function combines model fitting using tsrvfl() and forecasting in one step.
#' It's a convenience wrapper that fits a model and immediately generates forecasts.
#'
#' @param y A numeric vector or time series of the response variable.
#' @param h The number of steps ahead to forecast.
#' @param lags The number of lags to use in the model.
#' @param engine The engine to use for the model fitting.
#' @param level The confidence level for the prediction intervals.
#' @param n_hidden_features The number of hidden features to use in the model.
#' @param nodes_sim The method to use for simulating nodes.
#' @param activ The activation function to use in the model.
#' @param type_pi The type of prediction interval to use.
#' @param B The number of bootstrap samples to use.
#' @param agg The aggregation function to use.
#' @param seed The random seed to use.
#' @param coeffs The coefficients to use in the model.
#'
#' @return An object of class "forecast" containing the fitted model and forecasts
#'
#' @examples
#' y <- ts(rnorm(120,0,3) + 1:120 + 20*sin(2*pi*(1:120)/12), frequency=12)
#' fit <- tsrvflf(y, h=20)
#' #plot(fit)
#'
tsrvflf <- function(y, h = 5, 
                    lags = 15L,                
                    engine=stats::lm,
                    level = 95, 
                    n_hidden_features=5L, 
                    nodes_sim = c("sobol", "halton", "unif"),
                    activ = c("relu", "sigmoid", "tanh", 
                    "leakyrelu", "elu", "linear"),                    
                    type_pi = c("surrogate", "kde", "bootstrap"),
                    B = 250L, agg = c("mean", "median"), 
                    seed = 123, coeffs = NULL,
                    ...)
{
  set.seed(seed)  
  n <- length(y)
  freq_x <- frequency(y)
  half_n <- ceiling(n/2)
  idx_train <- seq_len(half_n)
  idx_calib <- setdiff(seq_len(n), idx_train)
  splitted_y <- misc::splitts(y)
  y_train <- splitted_y$training
  y_calib <- splitted_y$testing
  type_pi <- match.arg(type_pi)
  agg <- match.arg(agg)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  fit_func <- function(x, y, ...) {
    rvfl::rvfl(x, y, 
               engine = engine,
               n_hidden_features = n_hidden_features, 
               nodes_sim = nodes_sim, 
               activ = activ, 
               ...)
  }  
  y_pred_calibration <- ml_forecast(y = y_train, 
                            h = length(y_calib), 
                            lags=lags, 
                            fit_func = fit_func,
                            predict_func = predict.rvfl,
                            coeffs = coeffs,
                            ...)$mean
  #misc::debug_print(y_calib)
  #misc::debug_print(y_pred_calibration)
  preds_obj <- ml_forecast(y = y_calib, 
                           h = h,  
                           lags = lags, 
                           fit_func = fit_func,
                           predict_func = predict.rvfl,
                           coeffs = coeffs,
                           ...) 
  #misc::debug_print(preds_obj)
  preds <- preds_obj$mean

  tspx <- tsp(y_calib)
  start_preds <- tspx[2] + 1 / tspx[3]                         
  matrix_preds <- replicate(B, preds)     
  #misc::debug_print(y_calib)
  #misc::debug_print(y_pred_calibration)
  calibrated_raw_residuals <- y_calib - y_pred_calibration
  #misc::debug_print(calibrated_raw_residuals)
  scaled_calib_resids <- base::scale(calibrated_raw_residuals)
  xm <- attr(scaled_calib_resids, "scaled:center")
  xsd <- attr(scaled_calib_resids, "scaled:scale")
  scaled_calibrated_residuals <- base::scale(calibrated_raw_residuals,
                                             center = TRUE,
                                             scale = TRUE)
  #misc::debug_print(scaled_calibrated_residuals)
  if (type_pi == "kde") {        
        simulated_scaled_calibrated_residuals <-
            rgaussiandens(
            scaled_calibrated_residuals,
            n = h,
            p = B,
            seed = seed
            )
        sd_calibrated_residuals <- sd(calibrated_raw_residuals)
    }

    if (type_pi == "surrogate") {
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
            tseries::surrogate(scaled_calibrated_residuals,
                                    ns =
                                        B)[seq_along(h), ]

        sd_calibrated_residuals <- sd(calibrated_raw_residuals)      
    }

    if (type_pi == "bootstrap") {
      freq_calibrated_raw_residuals <- frequency(calibrated_raw_residuals)
      if (length(calibrated_raw_residuals) <= 2 * freq_calibrated_raw_residuals)
        freq_calibrated_raw_residuals <- 1L
      block_size <-
        ifelse(
          freq_calibrated_raw_residuals > 1,
          2 * freq_calibrated_raw_residuals,
          min(8, floor(
            length(calibrated_raw_residuals) / 2
          ))
        )
      block_size <-
        floor(min(
          max(3L, block_size),
          length(calibrated_raw_residuals) - 1L
        ))
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
          tseries::tsbootstrap(
            scaled_calibrated_residuals,
            nb =
              B,
            b = floor(block_size),
            type =
              "block"
          )[seq_along(h), ]
        sd_calibrated_residuals <- sd(calibrated_raw_residuals)      
    }

    sims <- matrix_preds 
    sims <- sims + sd_calibrated_residuals * simulated_scaled_calibrated_residuals

    sims <- ts(sims,
               start = start_preds,
               frequency = frequency(y_train))
    preds_lower <-
      apply(sims, 1, function(x)
        quantile(x, probs = (1 - level / 100) / 2))
    preds_upper <-
      apply(sims, 1, function(x)
        quantile(x, probs = 1 - (1 - level / 100) / 2))

    out <- vector("list", 8) 
    class(out) <- "forecast"
    out$mean <- ts(switch(
      agg,
      median = apply(sims, 1, median),
      mean = apply(sims, 1, mean)
    ),
    start = start_preds,
    frequency = freq_x)    
    out$lower <- ts(preds_lower,
                    start = start_preds,
                    frequency = freq_x)
    out$upper <- ts(preds_upper,
                    start = start_preds,
                    frequency = freq_x)
    out$x <- y_calib
    out$level <- level 
    out$method <- "conformalized ML"
    out$model <- preds_obj$model
    out$residuals <- ts(
      calibrated_raw_residuals,
      start = start(y_calib),
      frequency = frequency(y_train)
    ) # /!\ not the same residuals, beware
    out$fitted <- ts(
      y_pred_calibration,
      start = start(y_calib),
      frequency = frequency(y_train)
    ) # /!\ not the same fitted, beware
    out$sims <- sims
    return(out)
}

ml_forecast <- function(y, h, 
                        lags=1, 
                        fit_func = rvfl::rvfl,
                        predict_func = predict.rvfl,
                        coeffs = NULL, 
                        ...)
{
    df <- as.data.frame(embed(rev(as.numeric(y)), lags + 1L))
    ncol_df <- ncol(df)
    colnames(df) <- c(paste0("lag", rev(seq_len(lags))), "y") 
    
    if(is.null(coeffs))
    {            
      # Fit model using formula interface
      fit <- try({
          model <- fit_func(y ~ ., data = df, ...)
          if (is.null(model$coefficients) && is.null(model$model)) {
              stop("Model fitting failed")
          }
          model
      }, silent=TRUE)
      
      if (inherits(fit, "try-error"))
      {
          # Try matrix interface if formula interface fails
          x_matrix <- as.matrix(df)[, -ncol_df]
          y_vector <- df$y
          fit <- try({
              model <- fit_func(x = x_matrix, y = y_vector, ...)
              if (is.null(model$coefficients) && is.null(model$model)) {
                  stop("Model fitting failed")
              }
              model
          }, silent=TRUE)
          
          if (inherits(fit, "try-error")) {
              stop("Both formula and matrix interfaces failed to fit model")
          }
      }
      
      # Get the most recent lags and their differences
      latest_lags <- rev(as.numeric(y)[1:lags])
      last_diff <- diff(tail(y, 2))[1]
      
      predictions <- numeric(h)
      for (i in 1:h)
      {
          # Create prediction data frame
          newdata <- data.frame(matrix(latest_lags, nrow=1))
          colnames(newdata) <- paste0("lag", rev(seq_len(lags)))
          
          # Make prediction
          pred <- try({
              if (inherits(fit, "lm")) {
                  predict(fit, newdata)
              } else {
                  predict_func(fit, as.matrix(newdata))
              }
          }, silent=TRUE)
          
          if (inherits(pred, "try-error") || is.na(pred)) {
              warning("Prediction error at step ", i)
              # Instead of using mean, extrapolate using last difference
              predictions[i] <- latest_lags[1] + last_diff
          } else {
              predictions[i] <- as.numeric(pred)[1]
          }
          
          # Update lags for next prediction
          if (i < h) {
              latest_lags <- c(latest_lags[-1], predictions[i])
              # Update last_diff based on new prediction
              last_diff <- predictions[i] - latest_lags[1]
          }
      }
    } else {
      latest_lags <- rev(as.numeric(y)[1:lags])
      predictions <- numeric(h)
      
      for (i in 1:h)
      {
          newdata <- matrix(latest_lags, nrow=1)
          colnames(newdata) <- paste0("lag", rev(seq_len(lags)))
          
          predictions[i] <- as.numeric(newdata %*% coeffs)
          
          if (i < h) {
              latest_lags <- c(latest_lags[-1], predictions[i])
          }
      }
    }    
    
    return(list(model = fit, mean = predictions))
}
