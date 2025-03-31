#' @export
create_data_partition <- function(y, times = 1, p = 0.5, list = TRUE, groups = min(5, length(y))) {
  # Input validation
  if (length(y) < 2) stop("y must have at least 2 data points")
  if (groups < 2) groups <- 2
  if (p <= 0 || p >= 1) stop("p must be between 0 and 1")
  if (times <= 0 || !is.numeric(times)) stop("times must be a positive integer")  
  # Handle numeric y by binning into groups
  if (is.numeric(y)) {
    y <- cut(y, breaks = unique(quantile(y, probs = seq(0, 1, length.out = groups))), include.lowest = TRUE)
  } else if (is.factor(y)) {
    # Handle factor y by checking for empty or single-record classes
    xtab <- table(y)
    if (any(xtab == 0)) {
      warning(paste("Some classes have no records (", paste(names(xtab)[xtab == 0], sep = "", collapse = ", "), ") and these will be ignored"))
      y <- factor(as.character(y))
    }
    if (any(xtab == 1)) {
      warning(paste("Some classes have a single record (", paste(names(xtab)[xtab == 1], sep = "", collapse = ", "), ") and these will be selected for the sample"))
    }
  }  
  # Function to subsample indices for each group
  subsample <- function(dat, p) {
    if (nrow(dat) == 1) {
      return(dat$index)
    } else {
      num <- ceiling(nrow(dat) * p)
      return(sample(dat$index, size = num))
    }
  }  
  # Create partitions
  out <- vector("list", times)
  for (j in 1:times) {
    # Split data by y and subsample
    tmp <- split(data.frame(y = y, index = seq_along(y)), y)
    tmp <- lapply(tmp, subsample, p = p)
    tmp <- sort(unlist(tmp))
    out[[j]] <- tmp
  }  
  # Format output
  if (!list) {
    out <- matrix(unlist(out), ncol = times)
    colnames(out) <- paste0("Resample", 1:times)
  } else {
    names(out) <- paste0("Resample", 1:times)
  }  
  return(out)
}

remove_zero_cols <- function(x)
{
  x[, colSums(x == 0) != nrow(x)]
}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

my_sd <- function(x) {
  n <- dim(x)[1]
  return(drop(rep(1 / (n - 1), n) %*% (x - tcrossprod(
    rep.int(1, n), colMeans(x)
  )) ^ 2) ^ 0.5)
}
my_sd <- compiler::cmpfun(my_sd)

my_scale <- function(x, xm = NULL, xsd = NULL) {
  rep_1_n <- rep.int(1, dim(x)[1])
  
  # centering and scaling, returning the means and sd's
  if (is.null(xm) && is.null(xsd)) {
    xm <- colMeans(x)
    xsd <- my_sd(x)
    return(list(
      res = (x - tcrossprod(rep_1_n, xm)) / tcrossprod(rep_1_n, xsd),
      xm = xm,
      xsd = xsd
    ))
  }
  
  # centering and scaling
  if (is.numeric(xm) && is.numeric(xsd)) {
    return((x - tcrossprod(rep_1_n, xm)) / tcrossprod(rep_1_n, xsd))
  }
  
  # centering only
  if (is.numeric(xm) && is.null(xsd)) {
    return(x - tcrossprod(rep_1_n, xm))
  }
  
  # scaling only
  if (is.null(xm) && is.numeric(xsd)) {
    return(x / tcrossprod(rep_1_n, xsd))
  }
}
my_scale <- compiler::cmpfun(my_scale)

# create new predictors
#' @export
create_new_predictors <- function(x, nb_hidden = 5,
                                  nodes_sim = c("sobol", "halton", "unif"),
                                  activ = c("relu", "sigmoid", "tanh",
                                            "leakyrelu", "elu", "linear"),
                                  nn_xm = NULL, nn_scales = NULL)
{
  x <- as.matrix(x)
  if (identical(nb_hidden, 0))
  {
    return(x)
  }
  g <- switch(match.arg(activ),
              "relu" = function(x) x*(x>0),
              "sigmoid" = function(x) (1/(1 + exp(-x)))*(x >= 0) + (exp(x)/(1 + exp(x)))*(x < 0),
              "tanh" = function(x) tanh(x),
              "leakyrelu" = function(x) x*(x > 0) + 0.01*x*(x <= 0),
              "elu" = function(x) x*(x >= 0) + 0.01*(exp(x)-1)*(x < 0),
              "linear" = function(x) x)
  
  p <- ncol(x)
  set.seed(1)
  w <- remove_zero_cols(switch(match.arg(nodes_sim),
                               "sobol" = 2*t(randtoolbox::sobol(nb_hidden + 1, p)) - 1,
                               "halton" = 2*t(randtoolbox::halton(nb_hidden, p)) - 1,
                               "unif" = matrix(stats::runif(nb_hidden*p, min = -1, max = 1),
                                               nrow = p, ncol = nb_hidden)))
  
  if((!is.null(nn_xm) && is.null(nn_scales)) || (is.null(nn_xm) && !is.null(nn_scales)))
    stop("either nn_xm and nn_scales provided, or both left to NULL")
  
  if (is.null(nn_xm) && is.null(nn_scales))
  {
    scaled_x <- my_scale(x)
    hidden <- g(scaled_x$res%*%w)
    res <- cbind(x, hidden)
    
    if (length(colnames(x)) > 0){
      colnames(res) <- c(colnames(x),
                         paste0("h", 1:ncol(hidden)))
    } else {
      colnames(res) <- c(paste0("x", 1:ncol(x)),
                         paste0("h", 1:ncol(hidden)))
    }
    
    return(list(nn_xm = scaled_x$xm,
                nn_scales = scaled_x$xsd,
                w = w, predictors = res))
  }
  
  if (!is.null(nn_xm) && !is.null(nn_scales))
  {
    stopifnot(length(nn_xm) == ncol(x) || length(nn_scales) == ncol(x))
    scaled_x <- my_scale(as.matrix(x),
                         xm = as.vector(nn_xm),
                         xsd = as.vector(nn_scales))
    hidden <- g(as.matrix(scaled_x)%*%w)
    res <- cbind(x, hidden)
    
    if (length(colnames(x)) > 0){
      colnames(res) <- c(colnames(x),
                         paste0("h", 1:ncol(hidden)))
    } else {
      colnames(res) <- c(paste0("x", 1:ncol(x)),
                         paste0("h", 1:ncol(hidden)))
    }
    
    return(list(w = w, predictors = res))
  }
}

toMat <- function(x) {
  if (NCOL(x) > 1 && !is.matrix(x)) {
    x <- matrix(x, ncol = NCOL(x))
  }
  return(x)
}

# Converts arguments into data.frame, whilst retaining mts/ts/matrix properties
datamat <- function(..., flatten=TRUE, functions=TRUE) {
  vars <- list(...)
  if (length(vars) == 0) {
    return(data.frame())
  }
  if (!is.null(names(vars))) {
    names(vars)[!nzchar(names(vars))] <- as.character(substitute(list(...))[-1])[!nzchar(names(vars))]
  }
  else {
    names(vars) <- as.character(substitute(list(...))[-1])
  }
  if (flatten) {
    i <- 1
    while (i <= length(vars)) {
      if (is.data.frame(vars[[i]])) {
        vars <- c(vars, c(vars[[i]])) # Append data.frame components
        vars[[i]] <- NULL # Remove data.frame
      }
      else if (is.matrix(vars[[i]])) {
        for (j in 1:NCOL(vars[[i]])) {
          vars[[length(vars) + 1]] <- vars[[i]][, j]
          names(vars)[length(vars)] <- make.names(colnames(vars[[i]])[j])
        }
        i <- i + 1
      }
      else {
        i <- i + 1
      }
    }
  }
  class(vars) <- "data.frame"
  row.names(vars) <- 1:max(sapply(vars, NROW))
  #   if(is.ts(vars[,1])){
  #     if(NCOL(vars)>1){
  #       class(vars) <- c(class(vars),"mts")
  #     }
  #     class(vars) <- c(class(vars),"ts")
  #     tspx <- unique(sapply(vars,tsp), MARGIN = 2)
  #     if(length(tspx)==3){
  #       attr(vars, "tsp") <- tspx
  #     }
  #   }
  return(vars)
}

recoverTSP <- function(times.x) {
  freq <- sort(unique(round(times.x %% 1, digits = 6))) # The subset cannot increase frequency
  freq <- length(freq)
  return(c(min(times.x), min(times.x) + (length(times.x) - 1) / freq, freq))
}


# simulate Gaussian density -----
rgaussiandens <- function(x,
                          n = length(x),
                          p = 1,
                          seed = 32124) {
  if(any(is.na(x)))
  {
    x <- approx(x, method = "linear", n = length(x))$y
  }
  z <- try(stats::density(x, bw = "sj",
                          kernel = "gaussian"),
           silent = TRUE)
  if (inherits(z, "try-error")) {
    z <- density(x, kernel = "gaussian")
  }

  width <- z$bw                              # Kernel width
  rkernel <-
    function(n, seed) {
      set.seed(seed)
      stats::rnorm(n, sd = width)
    }  # Kernel sampler
  if (p <= 1)
  {
    set.seed(seed)
    return(sample(x, n, replace = TRUE) + rkernel(n, seed))    # Here's the entire algorithm
  } else {
    return(sapply(1:p,
                  function(i) {
                    set.seed(seed + i - 1)
                    sample(x, n, replace = TRUE) + rkernel(n, seed + i - 1)
                  }))
  }
}
