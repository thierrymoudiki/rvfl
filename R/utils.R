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