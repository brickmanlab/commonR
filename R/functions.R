#' Non-zero mean of data frames, matrices or vectors
#'
#' Calculate non zero mean of a data frame, matrix or vector. If the parameter \code{by} is applied,
#' mean will be calculated per groups of columns, not rows.
#' @param x a numerical vector, matrix or data frame.
#' @param by a vector or list of grouping elements, each as long as the variables in the data frame
#' \code{x}. The elements are coerced to factors before use.
#' @return Non zero mean of \code{x} or non zero mean of groups of columns in \code{x}.
#' @examples
#' vector_test <- c(rep(0,10),1:100)
#' nz_mean(vector_test)
#' matrix_test <- matrix(c(rep(0,10),1:100),nrow=10)
#' nz_mean(matrix_test, by = c(1,2,3,1,2,3,1,2,3,1,2))
#' df_test <- as.data.frame(matrix_test)
#' nz_mean(df_test, by = c(1,2,3,1,2,3,1,2,3,1,2))
#'
#' @export
nz_mean <- function(x, by = NULL) {
  # Function to calculate the non zero mean of a vector or data frame/matrix.
  # In case of a data frame/matrix, it can be the non zero mean can be calculated by
  # groups of columns using a vector that specify which column belongs to which group

  calculate_nz_mean <- function(vector) {
    # Actual function to calculate the non zero mean
    return(mean(vector[vector != 0]))
  }

  if (is.vector(x)) {
    return(calculate_nz_mean(x))

    } else if (is.data.frame(as.data.frame(x))) {

    x <- as.data.frame(x)
    if (!is.null(by)) {

      if (is.vector(by)) {
        by <- list(by)
      }

      x <- stats::aggregate(t(x), by, calculate_nz_mean)
      rownames(x) <- x$Group.1
      x <- t(x[, -1])
      return(x)

    } else {
      return(calculate_nz_mean(x))
    }
  }
}

#' Non-zero median of data frames, matrices or vectors
#'
#' Calculate non zero median of a data frame, matrix or vector. If the parameter \code{by} is
#' applied, median will be calculated per groups of columns, not rows.
#'
#' @param x a numerical vector, matrix or data frame.
#' @param by a vector or list of grouping elements, each as long as the variables in the data frame
#' \code{x}. The elements are coerced to factors before use.
#' @return Non zero median of \code{x} or non zero median of groups of columns in \code{x}.
#' @examples
#' vector_test <- c(rep(0,10),1:100)
#' nz_median(vector_test)
#' matrix_test <- matrix(c(rep(0,10),1:100),nrow=10)
#' nz_median(matrix_test, by = c(1,2,3,1,2,3,1,2,3,1,2))
#' df_test <- as.data.frame(matrix_test)
#' nz_median(df_test, by = c(1,2,3,1,2,3,1,2,3,1,2))
#'
#' @export
nz_median <- function(x, by = NULL) {
  # Function to calculate the non zero mean of a vector or data frame/matrix.
  # In case of a data frame/matrix, it can be the non zero mean can be calculated by
  # groups of columns using a vector that specify which column belongs to which group
  # x = data frame, matrix or vector

  calculate_nz_median <- function(vector) {
    # Actual function to calculate the non zero mean
    return(stats::median(vector[vector != 0]))
  }

  if (is.vector(x)) {
    return(calculate_nz_median(x))

  } else if (is.data.frame(as.data.frame(x))) {

    x <- as.data.frame(x)
    if (!is.null(by)) {

      if (is.vector(by)) {
        by <- list(by)
      }

      x <- stats::aggregate(t(x), by, calculate_nz_median)
      rownames(x) <- x$Group.1
      x <- t(x[, -1])
      return(x)

    } else {
      return(calculate_nz_median(x))
    }
  }
}

#' Multiple substitution of values
#'
#' Makes a pair-wise substitions of a vector of values by another vector of values of same length.
#'
#' @param x a numerical vector, matrix or data frame.
#' @param old a vector of values to substitute.
#' @param new a vector of values that will substitute old values.
#' @return a vector with new values.
#' @examples
#' # Replace 1 for 2, 2 for 3 and 3 for 4.
#' msub(1:3, 1:3, 2:4)
#' @export
msub <- function(x, old, new) {
  return(c(new, x)[match(x, c(old, x))])
}


#' Combine list of vectors filling with NA's
#'
#' Combines a list of vectors, filling with NA's to reach the longest vector. Returns a data.frame
#'
#' @param x a list of vectors.
#' @return a data.frame with vectors as columns.
#' @examples
#' # Combine vectors a, b, and c:
#' a <- rnorm(10)
#' b <- rnorm(5)
#' c <- rnorm(9)
#' fill_na(list(a,b,c))
#' @export
fill_na <- function(x) {
  maximum <- max(sapply(x, length))
  add_na <- function(vector, maximum){
    return(c(vector, rep(NA,  maximum - length(vector))))
  } 
  filled_na <- sapply(x, add_na, maximum)
  return(as.data.frame(filled_na))
}
