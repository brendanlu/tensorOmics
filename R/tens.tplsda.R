# tensorOmics - Temporal Omics Data Analysis Package
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ==============================================================================
# Tensor plsda generalization; developed @ Melbourne Integrative Genomics based
# on Kilmer's tensor m-product algebra
# ==============================================================================

#' Check if y input is a vector, matrix or tensor that can be inferred as fixed
#' class labels for each sample across all time points; if so return these
#' classes in consistent form as a vector.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.y_to_vec <- function(y) {
  # input is already a n x 1 vector
  if (is.null(dim(y))) {
    return(y)
  }
  # input is a n x 1 matrix
  if (length(dim(y)) == 2 && dim(y)[2] == 1) {
    return(c(y))
  }
  # input is a n x 1 x 1 tensor
  if (length(dim(y)) == 3 && dim(y)[2] == 1 && dim(y)[3] == 1) {
    return(c(y))
  }
  stop("
    Please ensure y input is either a vector, or a matrix / tensor with one
    vertical column representing outcome categories for each sample

    If specifying categories across time points, please specify
    `multi_label = TRUE` in function call and pass in a matrix or tensor with
    the number of columns corresponding to the number of time points
  ")
}

#' Check if y input is a matrix or tensor that can be inferred as repeated class
#' measurements for each sample; if so, return these classes in a consistent
#' form as a n x 1 x t tensor.
#'
#' Note: We do not have to check that the dimensions are compatible with x
#' tensor input, this check will be done in the call to tpls.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.y_to_tens <- function(y) {
  # input is a n x t matrix
  if (length(dim(y)) == 2) {
    n <- dim(y)[1]
    t <- dim(y)[2]
    return(array(y, dim = c(n, 1, t)))
  }
  # input is a n x 1 x t tensor
  if (length(dim(y)) == 3 && dim(y)[2] == 1) {
    return(y)
  }
  stop("
    When `multi_label = TRUE` in tplsda, please ensure that y input is either a
    n x t matrix, or n x 1 x t column tensor, representing the classes for
    each sample across the t repeated measurements.
  ")
}

#' Preprocess a categorical "Y" inputs tensor for discriminant analysis. Returns
#' in a tensor form.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.create_indicator_tensor <- function(y, multi_label, t_from_x) {
  if (multi_label) {
    y_tens <- .y_to_tens(y)
    n <- dim(y_tens)[1]
    t <- dim(y_tens)[3]

    # we need to determine all the unique levels that show up across all
    # repeated observations for all samples
    factors_list <- lapply(array(seq_len(t)), function(i) factor(y_tens[, , i]))
    unique_levels <- unique(unlist(lapply(factors_list, levels)))
    if (length(unique_levels) == 1) {
      stop("y should contain 2 or more distinct classes")
    }

    y <- array(0, dim = c(n, length(unique_levels), t))
    for (i in seq_len(t)) {
      # unmap conveniently allows us to manually control the total sets of
      # groups (class measurements) our input sample is drawn from
      y[, , i] <- unmap(factors_list[[i]], unique_levels)
      # bltodo: preserve names, see below
    }
  } else {
    y_vec <- .y_to_vec(y)
    y_vec <- factor(y_vec)
    if (nlevels(y_vec) == 1) {
      stop("y should contain 2 or more distinct classes")
    }
    y_mat <- unmap(y_vec)
    # BLTODO: WANT colnames(y_mat) = unique_levels, and be able to propogate
    # names into tpls - need to work on preserving names for tpls and tpca
    # when mixOmics2 is rolled out to end users (see existing plsda)

    # fill y matrix into a tensor with the same t as x tensor input
    y <- array(y_mat, dim = c(dim(y_mat), t_from_x))
  }

  return(y)
}

#' Run tensor PLS-DA analysis
#'
#' Note: always returns a compressed-matrix form output.
#'
#' Developed @ Melbourne Integrative Genomics
#'
#' @inheritParams tpls
#'
#' @param y A vector / column matrix / column tensor with n class labels, or a
#' n x t matrix / n x 1 x t tensor if `multi_label == TRUE`.
#' @param multi_label Set to `TRUE` if y contains repeated class measurements
#' across the t timepoints specified in the x tensor.
#'
#' @inherit tpls return
#' @return Adds \code{y_original}: the original labels input
#' (vector/matrix/tensor) passed to \code{tplsda()}.
#'
#' @author Brendan Lu
#' @export
tplsda <- function(
  x,
  y,
  multi_label = FALSE,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  center = TRUE,
  bpparam = NULL
) {
  # most parameter checking is done in tpls call, only the checks that are
  # unique / necessary for tplsda are done here
  if (length(dim(x)) != 3) {
    stop("Please ensure x input tensor is an order-3 array")
  } else {
    t <- dim(x)[3]
  }

  # make a copy of the original y
  y_original_copy <- y

  # note the rest of parameter validation will be done in tpls call
  output <- tpls(
    x = x,
    y = .create_indicator_tensor(y, multi_label, t),
    ncomp = ncomp,
    m = m,
    minv = minv,
    mode = "regression",
    center = center,
    bpparam = bpparam
  )
  output$y_original <- y_original_copy
  class(output) <- "tplsda"
  return(invisible(output))
}
