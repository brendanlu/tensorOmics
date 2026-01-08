# ==============================================================================
# utilities and functions for Kilmer's m product
# ==============================================================================
# bltodo: use classes for m and minv?

#' Apply a function across the last dimension of an input vector,
#' matrix, or tensor. This function defines both a parallel algorithm using
#' BiocParrallel, and also a simple \code{apply} algorithm. Even on Windows
#' Machines, setting \code{bpparam = BiocParallel::SerialParam()} offers a
#' notable speedup for larger 3D array inputs.
#'
#' @param x Numerical array input.
#' @param mat Function which defines the tubal transform.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} object
#' indicating the type of parallelisation.
#' @return A tensor of the same size under the specified tubal transform,
#' denoted \eqn{\hat{x}}.
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.apply_mat_transform <- function(x, mat, bpparam) {
  # first cache names for x, using names if a vector otherwise dimnames for
  # array type inputs
  # search: ".apply_mat_transform based mode-3 product propogates names"
  cached_names <- if (is.vector(x)) names(x) else dimnames(x)

  if (length(dim(x)) == 1) {
    out <- mat %*% x
    names(out) <- cached_names
  } else if (length(dim(x)) == 2) {
    out <- mat %*% x
    dimnames(out) <- cached_names
  } else if (length(dim(x)) == 3) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    t <- dim(x)[3]
    if (is.null(bpparam)) {
      # apply algorithm --------------------------------------------------------
      transformed_slices <- apply(
        aperm(x, c(3, 1, 2)), c(2, 3), # permuted tensor with mode 3 in front
        function(slice) mat %*% slice # left multiply by transform mat
      )
      # permute back to original orientation
      out <- aperm(array(transformed_slices, dim = c(t, n, p)), c(2, 3, 1))
      # ------------------------------------------------------------------------
    } else {
      # BiocParallel algorithm -------------------------------------------------
      transformed_slices <- BiocParallel::bplapply(
        lapply(seq_len(p), function(i) t(x[, i, ])), # a list of t x n matrices
        function(slice) mat %*% slice, # left multiply by transform mat
        BPPARAM = bpparam
      )
      # unlist and cast into array with p facewise t x n matrices
      # then appropriately rotate to get original n x p x t matrix
      out <-
        aperm(array(unlist(transformed_slices), dim = c(t, n, p)), c(2, 3, 1))
      # ------------------------------------------------------------------------
    }
    dimnames(out) <- cached_names
  } else {
    # error: some array >3D has been inputted
    stop(
      "Only order 1 (vector), 2 (matrix), 3 (3D tensor) arrays are supported"
    )
  }
  return(out)
}

#' Validate appropriate 'null-ness' of m, minv inputs, and apply default
#' transform as needed.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.stop_invalid_transform_input <- function(m, minv, t, bpparam) {
  if (
    xor(is.function(m), is.function(minv)) ||
      xor(is.null(m), is.null(minv))
  ) {
    stop(
      "If explicitly defined, both m and its inverse must be defined as 
      functions"
    )
  }

  # due to above checks, the code below executes when both m and minv are NULL
  # and returns the calling state the default transform used in this library
  # which is currently the dct-ii transform
  if (is.null(m)) {
    transforms <- dctii_m_transforms(t, bpparam = bpparam)
    m <- transforms$m
    minv <- transforms$minv
  }

  return(list(
    m = m,
    minv = minv
  ))
}

#' Create m-transform functions as defined by a matrix
#'
#' Returns functions \code{m} and \code{m_inv} which apply tubal transforms
#' defined by the matrix \code{m_mat}.
#'
#' @param m_mat Matrix which defines the tubal transform.
#' @param m_inv_mat Matrix which defines inverse tubal transform
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} object
#' indicating the type of parallelisation.
#' @return
#' \item{m}{A function which applies the matrix m_mat along the last dimension
#' of a given numerical input array. For 3D tensors it performs the mode-3
#' product.}
#' \item{minv}{The inverse of m.}
#' @author Brendan Lu
#' @export
matrix_to_m_transforms <- function(
  m_mat,
  m_inv_mat = NULL,
  bpparam = NULL
) {
  # error: non-matrix input
  stopifnot(length(dim(m_mat)) == 2)
  # error: non-square matrix input
  stopifnot(dim(m_mat)[1] == dim(m_mat)[2])

  # invert m_mat if m_inv_mat not specified
  if (is.null(m_inv_mat)) {
    # error: non-invertible input
    m_inv_mat <- solve(m_mat)
  } else {
    # error: specified matrices are not the same size
    stopifnot(identical(dim(m_mat), dim(m_inv_mat)))
  }

  return(list(
    m = function(x) .apply_mat_transform(x, m_mat, bpparam = bpparam),
    minv = function(x) .apply_mat_transform(x, m_inv_mat, bpparam = bpparam)
  ))
}

#' Create m-transform functions to apply the Discrete Cosine Transform 2
#'
#' Returns functions \code{m} and \code{m_inv} which apply tubal transforms
#' defined by the Discrete Cosine Transform (DCT-II variant). This is equivalent
#' to Scipys DCTI-ii algorithm with \code{norm='ortho'}.
#'
#' @param t The length of the transform.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} object
#' indicating the type of parallelisation.
#' @return
#' \item{m}{A function which applies the dct-ii along the last dimension of a
#' given numerical input array. For 3D tensors it performs the mode-3 product
#' with the DCT matrix.}
#' \item{minv}{The inverse of m,}
#' @author Brendan Lu
#' @export
dctii_m_transforms <- function(t, bpparam = NULL) {
  # handle edge case when t = 1, because gsignal::dctmtx(1) returns a scalar
  # note that mathematically dctmtx(n = 1) is indeed the identity
  if (t == 1) {
    m_mat <- diag(1)
  } else {
    m_mat <- gsignal::dctmtx(t)
  }
  return(matrix_to_m_transforms(m_mat = m_mat, bpparam = bpparam))
}

#' For two name vectors, return a named vector which "summarises" them by:
#' - if both name vectors are NULL, return NULL
#' - if one name vector is NULL, return the other
#' - if both name vectors are not fully NULL, return a named vector which has
#' non-NULL elements only if the names are identical in both inputs, OR one of
#' the inputs has NULL at that position (where the other has a non-NULL value
#' at that position which we can use)
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.resolve_names <- function(names1, names2) {
  stopifnot(length(names1) == length(names2))

  # both NULL case handled appropriately from below
  if (is.null(names1)) return(names2)
  if (is.null(names2)) return(names1)

  # both non-NULL case, logic tested in unit tests
  # search: "facewise product resolves names correctly, and propogates them"
  out <- rep(NA, length(names1))
  match_mask <- names1 == names2 & !is.na(names1) & !is.na(names2)
  out[match_mask] <- names1[match_mask]
  out[is.na(names1)] <- names2[is.na(names1)]
  out[is.na(names2)] <- names1[is.na(names2)]
  return(out)
}

#' Compute Kilmer's facewise product. Note that the for-loop implementation is
#' relatively fast, and very readable. There's also a BiocParralel
#' implementation here, but it lacks significant benchmarking results.
#' bltodo: the parallel algorithm is probably stupid remove sometime?
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.binary_facewise <- function(a, b, bpparam) {
  # require 3d array inputs
  stopifnot(length(dim(a)) == 3 && length(dim(b)) == 3)

  na <- dim(a)[1]
  pa <- dim(a)[2]
  ta <- dim(a)[3]

  nb <- dim(b)[1]
  pb <- dim(b)[2]
  tb <- dim(b)[3]

  # error: different t for each input
  stopifnot(ta == tb)
  t <- ta
  # error: non-conforming facewise dimensions
  stopifnot(pa == nb)

  if (!is.null(dimnames(a)) || !is.null(dimnames(b))) {
    # attempt to cache some sensible names
    out_dimnames <- list(
      dimnames(a)[[1]],
      dimnames(b)[[2]],
      .resolve_names(dimnames(a)[[3]], dimnames(b)[[3]])
    )
  } else {
    # otherwise no point having a list of NULL's
    out_dimnames <- NULL
  }

  if (is.null(bpparam)) {
    # for-loop algorithm -------------------------------------------------------
    out <- array(0, dim = c(na, pb, t))
    for (i in seq_len(t)) {
      out[, , i] <- a[, , i] %*% b[, , i]
    }
    # --------------------------------------------------------------------------
  } else {
    # BiocParallel algorithm ---------------------------------------------------
    # bltodo: benchmark / investigate preallocation here?
    out <- simplify2array(
      BiocParallel::bplapply(
        array(seq_len(t)),
        FUN = function(i) a[, , i] %*% b[, , i],
        BPPARAM = bpparam
      )
    )
    # --------------------------------------------------------------------------
  }
  dimnames(out) <- out_dimnames
  return(out)
}

#' Tensor facewise product
#'
#' Compute Kilmer's facewise product cumulatively across any arbitrary number of
#' tensor inputs.
#'
#' @param ... Arbitrary number of numerical tensor inputs.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} object
#' indicating the type of parallelisation.
#' @return Cumulative facewise product.
#' @author Brendan Lu
#' @export
facewise_product <- function(..., bpparam = NULL) {
  return(
    Reduce(
      function(a, b) .binary_facewise(a, b, bpparam = bpparam),
      list(...)
    )
  )
}

#' @describeIn facewise_product Custom facewise product operator
#' @param a Tensor input.
#' @param b Tensor input.
#' @export
`%fp%` <- function(a, b) .binary_facewise(a, b, bpparam = NULL)

#' Facewise transpose an order-3 tensor
#'
#' @param tensor Numerical 3D array input.
#' @return Facewise transpose of \code{tensor}
#' @author Brendan Lu
#' @aliases ft facewise_transpose
#' @export
facewise_transpose <- function(tensor) {
  return(aperm(tensor, c(2, 1, 3)))
}

#' @describeIn facewise_transpose Alias for \code{\link{facewise_product}}
#' @export
ft <- function(tensor) facewise_transpose(tensor)

#' Kilmer's tensor-tensor m-product
#'
#' This function works cumulatively across any arbitrary number of tensor
#' inputs.
#'
#' @param ... Arbitrary number of numerical tensor inputs.
#' @param m A function which applies an orthogonal tensor tubal transform.
#' @param minv The inverse of m.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} object
#' indicating the type of parallelisation. Does not have any effect if transform
#' functions explicitly set using \code{m}, \code{minv}.
#' @return Cumulative m-product.
#' @author Brendan Lu
#' @export
m_product <- function(
  ...,
  m = NULL,
  minv = NULL,
  bpparam = NULL
) {
  tensors <- list(...)
  if (length(tensors) == 0) {
    stop("No input tensors provided.")
  } else {
    t <- dim(tensors[[1]])[3]
  }
  if (
    xor(is.function(m), is.function(minv)) ||
      xor(is.null(m), is.null(minv))
  ) {
    stop(
      "If explicitly defined, both m and its inverse must be defined as 
      functions."
    )
  }
  # use dctii as default transform if user does not specify an explicit one
  if (is.null(m)) {
    transforms <- dctii_m_transforms(t, bpparam = bpparam)
    m <- transforms$m
    minv <- transforms$minv
  }
  return(minv(
    Reduce(
      # bpparam MUST be NULL here to prevent double parallelisation!
      function(a, b) .binary_facewise(a, b, bpparam = NULL),
      lapply(list(...), m)
    )
  ))
}

#' Tensor cross product
#'
#' Compute the equivalent of ft(a) %fp% b
#'
#' @param a Tensor input.
#' @param b Tensor input.
#' @return Tensor facewise cross product.
#' @author Brendan Lu
#' @export
facewise_crossproduct <- function(a, b) {
  na <- dim(a)[1]
  pa <- dim(a)[2]
  ta <- dim(a)[3]

  nb <- dim(b)[1]
  pb <- dim(b)[2]
  tb <- dim(b)[3]

  # error: different t for each input
  stopifnot(ta == tb)
  t <- ta
  # error: non-conforming facewise dimensions
  stopifnot(na == nb)

  fcp_ab <- array(0, dim = c(pa, pb, t))
  for (i in seq_len(t)) {
    fcp_ab[, , i] <- crossprod(a[, , i], b[, , i])
  }
  return(fcp_ab)
}
