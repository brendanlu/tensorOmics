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
# Tensor pls generalization; developed @ Melbourne Integrative Genomics based
# on Kilmer's tensor m-product algebra
# ==============================================================================

#' Convert a singular values tensor (in compressed matrix form) to a set of
#' indices corresponding to the (column,face) pairs of the top `ncomp` singular
#' values. NEEDS singular values to be in matrix form.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.obtain_k_t_flatten_sort <- function(s_mat, ncomp) {
  # bltodo: if ultimately only use once just place it in function body directly
  return(
    .unravel_index(
      order(as.vector(s_mat), decreasing = TRUE)[1:ncomp],
      dim(s_mat)
    )
  )
}

#' Get the k, t index corresponding to the largest singular value. Mirrors
#' .unravel_index() but more efficient as we only care about the largest
#' singular value in tpls.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.obtain_k_t_top <- function(s_mat) {
  flat_index <- which.max(as.vector(s_mat))
  nrows <- dim(s_mat)[1]
  return(c(
    (flat_index - 1) %% nrows + 1, # transformed k tensor position
    (flat_index - 1) %/% nrows + 1 # transformed t tensor position
  ))
}

#' Run some sense checks on x and y tensor inputs for tpls.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.validate_tpls_x_y_dim <- function(x, y) {
  if (length(dim(x)) != 3) {
    stop("Please ensure x input tensor is an order-3 array")
  } else {
    n <- dim(x)[1]
    p <- dim(x)[2]
    t <- dim(x)[3]
  }

  if (length(dim(y)) != 3) {
    stop("Please ensure y input tensor is an order-3 array")
  } else {
    n2 <- dim(y)[1]
    q <- dim(y)[2]
    t2 <- dim(y)[3]
  }

  if (n != n2) {
    stop("Please ensure x and y tensor inputs have matching number of samples")
  }

  if (t != t2) {
    stop("Please ensure x and y tensor inputs have matching number of time 
    points")
  }

  return(list(
    n = n,
    p = p,
    q = q,
    t = t
  ))
}

#' Lay a response with no time dimension (a vector or an n x q matrix), which
#' is fixed per subject, out over t time points as a tensor constant in time.
#'
#' @keywords internal
#' @noRd
.subject_level_to_tensor <- function(y, t) {
  y_mat <- as.matrix(y)
  y_dimnames <- dimnames(y_mat)
  if (!is.null(y_dimnames)) {
    y_dimnames <- c(y_dimnames, list(NULL))
  }
  return(array(y_mat, dim = c(dim(y_mat), t), dimnames = y_dimnames))
}

#' Stop if a response that is to be left untransformed changes over time: only
#' a response fixed per subject has a single matrix to meet every face.
#'
#' @keywords internal
#' @noRd
.stop_time_varying_untransformed <- function(y) {
  if (any(y != as.vector(y[, , 1, drop = FALSE]))) {
    stop(paste(
      "transform_y = FALSE needs a response that is fixed per subject: a",
      "vector, an n x q matrix, or a tensor that is the same at every time point"
    ))
  }
}

#' Check integer ncomp inputs. Used for tensor PLS-like analyses. Depends on a
#' sensible calculation of max_rank which is left to the calling state.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.validate_integer_ncomp <- function(ncomp, max_rank) {
  # process ncomp input, much simpler than tpca as we only accept integer input
  # bltodo: investigate non integer input options? explained variance semantics?
  if (is.null(ncomp)) {
    ncomp <- max_rank
  } else if (ncomp %% 1 == 0) {
    stopifnot(ncomp > 0 && ncomp <= max_rank)
  } else {
    stop("Please input an integer or NULL for ncomp parameter")
  }

  return(ncomp)
}

#' Tensor PLS-like regression
#'
#' Developed @ Melbourne Integrative Genomics. Intuition: Fourier-based m
#' transforms compact the data across time points, but preserve the information
#' association for each sample and each feature. In doing so, pls works by
#' calculating the singular vectors associated with a 'global' (across time
#' points) top singular value, and deflating just that relevant face.
#'
#' @param x Tensor input.
#' @param y Tensor input, or a vector / n x q matrix for a response fixed per
#' subject (see `transform_y`).
#' @param ncomp The estimated number of components. ncomp must be explicitly set
#' as an integer in tpls.
#' @param m A function which applies an orthogonal tensor tubal transform.
#' @param minv The inverse of m.
#' @param mode Currently supports tensor analogues of canonical ("canonical"),
#' regression ("regression"), and svd ("tsvdm") PLS variants. Defaults to
#' "regression".
#' @param center If set to false, the data tensor will not be centralized into
#' Mean Deviation Form (see Mor et al. 2022). By default, the mean horizontal
#' slice of the input tensor(s) are subtracted, so that all of the horizontal
#' slices sum to 0, analgous to centering matrix data.
#' @param matrix_output Note: FALSE ONLY AFFECTS THE OUTPUT IN "tsvdm" MODE.
#' TRUE by default. Collects the top singular vectors across the tensor,
#' organized by the magnitude of the corresponding singular vector, and place
#' them into the columns of a matrix. Corresponds to the 'matrix compression'
#' type of scheme described in Mor et al. 2022.
#' @param transform_y Whether to m-transform the response with x. NULL by
#' default, which decides from the shape of y: a vector or an n x q matrix is a
#' response fixed per subject, such as an endpoint measured once, and is not
#' transformed, so that every face of x is compared with the response itself;
#' an order-3 tensor is transformed. A response that is constant over time and
#' transformed is non-zero on the first face only, which under the default DCT
#' holds the subject time averages. An untransformed response is deflated on
#' every face after each component. TRUE copies a subject-fixed response over
#' time and transforms it; FALSE needs a response that is the same at every
#' time point.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} object
#' indicating the type of parallelisation. Does not have any effect if transform
#' functions explicitly set using \code{m}, \code{minv}.
#' @return A list containing:
#' \describe{
#'   \item{ncomp}{The number of components used.}
#'   \item{x}{The centered input tensor X.}
#'   \item{y}{The centered input tensor Y.}
#'   \item{mode}{The PLS mode used: "canonical", "regression", or "tsvdm".}
#'   \item{transform_y}{Whether the response was m-transformed.}
#'   \item{x_loadings}{Matrix of X loadings (p x ncomp), or loadings tensor
#'   if `mode = "tsvdm"` AND `matrix_output = FALSE`.}
#'   \item{y_loadings}{Matrix of Y loadings (q x ncomp), or loadings tensor
#'   if `mode = "tsvdm"` AND `matrix_output = FALSE`.}
#'   \item{x_projected}{Matrix of X projections (n x ncomp), or proections
#'   tensor when `mode = "tsvdm"` AND `matrix_output = FALSE`.}
#'   \item{y_projected}{Matrix of Y projections (n x ncomp), or proections
#'   tensor when `mode = "tsvdm"` AND `matrix_output = FALSE`.}
#'   \item{features}{Vector of selected feature indices per component (only
#'   when `matrix_output = TRUE`).}
#'   \item{faces}{Vector of selected tensor face indices per component (only
#'   when `matrix_output = TRUE`).}
#' }
#' @author Brendan Lu
#' @export
tpls <- function(
  x,
  y,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  mode = "regression",
  center = TRUE,
  matrix_output = TRUE, # FALSE only takes effect for tsvdm method
  transform_y = NULL,
  bpparam = NULL
) {
  # BLTODO: add option to minv back into original space...? currently it does
  # nothing

  # allowed modes mirror sklearn's 2D PLS models
  # see https://scikit-learn.org/stable/modules/cross_decomposition.html
  supported_modes <- c("canonical", "regression", "tsvdm")
  if (!(mode %in% supported_modes)) {
    stop(paste(
      "Please ensure mode is one of:",
      paste(supported_modes, collapse = ", ")
    ))
  }

  # a response with no time dimension is fixed per subject, and is not
  # transformed unless asked; lay it out over the time points of x
  subject_level_y <- length(dim(y)) < 3
  if (is.null(transform_y)) {
    transform_y <- !subject_level_y
  }
  if (subject_level_y && length(dim(x)) == 3) {
    y <- .subject_level_to_tensor(y, dim(x)[3])
  }

  dims <- .validate_tpls_x_y_dim(x, y)
  if (!transform_y) {
    .stop_time_varying_untransformed(y)
  }
  n <- dims$n
  p <- dims$p
  q <- dims$q
  t <- dims$t
  k <- min(n, p, q)
  # maximum number of non-zero entries in the f-diagonal singular values tensor
  # of tsvdm(XtY) is k * t
  # BLTODO: of course the rank could be smaller than this, is it worth
  # calculating to give user a more informative error if they pass rank
  # deficient data they are not aware of...?
  max_rank <- k * t

  # check user ncomp is sensible, or impute NULL input; remember whether the
  # user asked for a number, to warn if fewer components can be returned
  ncomp_requested <- !is.null(ncomp)
  ncomp <- .validate_integer_ncomp(ncomp, max_rank)

  # use dctii as default transform if user does not specify an explicit one
  validated_transforms <- .stop_invalid_transform_input(m, minv, t, bpparam)
  m <- validated_transforms$m
  minv <- validated_transforms$minv

  if (center) {
    mean_slice_x <- apply(x, c(2, 3), mean)
    mean_slice_y <- apply(y, c(2, 3), mean)
    x <- sweep(x, c(2, 3), STATS = mean_slice_x, FUN = "-")
    y <- sweep(y, c(2, 3), STATS = mean_slice_y, FUN = "-")
  }

  # project x and y into 'hat-space', expensive computation so do it once and
  # save into variable
  xhat <- m(x)
  # an untransformed response stays as it is, so every face of xhat meets the
  # response itself
  yhat <- if (transform_y) m(y) else y

  # store the columns and faces of the top singular vector in each component
  # iteration
  features <- array(NA, dim = ncomp)
  faces <- array(NA, dim = ncomp)

  if (mode == "tsvdm") {
    # simplest algorithm - just uses everything from the tsvdm call of XtY based
    # without any deflation steps
    tsvdm_decomposition_xt_y <- tsvdm(
      # bltodo: change to facewise_crossproduct when this is improved
      ft(xhat) %fp% yhat,
      transform = FALSE,
      svals_matrix_form = TRUE
    )

    # loadings are just the u,v tensors from the tsvdm call
    x_loadings <- tsvdm_decomposition_xt_y$u
    y_loadings <- tsvdm_decomposition_xt_y$v

    # project the scaled / transformed x, y data onto the loadings
    x_projected <- xhat %fp% x_loadings
    y_projected <- yhat %fp% y_loadings

    if (matrix_output) {
      # singular values below rounding noise (the tolerance the component loop
      # below uses) become NA: order() puts NA last whichever way it sorts, so
      # they are never selected
      svals <- tsvdm_decomposition_xt_y$s
      noise <- svals < sqrt(.Machine$double.eps) * max(svals)
      if (ncomp > sum(!noise)) {
        if (ncomp_requested) {
          warning(sprintf(paste(
            "Only %d of the %d requested components were returned: after",
            "component %d, t(X) %%*%% Y has no structure left to extract"
          ), sum(!noise), ncomp, sum(!noise)))
        }
        ncomp <- sum(!noise)
      }
      svals[noise] <- NA
      # bltodo: tpls explained variance?
      k_t_flatten_sort <- .obtain_k_t_flatten_sort(
        svals,
        ncomp
      )
      # populate features and faces from k_t_flatten_sort
      features <- k_t_flatten_sort[1, ]
      faces <- k_t_flatten_sort[2, ]
      # fill return matrices
      x_loadings <- .extract_tensor_columns(x_loadings, k_t_flatten_sort)
      y_loadings <- .extract_tensor_columns(y_loadings, k_t_flatten_sort)
      x_projected <- .extract_tensor_columns(x_projected, k_t_flatten_sort)
      y_projected <- .extract_tensor_columns(y_projected, k_t_flatten_sort)
    }

  } else if (mode == "canonical" || mode == "regression") {
    # there is lots of literature on this, e.g. Kim-Anh's mixOmics book, or
    # Wegelin: "A Survey of Partial Least Squares Methods with Emphasis on the
    # Two-Block Case" (2000), which one can see all describe the same thing
    # albeit with different notation
    #
    # also see: https://scikit-learn.org/stable/modules/cross_decomposition.html
    # and Kim-Anh's book to understand canonical vs regression, as the
    # algorithm in Wegelin page 10 describes canonical only
    #
    # comments added to cross match code with Wegelin (pg10), Kim-Anh (pg169)
    # and scikit-learn (PLS codebase) notation for easier review

    # preallocate output arrays which we will fill during the iterative process
    x_loadings <- array(0, dim = c(p, ncomp))
    y_loadings <- array(0, dim = c(q, ncomp))
    x_projected <- array(0, dim = c(n, ncomp))
    y_projected <- array(0, dim = c(n, ncomp))
    ncomp_kept <- ncomp

    for (i in seq_len(ncomp)) {
      # compute tsvdm
      # BLTODO: tensor crossprod to speed up?
      # bltodo: minor as ncomp is usually small, but we technically do not need
      # to recalculate the whole tsvdm, we only need svd re-done on the face
      # that was deflated as everything else is still the same
      tsvdm_decomposition_xt_y <- tsvdm(
        # bltodo: change to facewise_crossproduct when this is improved
        ft(xhat) %fp% yhat,
        transform = FALSE,
        svals_matrix_form = TRUE
      )

      # get indices corresponding to largest singular value
      k_t_top <- .obtain_k_t_top(tsvdm_decomposition_xt_y$s)

      # once t(X) %*% Y has no structure left, its largest singular value is
      # rounding noise and its singular vectors are arbitrary, so stop there
      # (the relative tolerance is the one MASS::ginv uses)
      curr_singular_value <- tsvdm_decomposition_xt_y$s[k_t_top[1], k_t_top[2]]
      if (i == 1) {
        first_singular_value <- curr_singular_value
      } else if (
        curr_singular_value < sqrt(.Machine$double.eps) * first_singular_value
      ) {
        ncomp_kept <- i - 1
        break
      }

      features[i] <- k_t_top[1]
      faces[i] <- k_t_top[2]
      # u, v (Wegelin); a, b (Kim-Anh); weights (scikit-learn)
      # NOTE: due to the svd, these are normalized, but in scikit-learn's
      # regression mode, the curr_y_loadings are obtained using another method
      # and deliberately left unnormalized, whereas Kim-Anh and mixOmics
      # uses normalized y loadings
      curr_x_loadings <- tsvdm_decomposition_xt_y$u[, k_t_top[1], k_t_top[2]]
      curr_y_loadings <- tsvdm_decomposition_xt_y$v[, k_t_top[1], k_t_top[2]]

      # perform deflation ------------------------------------------------------
      # this is deliberately done in a fashion to match scikit-learn's code
      # and we could consider returning the reg_coef's in the future
      # using k_t_top[2] we just work in matrix world for this step

      # xi, omega (Wegelin); t, u (Kim-Anh); scores (scikit-learn)
      # note that only one face of xhat and yhat is relevant per iteration
      curr_x_projected <- xhat[, , k_t_top[2]] %*% curr_x_loadings
      # matrix() keeps the face n x q when q = 1, where `[` would drop it to a
      # vector and %*% would reject the vector-by-vector product
      curr_y_projected <- matrix(yhat[, , k_t_top[2]], nrow = n) %*%
        curr_y_loadings

      # gamma (Wegelin); c (Kim-Anh); loadings (scikit-learn)
      # the calculation of the x regression coefficient and deflation step of
      # xhat is the same regardless of pls mode
      curr_x_reg_coef <- crossprod(xhat[, , k_t_top[2]], curr_x_projected) /
        as.numeric(crossprod(curr_x_projected))

      xhat[, , k_t_top[2]] <- xhat[, , k_t_top[2]] -
        tcrossprod(curr_x_projected, curr_x_reg_coef)

      # the calculation of the y regression coefficient and deflation step of
      # yhat differs depending on the pls mode
      if (mode == "canonical") {
        # delta (Wegelin); e (Kim-Anh); loadings (scikit-learn canonical)
        curr_y_reg_coef <- crossprod(yhat[, , k_t_top[2]], curr_y_projected) /
          as.numeric(crossprod(curr_y_projected))

        yhat[, , k_t_top[2]] <- yhat[, , k_t_top[2]] -
          tcrossprod(curr_y_projected, curr_y_reg_coef)
      }

      if (mode == "regression") {
        # d (Kim-Anh); loadings (scikit-learn regression)
        # Wegelin pg 10 only describes canonical
        curr_y_reg_coef <- crossprod(yhat[, , k_t_top[2]], curr_x_projected) /
          as.numeric(crossprod(curr_x_projected, curr_x_projected))

        yhat[, , k_t_top[2]] <- yhat[, , k_t_top[2]] -
          tcrossprod(curr_x_projected, curr_y_reg_coef)
      }

      # an untransformed response is one matrix shared by every face, so it is
      # deflated on every face
      if (!transform_y) {
        yhat <- array(yhat[, , k_t_top[2]], dim = dim(yhat))
      }
      # ------------------------------------------------------------------------

      x_loadings[, i] <- curr_x_loadings
      y_loadings[, i] <- curr_y_loadings
      x_projected[, i] <- curr_x_projected
      y_projected[, i] <- curr_y_projected
    }

    if (ncomp_kept < ncomp) {
      if (ncomp_requested) {
        warning(sprintf(paste(
          "Only %d of the %d requested components were returned: after",
          "component %d, t(X) %%*%% Y has no structure left to extract"
        ), ncomp_kept, ncomp, ncomp_kept))
      }
      keep <- seq_len(ncomp_kept)
      x_loadings <- x_loadings[, keep, drop = FALSE]
      y_loadings <- y_loadings[, keep, drop = FALSE]
      x_projected <- x_projected[, keep, drop = FALSE]
      y_projected <- y_projected[, keep, drop = FALSE]
      features <- features[keep]
      faces <- faces[keep]
      ncomp <- ncomp_kept
    }

  } else {
    stop("Unexpected error in tpls, check 'mode' parameter input")
  }

  # set dimnames
  if (!matrix_output && mode == "tsvdm") {
    # only way for output to be in tensor form
    dimnames(x_loadings) <- .dimnames_for_loadings_tens(x)
    dimnames(x_projected) <- .dimnames_for_projections_tens(x)
    dimnames(y_loadings) <- .dimnames_for_loadings_tens(y)
    dimnames(y_projected) <- .dimnames_for_projections_tens(y)
  } else {
    dimnames(x_loadings) <- .dimnames_for_loadings_mat(x)
    dimnames(x_projected) <- .dimnames_for_projections_mat(x)
    dimnames(y_loadings) <- .dimnames_for_loadings_mat(y)
    dimnames(y_projected) <- .dimnames_for_projections_mat(y)
  }

  # prepare output
  output <- list(
    ncomp = ncomp,
    x = x,
    y = y,
    mode = mode,
    transform_y = transform_y,
    x_loadings = x_loadings,
    y_loadings = y_loadings,
    x_projected = x_projected,
    y_projected = y_projected,
    features = features,
    faces = faces
  )
  class(output) <- "tpls"
  return(invisible(output))
}
