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
      order(as.vector(s_mat))[1:ncomp],
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
#' @param y Tensor input.
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
#' @param matrix_output Note: ONLY AFFECTS THE OUTPUT IN "tsvdm" MODE.
#' Collects the top singular vectors across the tensor, organized by the
#' magnitude of the corresponding singular vector, and place them into the
#' columns of a matrix. Corresponds to the 'matrix compression' type of scheme
#' described in Mor et al. 2022.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation. Does not have any effect if transform functions
#' explicitly set using \code{m}, \code{minv}.
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

  dims <- .validate_tpls_x_y_dim(x, y)
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

  # check user ncomp is sensible, or impute NULL input
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
  yhat <- m(y)

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
      # bltodo: tpls explained variance?
      k_t_flatten_sort <- .obtain_k_t_flatten_sort(
        tsvdm_decomposition_xt_y$s,
        ncomp
      )
      x_loadings <- .extract_tensor_columns(x_loadings, k_t_flatten_sort)
      y_loadings <- .extract_tensor_columns(y_loadings, k_t_flatten_sort)
      x_projected <- .extract_tensor_columns(x_projected, k_t_flatten_sort)
      y_projected <- .extract_tensor_columns(y_projected, k_t_flatten_sort)
    }

    return(invisible(list(
      ncomp = ncomp,
      x = x,
      y = y,
      x_loadings = x_loadings,
      y_loadings = y_loadings,
      x_projected = x_projected,
      y_projected = y_projected
    )))

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
      curr_y_projected <- yhat[, , k_t_top[2]] %*% curr_y_loadings

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
      # ------------------------------------------------------------------------

      x_loadings[, i] <- curr_x_loadings
      y_loadings[, i] <- curr_y_loadings
      x_projected[, i] <- curr_x_projected
      y_projected[, i] <- curr_y_projected
    }

    return(invisible(list(
      ncomp = ncomp,
      x = x,
      y = y,
      mode = mode,
      x_loadings = x_loadings,
      y_loadings = y_loadings,
      x_projected = x_projected,
      y_projected = y_projected,
      features = features,
      faces = faces
    )))

  } else {
    stop("Unexpected error in tpls, check 'mode' parameter input")
  }
}
