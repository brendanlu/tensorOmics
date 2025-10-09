# ==============================================================================
# Tensor block pls generalization; developed @ Melbourne Integrative Genomics
# Authors:
#   Brendan Lu
#   Saritha Kodikara
# Based on Tenehaus' RGCCA and Kilmer's tensor m-product
# ==============================================================================

#' Return the function implementing the scoring method specified by the user
#' input string; this influences the face on which the algorithm applies the
#' 2D block pls iter function on.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.score_parser <- function(pairwise_scoring) {
  supported_pairwise_scoring <- c(
    "top_singular",
    "sum_singular",
    "frobenius"
  )
  if (!(pairwise_scoring %in% supported_pairwise_scoring)) {
    stop(paste(
      "Please ensure facewise_scoring is one of:",
      paste(supported_pairwise_scoring, collapse = ", ")
    ))
  }

  return(switch(
    pairwise_scoring,
    "top_singular" = .top_singular_score,
    "sum_singular" = .sum_singular_score,
    "frobenius" = .frobenius_score
  ))
}

#' Scores relationships based on the top singular value in the decomposition of
#' a' b, where a and b are the matrices.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.top_singular_score <- function(a, b) {
  svd_decomposition_atb <- svd(crossprod(a, b))
  return(max(svd_decomposition_atb$d))
}

#' Scores relationships based on the sum of the singular values in the
#' decomposition of a' b, where a and b are the matrices.
#'
#' @author Brendadn Lu
#' @keywords internal
#' @noRd
.sum_singular_score <- function(a, b) {
  svd_decomposition_atb <- svd(crossprod(a, b))
  return(sum(svd_decomposition_atb$d))
}

#' Scores relationships based on the frobenius norm approximating (upper bound)
#' top singular value (as computed in `top_singular_score`)
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.frobenius_score <- function(a, b) {
  return(norm(crossprod(a, b), type = "F"))
}

#' Validate params relevant only to the single iterations of a RGCCA
#' optimization, i.e. the parts from Tenehaus.
#'
#' NOTE: design input is validated separately as we use it for more than just
#' the single RGCCA iterations
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.validate_rgcca_params <- function(
  scheme, # to validate
  tau, # to validate
  a, # needed to compute valid tau
  j_total # needed to compute valid tau
) {
  # modes are from Tenehaus page 267, and mainly dictate slight modifications
  # to the objective function
  # provided tau is not changed to "optimal", i.e. tau_j is 1 for all j, and
  # design is "full", then the following hold...
  # "horst": SUMCOV
  # "factorial": SSQCOV
  # "centroid": SABSCOV
  supported_schemes <- c("horst", "factorial", "centroid")
  if (!(scheme %in% supported_schemes)) {
    stop(paste(
      "Please ensure scheme is one of:",
      paste(supported_schemes, collapse = ", ")
    ))
  }

  if (is.null(tau)) {
    # Tenehaus page 267, see SUMCOV, SSQCOV, SABSCOV
    tau <- rep(1, j_total)
  } else if (identical(tau, "optimal")) {
    # tau_hat_star_j from page 263 of Tenehaus, "Schafer and Strimmer formula"
    #
    # empirical estimation of optimal tau (shrinkage constant) from the data to
    # provide better estimate of covariance ("shrinkage estimate of sigma_jj")
    tau <- sapply(a, tau.estimate)
  } else if (is.numeric(tau) && length(tau) == 1) {
    if (tau >= 0 && tau <= 1) {
      tau <- rep(tau, j_total)
    } else {
      stop("For a numeric input, please ensure that 0 <= tau <= 1")
    }
  } else if (is.numeric(tau) && length(tau) >= 1) {
    if (length(tau) != j_total) {
      stop(paste0("Invalid tau input of length ", length(tau), " which needs to
      be length ", j_total, " instead to match the total number of your inputted
      blocks (counting Y as well)"))
    }
  } else {
    stop("Invalid tau input, please ensure it's either NULL, 'optimal', a
    scalar value 0 <= tau <= 1, or a vector of numeric values matching the
    total number of data blocks")
  }

  return(list(
    scheme = scheme,
    tau = tau
  ))
}

#' Dispatch to the appropriate svd solver to obtain top singular vector for
#' depending on size of data. This makes it much more transparent in the
#' calling state whether we are using right or left singular vectors.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.svd_solver <- function(mat, n, p) {
  if (n > 3 && p > 3) {
    return(svds(mat, k = 1))
  } else {
    return(svd(mat, nu = 1, nv = 1))
  }
}

#' Slightly modified and documented single iteration step of rgcca algorithms
#' for 2D array (matrix) blocks.
#'
#' See REGULARIZED GENERALIZED CANONICAL CORRELATION ANALYSIS, Psychometrika
#' vol 76, no 2, pgs 257-284 by Arthur and Michel Tenehaus
#' DOI: 10.1007/S11336-011-9206-8
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.matrix_block_pls_single_iter <- function(
  a, # list of matrices, not tensors
  y_ind,
  design,
  scheme,
  tau,
  nblocks, # note length(a) = nblocks + 1
  n,
  ps, # list of p's corresponding to the dimension p for each block in a
  q,
  tol,
  max_iter,
  solve_dual
) {
  # comments to match with notation in Tenehaus paper for reference
  # see algorithm on page 266
  # NOTE: `a` is capital X in Tenehaus
  # J from Tenehaus, note J = length(a)
  j_total <- nblocks + 1

  # merge all the p's together locally in this function for convenience
  # the Tenehaus method makes no distinction of a "Y" or "X" block, whereas
  # in the main calling function this matters because it affects how we deflate
  ps <- append(ps, q, after = y_ind - 1)

  # just some basic checks so we have valid `scheme` and `tau` to use below
  validated_params <- .validate_rgcca_params(scheme, tau, a, j_total)
  scheme <- validated_params$scheme
  tau <- validated_params$tau

  # copying the approach taken by block.pls, we partition `a` into 'tall' and
  # 'wide' matrices based on their dimensions
  # for the former, the primal problem on the original data X is solved (which
  # is identical to the optimization problem presented in Tenehaus)
  # for the latter, the dual formulation using Gram matrix XXt is solved instead
  if (solve_dual) {
    # more obs than features:
    primal_idxs <- which(ps <= n)
    # more features than obs, "high-dimensional data":
    dual_idxs <- which(ps > n)
    # dual formulation works with Gram matrices
    # be explicit about size here, to indicate we only intend to populate
    # `length(dual_idxs)` entries in the below lists in reality R will pad other
    # indices so the length will just be `j_total` probably
    gram_mats <- vector("list", length = length(dual_idxs))
    gram_loadings <- vector("list", length = length(dual_idxs))
  } else {
    # otherwise just solve the primal for all idxs
    primal_idxs <- seq_len(j_total)
    dual_idxs <- seq_len(0)
  }

  # lowercase a_j's from Tenehaus, termed "normalized outer weight vectors"
  # note must be a list as these are p_j x 1 loading vectors with p_j coming
  # from ps[[j]], so they are of different sizes
  loadings <- vector("list", length = j_total)

  # X_k %*% a_k in (B.) of Tenehaus, interpreted in PLS/PCA-like way as
  # projections of data in loadings space
  # called 'variates' in mixOmics
  projected <- array(0, dim = c(n, j_total))

  # shrinkage estimate of the covariance matrices from Tenehaus, see pg. 263
  # defined by S_hat_j = tau_j * I + (1-tau_j) S_j, where S_j is the normal
  # sample covariance matrix
  # these are all the terms on pg. 266 which are inverted
  s_hat_inv <- vector("list", length = j_total)

  # Part (A.) of Tenehaus: initialization --------------------------------------
  # NOTE: in this loop, we also calculate the inverted shrinkage estimates of
  # the covariance matrices, as they do not change between iterations
  #
  # in terms of initialization, we are copying "svd.single" approach (seems to
  # be the only approach anyway) used in sparse.rgcca_iteration() for
  # 2D block.tpls() function in mixOmics, but this is not necessary and just
  # produces potentially more efficient startings points for a-tildes in A1.
  # of Tenehaus
  for (j in primal_idxs) {
    # inverse of shrinkage estimate of covariance matrix
    s_hat_inv[[j]] <- ginv(
      tau[j] * diag(ps[[j]]) + (1 - tau[j]) * cov2(a[[j]])
    )
    # (A1.) 'sensible' starting point for a-tilde
    init <- .svd_solver(a[[j]], n, ps[[j]])$v
    # (A2.) normalized outer weight vector
    # store the s_hat_inv %*% a_j-tilde term as it appears twice
    temp <- s_hat_inv[[j]] %*% init # this is kind of like the direction
    loadings[[j]] <- as.numeric(1 / sqrt(t(init) %*% temp)) * temp
    # calculate projection
    projected[, j] <- a[[j]] %*% loadings[[j]]
  }

  for (j in dual_idxs) {
    # BLTODO: would be real nice have a man page displaying the maths here
    # calculate Gram matrix over samples (XXt), NOT features (XtX)
    gram_mats[[j]] <- a[[j]] %*% t(a[[j]])
    # for dual, this is more like a covariance across samples, we temporarily
    # need the non-inverted s_hat to calculate the initial normalized
    # "outer weight vector" on XXt, `gram_loadings`
    s_hat <- tau[j] * diag(n) + ((1 - tau[j]) / (n - 1)) * gram_mats[[j]]
    s_hat_inv[[j]] <- ginv(s_hat)
    # (A1.) 'sensible' starting point for a-tilde for XXt instead of X
    init <- .svd_solver(a[[j]], n, ps[[j]])$u
    # (A2.) normalized outer weight vector for XXt instead of X
    gram_loadings[[j]] <- as.numeric(
      1 / sqrt(t(init) %*% s_hat %*% gram_mats[[j]] %*% init)
    ) * init
    # premultiply by Xt to get a p x 1 loadings matrix in original space
    loadings[[j]] <- t(a[[j]]) %*% gram_loadings[[j]]
    # calculate projection back into original space
    # should be more efficient to project gram matrix onto gram loadings,
    # but projecting original data onto loadings with a[[j]] %*% loadings[[j]]
    # should be fine too
    projected[, j] <- gram_mats[[j]] %*% gram_loadings[[j]]
  }
  #-----------------------------------------------------------------------------

  # iteration ------------------------------------------------------------------
  iters <- 0
  # preallocate array to store inner component each iteration
  inner_components <- array(0, dim = c(n, j_total))

  while (TRUE) {
    loadings_prev_iter <- loadings # to compare for convergence
    for (j in seq_len(j_total)) {
      # BLTODO: the current mixOmics implementation does not seem to iterate in
      # ascending order of j. from inspection of Tenehaus algorithm, this seems
      # like it would affect the calculation of the inner component, it might be
      # wrong?

      # Part (B.) of Tenehaus: compute inner component - this is not affected by
      # primal/dual formulation

      # inner term w(Cov(Xjaj, Xkak)) changes depending on what w is
      # store as 1 x j_total array s.t. w_term[k] = cov(Xjaj, Xkak)
      w_term <- switch(
        # bltodo: we could consider a switch that selects the function w()
        # instead, but the benefit of the below is that for the default "horst"
        # scheme we currently save an unecessary computation of the covariance
        scheme,
        "horst" = array(1, dim = c(1, j_total)),
        "factorial" = cov2(projected, projected[, j]),
        "centroid" = sign(cov2(projected, projected[, j]))
      )
      # z_j in Tenehaus, n x 1 vector stored in column of `inner_components`
      inner_components[, j] <- rowSums(sweep(
        # the sweep returns n x j_total array where each column contains a term
        # in the summation of the inner product, therefore taking row sums
        # takes the sum of these columns, i.e. the result of whole summation
        projected, # n x j_total array where each column contains Xkak
        2, # scales column k of projected by cjk * w(Cov(Xjaj, Xkak))
        design[j, ] * w_term, # 1 x j_total array, cjk * w(Cov(Xjaj, Xkak))
        "*"
      ))

      # Part (C.) of Tenehaus: compute weight vector - different approach
      # depending on if we are solving the primal or dual formulation
      if (j %in% primal_idxs) {
        # t(X_j) %*% z_j term, like the a-tilde term in (A2.)
        # kind of like temporarily projecting data onto updated inner component
        proj_to_inner <- t(a[[j]]) %*% inner_components[, j]
        # s_hat_inv %*% t(X_j) %*% z_j term, similar to the temp term in (A2.)
        # kind of gives us the s_hat_inv corrected direction of new loading
        temp <- s_hat_inv[[j]] %*% proj_to_inner
        # calculate the whole normalized outer weight vector
        loadings[[j]] <- as.numeric(1 / sqrt(t(proj_to_inner) %*% temp)) * temp
        # update projections, i.e. Xkak terms in (B.) to iteration s + 1
        projected[, j] <- a[[j]] %*% loadings[[j]]
      }

      if (j %in% dual_idxs) {
        # in the dual formulation, iteration step updates the gram loadings
        temp <- s_hat_inv[[j]] %*% inner_components[, j]
        gram_loadings[[j]] <- as.numeric(
          1 / sqrt(t(inner_components[, j]) %*% gram_mats[[j]] %*% temp)
        ) * temp
        # premultiply by Xt to get a p x 1 loadings matrix in original space
        loadings[[j]] <- t(a[[j]]) %*% gram_loadings[[j]]
        # calculate projection back into original space
        # should be more efficient to project gram matrix onto gram loadings,
        # but projecting original data onto loadings with
        # a[[j]] %*% loadings[[j]] should be fine too
        projected[, j] <- gram_mats[[j]] %*% gram_loadings[[j]]
      }
    }

    iters <- iters + 1
    # BLTODO: think this through - is this a reasonable measure of convergence?
    # what happens if some loadings are just bigger vectors because of ps[[j]]
    # being massive, should we normalize at all by vec length?
    # should it be monotonically decreasing?
    max_loadings_change <- max(sapply(
      seq_len(j_total),
      function(j) {
        crossprod(loadings[[j]] - loadings_prev_iter[[j]])
      }
    ))
    if (max_loadings_change < tol) {
      break
    }
    if (iters >= max_iter) {
      stop("The RGCCA algorithm did not converge, please try increasing
      `max_iter`, increasing `tol`, or toggling `solve_dual`")
    }
  }

  return(list(
    loadings = loadings,
    projected = projected,
    iters = iters
  ))
  #-----------------------------------------------------------------------------
}

#' Boolean: are all elements of the list the same?
#' OPTIONAL: are they all the same as the specific value `compare_with`?
#' Helper for validating tpls inputs.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.all_same <- function(l, compare_with = NULL) {
  if (is.null(compare_with)) {
    return(all(sapply(l, function(x) x == l[[1]])))
  } else {
    return(all(sapply(l, function(x) x == compare_with)))
  }
}

#' Run some sense checks on the block inputs. Note it splits the inputs into
#' "X" and "Y" blocks to produce some more informative error messages. Imagine
#' in research context user has probably gone to different datasets to merge
#' into one input list, so we give unique error messages if specifically the
#' target block has incompatible dimensions to the other blocks.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.validate_block_tpls_x_y <- function(x_blocks, y) {
  # BLTODO: improve these error messages to be more informative in the y_ind
  # input case which is currently being supported first
  if (!is.list(x_blocks)) {
    stop("x must be a list of longitudinal datasets")
  }
  nblocks <- length(x_blocks)

  # check that all inputs in the x list are order 3 arrays (tensors)
  orders <- lapply(x_blocks, function(tens) length(dim(tens)))
  if (!.all_same(orders, 3)) {
    stop("Please ensure all blocks are order-3 arrays")
  }

  # check the compatible dimensions across the list of tensors
  ns <- lapply(x_blocks, function(tens) dim(tens)[1])
  ps <- lapply(x_blocks, function(tens) dim(tens)[2])
  ts <- lapply(x_blocks, function(tens) dim(tens)[3])

  if (!(.all_same(ns) && .all_same(ts))) {
    stop("Please ensure all blocks have the same number of samples and
    timepoints")
  }

  # check compatibility with y input
  if (length(dim(y)) != 3) {
    stop("Please ensure target input tensor is an order-3 array")
  } else {
    n <- dim(y)[1]
    q <- dim(y)[2]
    t <- dim(y)[3]
  }

  if (n != ns[[1]]) {
    stop("Please ensure target input tensor matches the number of samples found
    within the other blocks")
  }

  if (t != ts[[1]]) {
    stop("Please ensure target input tensor matches the number of time points
    found within the other blocks")
  }

  return(list(
    nblocks = nblocks,
    n = n,
    ps = ps,
    q = q,
    t = t
  ))
}

#' Check the design matrix input ensuring it is valid, or parse keyword inputs.
#' See .create_design() in utils.R for reference. Rewriting this new function
#' was easier than refactoring that one to work well in block_tpls.
#'
#' Apologies in advance for the maintenence burden this might cause...
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.stop_invalid_design_input <- function(design, nblocks, y_ind) {
  # RECALL: nblocks is the number of "x" blocks, meaning we have nblocks + 1
  # if we include our response tensor y
  # this also serves as a flag as to whether or not we create the design matrix
  off_diag <- NULL

  # parse keyword input, or numeric input
  if (identical(design, "full")) {
    off_diag <- 1
  } else if (identical(design, "null")) {
    off_diag <- 0
  } else if (is.numeric(design) && length(design) == 1) {
    # is design scalar?, R must have better way to check than the crap above...?
    if (design > 0 && design < 1) {
      off_diag <- design
    } else {
      stop("For a numeric input, please ensure that 0 < design < 1")
    }
  }

  if (!is.null(off_diag)) {
    # create and return design matrix
    design <- matrix(off_diag, nrow = nblocks + 1, ncol = nblocks + 1)
    design[, y_ind] <- 1
    design[y_ind, ] <- 1
    diag(design) <- 0
    return(design)
  } else {
    # run a whole bunch of checks to make sure user inputted design matrix
    # makes sense
    # IF YOU CHANGE THESE CHECKS BELOW, YOU MUST CONSIDER THE IMPACT ON THE
    # SCORING FUNCTION (AND ADJUST IF NEEDED), AS IT ASSUMES A SYMMETRIC
    # 0-DIAGONAL STRUCTURE IN THE DESIGN MATRIX
    if (ncol(design) != (nblocks + 1) || nrow(design) != nblocks + 1) {
      stop(paste0(
        "Based on the number of inputted tensor blocks, `design` must be a
        square matrix with ",
        nblocks + 1,
        " columns and rows"
      ))
    }
    if (!isSymmetric(design)) {
      stop("Custom inputted design matrix must be symmetric, i.e. relationship
      between block Ai and Aj is the same as relationship between block
      Aj and Ai")
    }
    if (!(all(diag(design) == 0))) {
      # bltodo: actually could remove this for some more freedom
      stop("Diagonals of the design matrix must be 0's, i.e. do not aim to 
      model block relationships with themself - otherwise consider tpca method")
    }
    return(design)
  }
}

#' Deflate a matrix based on a vector of projections.
#'
#' @author Brendan Lu
#' @keywords internal
#' @noRd
.deflate <- function(mat, projected) {
  # see tpls.R code for details to match literature
  reg_coef <- crossprod(mat, projected) / as.numeric(crossprod(projected))
  return(mat - tcrossprod(projected, reg_coef))
}

#' Tensor block PLS
#'
#' Developed @ Melbourne Integrative Genomics.
#'
#' @param a A list of tensor inputs (termed 'blocks') measured on the same
#' samples.
#' @param y_ind Index of "target" block in `a`.
#' @param ncomp The estimated number of components. ncomp must be explicitly set
#' as an integer in tpls.
#' @param design Numeric matrix of size length(x) x length(x) with values
#' between 0 and 1; the i,j entry in the matrix indicates the strength of the
#' relationship to be modelled between the i-th and j-th blocks. A value of 0
#' indicates no relationship, with 1 being the maximum. Alternatively, one can
#' input "null" for a fully disconnected design (feature blocks are only
#' connected to the target block at strength 1, but not to each other), or
#' "full" for a fully connected design (target and feature blocks are all
#' connected to each other equally at strength 1), or a single scalar value
#' between 0 and 1 which will designate the relationship between feature blocks,
#' with the relationships to the target block being 1. "full" by default.
#' @param scheme One of "horst", "factorial" or "centroid", "horst" by default.
#' @param tau Shrinkage constant to provide better estimate of covariance,
#' introduced by Ledoit and Wolf (2004). NULL by default, but can set to
#' "optimal" to use optimal shrinkage constants from Schafer and Strimmer
#' (2005), or a single scalar value between 0 and 1.
#' @param m A function which applies an orthogonal tensor tubal transform.
#' @param minv The inverse of m.
#' @param mode Currently supports tensor analogues of canonical, regression,
#' and svd PLS modes. Defaults to "regression" mode.
#' @param pairwise_scoring The method on which to select the face of each tensor
#' for each iteration of block.pls. Currently supports "top_singular" (pick
#' based on weighted sum of top svd diagonal value within all 2-block
#' combinatations), "sum_singular" (sum of diagonal values) and
#' "frobenius" (top_singular but approximated with the Frobenius norm
#' which is an upper bound). "top_singular" by default, but consider selecting
#' "frobenius" for speed as it avoids a full svd decomposition.
#' @param center If set to FALSE, the data tensor will not be centralized into
#' Mean Deviation Form (see Mor et al. 2022). By default, the mean horizontal
#' slice of the input tensor(s) are subtracted, so that all of the horizontal
#' slices sum to 0, analgous to centering matrix data.
#' @param tol Positive scalar used as convergence criteria/tolerance during the
#' RGCCA iterative process. Defaults to \code{1e-06}.
#' @param max_iter Integer specifying the maximum number of iterations to be run
#' for each 2D RGCCA call.
#' @param solve_dual For high dimensional data where n > p, use a dual
#' formulation in the matrix RGCCA step that works with the gram matrix instead
#' of the raw data. TRUE by default, only recommended to set to FALSE for
#' testing purposes.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation. Does not have any effect if transform functions
#' explicitly set using \code{m}, \code{minv}.
#' @note When `design = "full"`, `tau = 1`, setting scheme to "horst",
#' "factorial" or "centroid" yields method SUMCOV, SSQCOV, SABSCOV respectively.
#'
#' @author Brendan Lu
#' @export
block_tpls <- function(
  a,
  y_ind,
  ncomp = 1,
  design = "full",
  scheme = "horst",
  tau = 1,
  m = NULL,
  minv = NULL,
  mode = "regression",
  pairwise_scoring = "top_singular",
  center = TRUE,
  tol = 1e-06,
  max_iter = 100,
  solve_dual = TRUE,
  bpparam = NULL
) {
  # bltodo: add docs link to existing block.pls for user reference
  # bltodo: need to validate y_ind input, and also support y input later

  # currently supports canonical and regression modes, which affects how we
  # deflate the data if ncomp > 1
  supported_modes <- c("canonical", "regression")
  if (!(mode %in% supported_modes)) {
    stop(paste(
      "Please ensure mode is one of:",
      paste(supported_modes, collapse = ", ")
    ))
  }

  # do some rough checks and return relevant dimensions
  dims <- .validate_block_tpls_x_y(a[-y_ind], a[[y_ind]])
  nblocks <- dims$nblocks # PLEASE NOTE: this is the number of "x" blocks
  n <- dims$n
  ps <- dims$ps # probably a bit confusing but this is only for "x" blocks
  q <- dims$q # this is nfeatures for the target block
  t <- dims$t
  min_k <- min(n, unlist(ps), q)
  max_rank <- min_k * t

  # check user ncomp is sensible, or impute NULL input
  ncomp <- .validate_integer_ncomp(ncomp, max_rank)

  # use dctii as default transform if user does not specify an explicit one
  validated_transforms <- .stop_invalid_transform_input(m, minv, t, bpparam)
  m <- validated_transforms$m
  minv <- validated_transforms$minv

  # check / parse the design input
  validated_design <- .stop_invalid_design_input(design, nblocks, y_ind)

  # bltodo: `facewise_scoring` parameter needs a better docstring and probably
  # some documentation in a man file as this is pretty new stuff. also,
  # parameter name might be a bit crap...?
  # bltodo: it would be cool to let the user pass in any binary function that
  # matches the correct signature to implement their own score, but implementing
  # this probably creates more headache than realistic user utility
  pairwise_scoring_func <- .score_parser(pairwise_scoring)

  if (center) {
    # almost identical to centering in tpca / tpls but applied on a list
    mean_slice_a <- lapply(a, FUN = function(tens) apply(tens, c(2, 3), mean))
    a <- Map(
      function(tens, mean_slice) sweep(tens, c(2, 3), mean_slice, "-"),
      a, mean_slice_a
    )
  }

  # transform into hat-space
  # BLTODO: this is no longer centered?!
  ahat <- lapply(a, m)

  # block_combination_idxs$i_idx looks like 1, 2, 3, 1, 2, 3, ...
  # block_combination_idxs$j_idx looks like 1, 1, 1, 2, 2, 2, ...
  #
  # if we apply a bivariate transformation over these indicies
  #   i.e. (1, 1), (2, 1), (3, 1), (1, 2), ...
  # we are effectively applying on the indices of a matrix specified in Fortran
  # style, so if we reshape the output of our transformation via matrix()
  # in R (which fills the matrix in column-major order) we will get the same
  # type of output as order(), such that element Aij = F(i, j)
  block_combination_idxs <- expand.grid(
    i_idx = seq_len(nblocks + 1),
    j_idx = seq_len(nblocks + 1)
  )

  # things to store for each component iteration
  tensor_blocks_loadings <- lapply(
    append(ps, q, after = y_ind - 1),
    function(nfeatures) array(0, dim = c(nfeatures, ncomp))
  )
  tensor_blocks_projected <- lapply(
    seq_len(nblocks + 1),
    function(i) array(0, dim = c(n, ncomp))
  )
  names(tensor_blocks_loadings) <- names(tensor_blocks_projected) <- names(a)
  iters <- array(0, dim = ncomp)
  faces <- array(0, dim = ncomp)

  # for each component, find a "maximum" tensor face based on the pairwise
  # scoring method and then perform 2D RGCCA methods to obtain loadings and
  # projections for each block
  for (i in seq_len(ncomp)) {
    # pick which face to perform 2D block tpls iteration on
    #
    # alright, the below is a big chunk of nested function calls, but it's the
    # only way found so far to avoid clunky R for-loops
    #
    # so, the outer lapply obtains 'scores' for each face across t's
    # the function basically computes, for each t ("tt"):
    # a "weighted" pairwise_scores matrix = design %*% t(pairwise_scores),
    # where pairwise_scores is a matrix such that element i,j contains the
    # pairwise score between matrices a[[i]][, , tt] and a[[j]][, , tt], and
    # pre-multiplying effectively weights these values according to design
    #
    # it then computes a single score by taking the total sum of all elements in
    # the "weighted" (by design) pairwise_scores matrix
    facewise_scores <- lapply(seq_len(t),
      function(tt) {
        sum(
          validated_design %*% t(
            matrix(
              # the conversion to matrix works perfectly with the output from
              # block_combination_idxs as it's filled in column major order
              # and block_combination_idxs gives the matrix indices in this
              # order as well
              mapply( # bltodo: the 'outer' function does not work somehow?
                # THE BELOW ASSUMES THAT THE `validated_design` MATRIX IS A
                # SYMMETRIC 0-DIAGONAL MATRIX, AND `pairwise_scoring_func` IS A
                # COMMUTATIVE FUNCTION SO IT ONLY CALCULATES THE LOWER TRIANGLE
                # OF PAIRWISE SCORES
                function(i_idx, j_idx) {
                  if (i_idx < j_idx) {
                    return(
                      pairwise_scoring_func(
                        ahat[[i_idx]][, , tt],
                        ahat[[j_idx]][, , tt]
                      )
                    )
                  } else {
                    return(0)
                  }
                },
                block_combination_idxs$i_idx, block_combination_idxs$j_idx
              ),
              nrow = nblocks + 1, ncol = nblocks + 1
            )
          )
        )
      }
    )

    # now we obtain the index of the face we will perform 2D block pls iteration
    t_curr <- which.max(unlist(facewise_scores))
    faces[i] <- t_curr

    # pick out the matrices for each tensor block corresponding to t_curr face
    mats_curr <- lapply(
      seq_len(nblocks + 1),
      function(j) ahat[[j]][, , t_curr]
    )

    # obtain loadings and projections that maximise a network covariance
    # function introduced by Tenehaus, between the 2D matrices
    matrix_block_pls_single_iter <- .matrix_block_pls_single_iter(
      a = mats_curr,
      y_ind = y_ind,
      design = validated_design,
      scheme = scheme,
      tau = tau,
      nblocks = nblocks,
      n = n,
      ps = ps,
      q = q,
      tol = tol,
      max_iter = max_iter,
      solve_dual = solve_dual
    )

    # list of p_j x 1 loading vectors
    curr_all_mats_loadings <- matrix_block_pls_single_iter$loadings
    # n x j_total matrix where each column contains the projections onto the
    # feature loadings for each block
    curr_all_mats_projected <- matrix_block_pls_single_iter$projected

    # store iters for user info
    iters[i] <- matrix_block_pls_single_iter$iters

    # perform deflation --------------------------------------------------------
    # just modify the `t_curr` faces of each tensor in-place instead of
    # rewriting elements inside `ahat` completely
    for (j in seq_len(nblocks + 1)) {
      if (j == y_ind) {
        if (mode == "canonical") {
          ahat[[j]][, , t_curr] <- .deflate(
            mats_curr[[j]],
            curr_all_mats_projected[, j]
          )
        }
        if (mode == "regression") {
          # deflates the "Y" matrix against all other projected values, and then
          # takes the mean
          ahat[[j]][, , t_curr] <- Reduce(
            "+",
            lapply(
              seq_len(nblocks + 1)[-y_ind],
              # calculate "Y" matrix deflated against all other projections
              function(jj) {
                return(.deflate(
                  mats_curr[[j]],
                  curr_all_mats_projected[, jj]
                ))
              }
            )
          ) / nblocks
        }
      } else {
        ahat[[j]][, , t_curr] <- .deflate(
          mats_curr[[j]],
          curr_all_mats_projected[, j]
        )
      }

      # also fill in the loadings and projections in this for-loop
      tensor_blocks_loadings[[j]][, i] <- curr_all_mats_loadings[[j]]
      tensor_blocks_projected[[j]][, i] <- curr_all_mats_projected[, j]
    }
  }

  # BLTODO: need to decide if it is worth splitting X and Y blocks in the output
  # what will be more convenient for user...?
  return(list(
    ncomp = ncomp,
    a = a,
    design = design,
    scheme = scheme,
    mode = mode,
    pairwise_scoring = pairwise_scoring,
    loadings = tensor_blocks_loadings,
    projected = tensor_blocks_projected,
    iters = iters,
    faces = faces
  ))
}
