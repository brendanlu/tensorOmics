# ==============================================================================
# Tensor block pls generalization; developed @ Melbourne Integrative Genomics
# Authors:
#   Brendan Lu
#   Saritha Kodikara
# Based on Tenehaus' RGCCA and Kilmer's tensor m-product
# ==============================================================================

#' Tensor block PLSDA-like analysis.
#'
#' Developed @ Melbourne Integrative Genomics.
#'
#' @inheritParams block_tpls
#'
#' @inherit block_tpls note
#'
#' @param x A list of tensor inputs (termed 'blocks') measured on the same
#' samples that describe the output.
#' @param y A vector / column matrix / column tensor with n class labels, or a
#' n x t matrix / n x 1 x t tensor if `multi_label == TRUE`.
#' @param multi_label Set to TRUE if y contains repeated class measurements
#' across the t timepoints specified in the x tensor.
#'
#' @author Brendan Lu
#' @export
block_tplsda <- function(
  x,
  y,
  multi_label = FALSE,
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
  # just get a provisional value for t here, proper parameter checking is done
  # in the call to block_tpls
  if (length(dim(x[[1]])) != 3) {
    stop("Please ensure all blocks are order-3 arrays")
  } else {
    t <- dim(x[[1]])[3]
  }

  # create and add the "Y" target tensor into the end of the list
  a <- x
  y_ind <- length(x) + 1
  a[[y_ind]] <- .create_indicator_tensor(y, multi_label, t)

  # if the user looks like they've passed a design matrix for the "X" blocks,
  # be nice and append a "Y" column for them
  #
  # this is probably unlikely, and any dodgy inputs will causes error messages
  # in the block_tpls call anyway
  if (
    length(dim(design)) == 2 && # passed a matrix-looking input
      ncol(design) == length(x) && # ncols matches the length of their `x` input
      nrow(design) == length(x) # nrows matches the length of their `x` input
  ) {
    design <- cbind(design, 1)
    design <- rbind(design, 1)
    # we will not overwrite the other diagonal values, if the user inputted it
    # wrong it will cause an informative error message in block_tpls
    design[y_ind, y_ind] <- 0
  }

  # note the rest of parameter validation will be done in the block_tpls call
  return(invisible(block_tpls(
    a = a,
    y_ind = y_ind,
    ncomp = ncomp,
    design = design,
    scheme = scheme,
    tau = tau,
    m = m,
    minv = minv,
    mode = mode,
    pairwise_scoring = pairwise_scoring,
    center = center,
    tol = tol,
    max_iter = max_iter,
    solve_dual = solve_dual,
    bpparam = bpparam
  )))
}
