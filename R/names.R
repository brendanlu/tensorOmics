# ==============================================================================
# Standardized utils for propogating names into outputs of tpca, tpls,
# block_tpls based methods
# ==============================================================================

#' @author Brendan Lu
#' @keywords internal
#' @noRd
.dimnames_for_loadings_mat <- function(in_tens) {
  # the rows of this output 2D array correspond to the features of the original
  # tensor input
  if (!is.null(dimnames(in_tens))) {
    return(list(
      dimnames(in_tens)[[2]],
      NULL
    ))
  } else {
    return(NULL)
  }
}

#' @author Brendan Lu
#' @keywords internal
#' @noRd
.dimnames_for_projections_mat <- function(in_tens) {
  # the rows of this output 2D array correspond to the original rows (samples)
  # of the tensor input
  if (!is.null(dimnames(in_tens))) {
    return(list(
      dimnames(in_tens)[[1]],
      NULL
    ))
  } else {
    return(NULL)
  }
}

#' @author Brendan Lu
#' @keywords internal
#' @noRd
.dimnames_for_loadings_tens <- function(in_tens) {
  # the rows of this output 2D array correspond to the features of the original
  # tensor input
  if (!is.null(dimnames(in_tens))) {
    return(list(
      dimnames(in_tens)[[2]],
      NULL,
      dimnames(in_tens)[[3]]
    ))
  } else {
    return(NULL)
  }
}

#' @author Brendan Lu
#' @keywords internal
#' @noRd
.dimnames_for_projections_tens <- function(in_tens) {
  # the rows of this output 2D array correspond to the original rows (samples)
  # of the tensor input
  if (!is.null(dimnames(in_tens))) {
    return(list(
      dimnames(in_tens)[[1]],
      NULL,
      dimnames(in_tens)[[3]]
    ))
  } else {
    return(NULL)
  }
}
