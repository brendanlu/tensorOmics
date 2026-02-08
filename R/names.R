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
