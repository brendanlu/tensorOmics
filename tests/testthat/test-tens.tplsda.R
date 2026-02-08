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

context("tplsda")

#' Unnames and ensures the top row in a matrix-type output is all positive. This
#' allows for comparison between pls results that are the same but just differ
#' by a negative sign due to svd solver or other implementation detail.
.make_signs_consistent <- function(mat) {
  mat <- unname(mat)
  ncols <- dim(mat)[2]
  for (i in seq_len(ncols)) {
    if (mat[1, i] < 0) {
      mat[, i] <- -mat[, i]
    }
  }
  return(mat)
}

test_that(
  "tplsda agrees with mixOmics plsda",
  code = {
    n <- 6
    p <- 7
    q <- 9
    t <- 1
    k <- min(n, p, q)
    ncomp_input <- 4

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- letters[1:n]

    mixomics_plsda <- mixOmics::plsda(
      test_x[, , 1],
      test_y,
      ncomp = ncomp_input,
      scale = FALSE
    )

    # bltodo: we do not even need to specify identity transforms here right?
    transforms <- matrix_to_m_transforms(diag(1))
    tensor_plsda <- tplsda(
      test_x,
      test_y,
      ncomp = ncomp_input,
      m = transforms$m,
      minv = transforms$minv
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$loadings$X),
      .make_signs_consistent(tensor_plsda$x_loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$loadings$Y),
      .make_signs_consistent(tensor_plsda$y_loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$variates$X),
      .make_signs_consistent(tensor_plsda$x_projected)
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$variates$Y),
      .make_signs_consistent(tensor_plsda$y_projected)
    )
  }
)

test_that(
  "tplsda single time point y imputation is consistent",
  code = {
    n <- 6
    p <- 7
    t <- 3

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y_vec <- letters[1:n]
    test_y_mat <- array(test_y_vec, dim = c(n, 1))
    test_y_tens <- array(test_y_vec, dim = c(n, 1, 1))

    tplsda_vec <- tplsda(test_x, test_y_vec)
    tplsda_mat <- tplsda(test_x, test_y_mat)
    tplsda_tens <- tplsda(test_x, test_y_tens)

    expect_equal(tplsda_vec$y, tplsda_mat$y)
    expect_equal(tplsda_mat$y, tplsda_tens$y)
  }
)

test_that(
  "tplsda multi_label y transformation is appropriate",
  code = {
    # there's only 26 letters in the alphabet, stored in `letters` variable,
    # so this test will automatically fail if n * t > 26
    n <- 6
    p <- 7
    t <- 3

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y_mat <- array(letters[1:(n * t)], dim = c(n, t))
    test_y_tens <- array(letters[1:(n * t)], dim = c(n, 1, t))

    tplsda_mat <- tplsda(test_x, test_y_mat, multi_label = TRUE)
    tplsda_tens <- tplsda(test_x, test_y_tens, multi_label = TRUE)

    # both multi_label y inputs should lead to the same underlying tpls call
    expect_equal(tplsda_mat$y, tplsda_tens$y)

    # the y tensor that goes into the tpls call should have the appropriate
    # number of columns corresponding to all unique categories found across
    # all timepoints specified in y
    expect_equal(dim(tplsda_mat$y), c(n, n * t, t))
  }
)

test_that(
  "tplsda propagates names appropriately",
  code = {
    n <- 3
    p <- 8
    t <- 4

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y_vec <- letters[1:n]

    dimnames(test_x) <- list(
      paste0("r", seq_len(n)),
      paste0("c_x", seq_len(p)),
      paste0("t", seq_len(t))
    )

    tplsda = tplsda(test_x, test_y_vec)

    # rows of loadings output named after the features of the original data
    expect_equal(
      dimnames(tplsda$x_loadings),
      list(paste0("c_x", seq_len(p)), NULL)
    )

    # rows of the variates (projections) output named after the samples of the
    # original data
    expect_equal(
      dimnames(tplsda$x_projected),
      list(paste0("r", seq_len(n)), NULL)
    )
  }
)
