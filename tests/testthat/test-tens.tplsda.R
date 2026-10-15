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

test_that(
  "tplsda compares every face with a single label per subject",
  code = {
    n <- 20
    p <- 10
    t <- 4

    set.seed(1)
    # features 1 to 5 rise over time in group "a" and fall in group "b": the
    # same mean, a different shape
    test_y <- rep(c("a", "b"), each = n / 2)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    for (j in 1:5) {
      test_x[, j, ] <- test_x[, j, ] +
        outer(ifelse(test_y == "a", 1, -1), c(-1.5, -0.5, 0.5, 1.5))
    }

    # a transformed label is non-zero on face 1 only, the time averages
    expect_equal(
      as.vector(tplsda(test_x, test_y, ncomp = 2, transform_y = TRUE)$faces),
      c(1, 1)
    )

    tplsda_obj <- tplsda(test_x, test_y, ncomp = 2)
    expect_false(tplsda_obj$transform_y)
    expect_true(tplsda_obj$faces[1] != 1)
  }
)

test_that(
  "tplsda scales faces by default for a single label per subject",
  code = {
    n <- 30
    p <- 20
    t <- 4

    set.seed(1)
    # a difference in shape only (features 1 to 5), plus subject offsets that
    # are constant over time: they sit on face 1 and make it the largest face
    test_y <- rep(c("a", "b"), each = n / 2)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    for (j in 1:5) {
      test_x[, j, ] <- test_x[, j, ] +
        outer(ifelse(test_y == "a", 0.5, -0.5), c(-1.5, -0.5, 0.5, 1.5))
    }
    test_x <- test_x + array(rnorm(n * p, sd = 2), dim = c(n, p, t))

    tplsda_unscaled <- tplsda(test_x, test_y, ncomp = 1, face_scaling = "none")
    expect_equal(as.vector(tplsda_unscaled$faces), 1)

    tplsda_obj <- tplsda(test_x, test_y, ncomp = 1)
    expect_equal(tplsda_obj$face_scaling, "frobenius")
    expect_equal(as.vector(tplsda_obj$faces), 2)

    # labels repeated over time are transformed, and their faces not scaled
    tplsda_multi <- tplsda(
      test_x, matrix(test_y, n, t), multi_label = TRUE, ncomp = 1
    )
    expect_equal(tplsda_multi$face_scaling, "none")
  }
)

test_that(
  "tplsda returns the class names in the order of the indicator columns",
  code = {
    n <- 18
    p <- 10
    t <- 3

    set.seed(1)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    # labels whose sorted order is not their order of appearance
    test_y <- rep(c("b", "10", "a"), each = n / 3)

    tplsda_obj <- tplsda(test_x, test_y, ncomp = 2)
    expect_false(tplsda_obj$multi_label)
    expect_equal(tplsda_obj$levels, c("10", "a", "b"))
    # the centered indicator column of each class is positive exactly for its
    # samples
    for (j in seq_along(tplsda_obj$levels)) {
      expect_equal(tplsda_obj$y[, j, 1] > 0, test_y == tplsda_obj$levels[j])
    }

    # labels repeated over time: the levels of each time point in turn, so
    # "0", which appears only at the last one, comes last although it sorts
    # first
    test_y_time <- matrix(test_y, n, t)
    test_y_time[test_y == "a", t] <- "0"
    tplsda_multi <- tplsda(test_x, test_y_time, multi_label = TRUE, ncomp = 2)
    expect_true(tplsda_multi$multi_label)
    expect_equal(tplsda_multi$levels, c("10", "a", "b", "0"))
  }
)
