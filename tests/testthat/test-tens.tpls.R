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

context("tpls")

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
  "tsvdm-tpls mode agrees with tpca",
  code = {
    # it actually takes a lot of coercing to make these two outputs comparable
    # because their default configurations are aimed at completely different
    # problems. that being said, this test here is probably still useful as a
    # sense check.
    #
    # problem 1: squaring the singular values means that the k_t_flatten_sort
    # compression is different and picks out a different compressed matrix. Sort
    # order of the squared values is not necessarily the same.
    #
    # problem 2: loadings (and therefore projections) may differ in sign, think
    # this is dependent on the svd solver implementation?
    #
    # problem 3: if you use centering (on for tpca and tpls by default), the
    # faces of each tensor loses one rank, so only rank - 1 columns of each face
    # in the loadings tensor ("v") will match. for tsvdm mode this means we get
    # one less column to compare.
    #
    # problem 4: tpls computes the svd on XtX, which is rank-deficient if X is
    # non-square. this means that the last few columns of "v" are absolutely
    # crap and should not be compared. these columns may not even exist in
    # tpca depending on the input dimensions.
    n <- 4
    p <- 5
    t <- 3
    k <- min(n, p)

    # this test fails if `test_tensor` is rank deficient, as is the case if we
    # use array(1:24, dim = c(3, 4, 2))
    set.seed(1)
    test_tensor <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    transforms <- dctii_m_transforms(t)
    m <- transforms$m
    minv <- transforms$minv

    tpca_obj <- tpca(
      test_tensor,
      m = m,
      minv = minv,
      # turn off centering to get over problem 3
      center = FALSE,
      # return full tensors to get over problem 1
      matrix_output = FALSE
    )

    tpls_tsvdm <- tpls(
      test_tensor, test_tensor,
      m = m,
      minv = minv,
      mode = "tsvdm",
      # turn off centering to get over problem 3
      center = FALSE,
      # return full tensors to get over problem 1
      matrix_output = FALSE
    )

    # compare elementwise absolute values to get over problem 2
    expect_equal(
      abs(tpca_obj$loadings),
      # only compare k columns to get over problem 4
      abs(tpls_tsvdm$y_loadings[, 1:k, ])
    )

    expect_equal(
      abs(tpca_obj$variates),
      abs(tpls_tsvdm$y_projected)
    )
  }
)

test_that(
  "the first component of tpls regression and canonical match tpca",
  code = {
    n <- 4
    p <- 5
    t <- 3
    k <- min(n, p)

    set.seed(1)
    test_tensor <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    transforms <- dctii_m_transforms(t)
    m <- transforms$m
    minv <- transforms$minv

    tpca_obj <- tpca(
      test_tensor,
      m = m,
      minv = minv
    )

    tpls_regression <- tpls(
      test_tensor, test_tensor,
      m = m,
      minv = minv,
      mode = "regression",
    )

    expect_equal(tpls_regression$y_loadings[, 1], tpca_obj$loadings[, 1])

    # bltodo: i think this just passes on tolerance, they are not necessarily
    # mathematically equal for three columns...so whoever fails this test in the
    # future feel free to consider changing the indices from `1:3` to `1`
    expect_equal(
      .make_signs_consistent(tpls_regression$y_loadings[, 1:3]),
      .make_signs_consistent(tpca_obj$loadings[, 1:3])
    )

    tpls_canonical <- tpls(
      test_tensor, test_tensor,
      m = m,
      minv = minv,
      mode = "canonical",
    )

    expect_equal(tpls_canonical$y_loadings[, 1], tpca_obj$loadings[, 1])
  }
)

test_that(
  "canonical and regression modes: tpls agrees with MixOmics pls",
  code = {
    n <- 5
    p <- 8
    q <- 10
    t <- 1
    ncomp_input <- 4
    modes_to_test <- c("canonical", "regression")

    k <- min(n, p, q)

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- array(rnorm(n * q * t, mean = 0, sd = 3), dim = c(n, q, t))

    for (curr_mode in modes_to_test) {
      mixomics_pls <- mixOmics::pls(
        test_x[, , 1],
        test_y[, , 1],
        ncomp = ncomp_input,
        scale = FALSE,
        mode = curr_mode
      )

      tensor_pls <- tpls(
        test_x,
        test_y,
        ncomp = ncomp_input,
        mode = curr_mode
      )

      expect_equal(
        .make_signs_consistent(mixomics_pls$loadings$X),
        .make_signs_consistent(tensor_pls$x_loadings)
      )

      expect_equal(
        .make_signs_consistent(mixomics_pls$loadings$Y),
        .make_signs_consistent(tensor_pls$y_loadings)
      )

      expect_equal(
        .make_signs_consistent(mixomics_pls$variates$X),
        .make_signs_consistent(tensor_pls$x_projected)
      )

      expect_equal(
        .make_signs_consistent(mixomics_pls$variates$Y),
        .make_signs_consistent(tensor_pls$y_projected)
      )
    }
  }
)

test_that(
  "sense checks on the tpls returned features and faces items",
  code = {
    n <- 5
    p <- 8
    q <- 13
    t <- 4
    ncomp_input <- 3
    modes_to_test <- c("regression", "tsvdm", "regression")

    k <- min(n, p, q)

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- array(rnorm(n * q * t, mean = 0, sd = 3), dim = c(n, q, t))
    k <- min(n, p)

    for (curr_mode in modes_to_test) {
      tensor_pls <- tpls(
        test_x,
        test_y,
        ncomp = ncomp_input,
        mode = curr_mode
        # matrix_output must be TRUE for this test, should be by default
      )
      expect_equal(tensor_pls$ncomp, ncomp_input)
      expect_equal(length(tensor_pls$features), ncomp_input)
      expect_equal(length(tensor_pls$faces), ncomp_input)
      expect_equal(dim(tensor_pls$x_loadings), c(p, ncomp_input))
      expect_equal(dim(tensor_pls$y_loadings), c(q, ncomp_input))
      expect_equal(dim(tensor_pls$x_projected), c(n, ncomp_input))
      expect_equal(dim(tensor_pls$y_projected), c(n, ncomp_input))
    }
  }
)

test_that(
  "tpls error for non-supported mode input",
  code = {
    n <- 4
    p <- 5
    t <- 3
    k <- min(n, p)

    set.seed(1)
    test_tensor <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    unsupported_mode <- "crap"
    expect_error(
      tpls(test_tensor, test_tensor, mode = unsupported_mode),
      "Please ensure mode is one of: canonical, regression, tsvdm"
    )
  }
)

test_that(
  "tpls propagates names appropriately",
  code = {
    n <- 3
    p <- 8
    q <- 13
    t <- 4
    k <- min(n, p, q)

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- array(rnorm(n * q * t, mean = 0, sd = 3), dim = c(n, q, t))

    dimnames(test_x) <- list(
      paste0("r", seq_len(n)),
      paste0("c_x", seq_len(p)),
      paste0("t", seq_len(t))
    )
    dimnames(test_y) <- list(
      paste0("r", seq_len(n)),
      paste0("c_y", seq_len(q)),
      paste0("t", seq_len(t))
    )

    tpls_mat <- tpls(test_x, test_y, matrix_output = TRUE)
    tpls_tens <- tpls(test_x, test_y, mode = "tsvdm", matrix_output = FALSE)

    # rows of loadings output named after the features of the original data
    expect_equal(
      dimnames(tpls_mat$x_loadings),
      list(paste0("c_x", seq_len(p)), NULL)
    )
    expect_equal(
      dimnames(tpls_tens$x_loadings),
      list(paste0("c_x", seq_len(p)), NULL, paste0("t", seq_len(t)))
    )
    expect_equal(
      dimnames(tpls_mat$y_loadings),
      list(paste0("c_y", seq_len(q)), NULL)
    )
    expect_equal(
      dimnames(tpls_tens$y_loadings),
      list(paste0("c_y", seq_len(q)), NULL, paste0("t", seq_len(t)))
    )

    # rows of the variates (projections) output named after the samples of the
    # original data
    expect_equal(
      dimnames(tpls_mat$x_projected),
      list(paste0("r", seq_len(n)), NULL)
    )
    expect_equal(
      dimnames(tpls_tens$x_projected),
      list(paste0("r", seq_len(n)), NULL, paste0("t", seq_len(t)))
    )
    expect_equal(
      dimnames(tpls_mat$y_projected),
      list(paste0("r", seq_len(n)), NULL)
    )
    expect_equal(
      dimnames(tpls_tens$y_projected),
      list(paste0("r", seq_len(n)), NULL, paste0("t", seq_len(t)))
    )
  }
)

test_that(
  "tsvdm-tpls matrix output keeps the largest singular values",
  code = {
    n <- 6
    p <- 5
    q <- 4
    t <- 3
    ncomp <- 3

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- array(rnorm(n * q * t, mean = 0, sd = 3), dim = c(n, q, t))

    transforms <- dctii_m_transforms(t)
    m <- transforms$m

    # center = FALSE so the singular values can be recomputed directly
    tpls_tsvdm <- tpls(
      test_x, test_y,
      ncomp = ncomp,
      mode = "tsvdm",
      center = FALSE
    )

    svals <- tsvdm(
      ft(m(test_x)) %fp% m(test_y),
      transform = FALSE,
      svals_matrix_form = TRUE
    )$s
    expect_equal(
      svals[cbind(tpls_tsvdm$features, tpls_tsvdm$faces)],
      sort(svals, decreasing = TRUE)[1:ncomp]
    )
  }
)

test_that(
  "tpls accepts a single-column response",
  code = {
    n <- 12
    p <- 5
    t <- 4
    ncomp <- 2

    set.seed(1)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    test_y <- array(rnorm(n * t), dim = c(n, 1, t))

    for (mode in c("regression", "canonical", "tsvdm")) {
      tpls_obj <- tpls(test_x, test_y, ncomp = ncomp, mode = mode)
      expect_equal(dim(tpls_obj$x_projected), c(n, ncomp))
      expect_equal(dim(tpls_obj$y_loadings), c(1, ncomp))
      expect_true(all(is.finite(tpls_obj$y_projected)))
    }
  }
)

test_that(
  "tpls stops at the last component that t(X) Y can support",
  code = {
    n <- 20
    p <- 8
    t <- 4

    set.seed(1)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    # the same two columns at every time point: rank 2, so canonical mode has
    # nothing left to extract after two components
    test_y <- array(rep(rnorm(n * 2), t), dim = c(n, 2, t))

    expect_warning(
      tpls_obj <- tpls(test_x, test_y, ncomp = 4, mode = "canonical"),
      "Only 2 of the 4 requested components"
    )
    expect_equal(tpls_obj$ncomp, 2)
    expect_equal(dim(tpls_obj$x_loadings), c(p, 2))

    # tsvdm mode keeps only the two non-zero singular values; a noise component
    # would sit on a face where y is zero, so its y scores would be zero
    expect_warning(
      tpls_tsvdm <- tpls(test_x, test_y, ncomp = 4, mode = "tsvdm"),
      "Only 2 of the 4 requested components"
    )
    expect_true(all(sqrt(colSums(tpls_tsvdm$y_projected^2)) > 1e-8))

    # a response that varies over time supports all four components
    test_y_varying <- array(rnorm(n * 2 * t), dim = c(n, 2, t))
    expect_warning(
      tpls(test_x, test_y_varying, ncomp = 4, mode = "canonical"),
      NA
    )
  }
)

test_that(
  "tpls leaves a response fixed per subject untransformed",
  code = {
    n <- 10
    p <- 6
    t <- 3

    set.seed(1)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    # a rank 1 response with no time dimension
    z <- rnorm(n)
    test_y <- cbind(z, 2 * z)

    tpls_mat <- tpls(test_x, test_y, ncomp = 1, mode = "canonical")
    expect_false(tpls_mat$transform_y)

    # the same response as a tensor constant in time, not to be transformed
    tpls_tens <- tpls(
      test_x, array(test_y, dim = c(n, 2, t)),
      ncomp = 1, mode = "canonical", transform_y = FALSE
    )
    expect_identical(tpls_mat$x_loadings, tpls_tens$x_loadings)

    # the first component exhausts a rank 1 response in canonical mode; as it is
    # deflated on every face, no face has any of it left for a second
    expect_warning(
      tpls(test_x, test_y, ncomp = 2, mode = "canonical"),
      "Only 1 of the 2 requested components"
    )

    test_y_varying <- array(rnorm(n * 2 * t), dim = c(n, 2, t))
    expect_error(
      tpls(test_x, test_y_varying, transform_y = FALSE),
      "transform_y = FALSE needs a response that is fixed per subject"
    )
  }
)

test_that(
  "tpls with scaled faces chooses them on a common scale",
  code = {
    n <- 20
    p <- 10
    t <- 4

    set.seed(1)
    groups <- rep(c("a", "b"), each = n / 2)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    for (j in 1:5) {
      test_x[, j, ] <- test_x[, j, ] +
        outer(ifelse(groups == "a", 1, -1), c(-1.5, -0.5, 0.5, 1.5))
    }
    test_x <- test_x + array(rnorm(n * p, sd = 3), dim = c(n, p, t))
    test_y <- cbind(groups == "a", groups == "b") * 1

    # the same data with face 1 multiplied by 10, after the transform
    transforms <- dctii_m_transforms(t)
    xhat <- transforms$m(test_x)
    xhat[, , 1] <- 10 * xhat[, , 1]
    test_x_10 <- transforms$minv(xhat)

    # scaled faces choose the same faces; raw scores move to face 1
    expect_equal(
      tpls(test_x, test_y, ncomp = 2)$faces,
      tpls(test_x_10, test_y, ncomp = 2)$faces
    )
    expect_false(identical(
      tpls(test_x, test_y, ncomp = 2, face_scaling = "none")$faces,
      tpls(test_x_10, test_y, ncomp = 2, face_scaling = "none")$faces
    ))

    # a response that changes over time keeps raw scores by default
    test_y_varying <- array(rnorm(n * 2 * t), dim = c(n, 2, t))
    expect_equal(tpls(test_x, test_y_varying, ncomp = 1)$face_scaling, "none")

    expect_error(
      tpls(test_x, test_y, face_scaling = "max"),
      "Please ensure face_scaling is one of: frobenius, none"
    )
  }
)

test_that(
  "tpls returns what is needed to project new samples",
  code = {
    n <- 15
    p <- 20
    # canonical mode deflates a subject-level y by its own variates, so it
    # supports at most q components
    q <- 5
    t <- 4

    set.seed(1)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))
    test_y_tens <- array(rnorm(n * q * t), dim = c(n, q, t))
    test_y_mat <- matrix(rnorm(n * q), n, q)

    for (curr_mode in c("regression", "canonical")) {
      for (test_y in list(test_y_tens, test_y_mat)) {
        tpls_obj <- tpls(test_x, test_y, ncomp = 5, mode = curr_mode)
        # components that share a face test the deflation coefficients, not
        # just the loadings
        expect_true(anyDuplicated(tpls_obj$faces) > 0)

        expect_true(tpls_obj$center)
        expect_equal(tpls_obj$mean_slice_x, apply(test_x, c(2, 3), mean))
        # the transform pair, minv undoing m
        expect_equal(tpls_obj$minv(tpls_obj$m(test_x)), test_x)
        expect_equal(
          tpls_obj$mean_slice_y,
          if (is.matrix(test_y)) matrix(colMeans(test_y), q, t) else
            apply(test_y, c(2, 3), mean)
        )
        expect_equal(dim(tpls_obj$x_reg_coefs), c(p, 5))
        # in canonical mode Y is deflated by its own variates, so there are no
        # coefficients that predict it from X
        if (curr_mode == "regression") {
          expect_equal(dim(tpls_obj$y_reg_coefs), c(q, 5))
        } else {
          expect_null(tpls_obj$y_reg_coefs)
        }

        # from the raw data and the stored quantities alone, the components of
        # each face project with its rotation W (t(P) W)^-1
        xhat <- tpls_obj$m(
          sweep(test_x, c(2, 3), STATS = tpls_obj$mean_slice_x, FUN = "-")
        )
        x_projected <- matrix(0, n, tpls_obj$ncomp)
        for (f in unique(tpls_obj$faces)) {
          idx <- which(tpls_obj$faces == f)
          w <- tpls_obj$x_loadings[, idx, drop = FALSE]
          x_projected[, idx] <- matrix(xhat[, , f], nrow = n) %*%
            w %*% solve(crossprod(tpls_obj$x_reg_coefs[, idx, drop = FALSE], w))
        }
        expect_equal(x_projected, unname(tpls_obj$x_projected))
      }
    }

    tpls_tsvdm <- tpls(test_x, test_y_tens, ncomp = 2, mode = "tsvdm")
    expect_null(tpls_tsvdm$x_reg_coefs)
    expect_null(tpls_tsvdm$y_reg_coefs)
  }
)

test_that(
  "tpls deflation coefficients agree with mixOmics pls at a single time point",
  code = {
    n <- 12
    p <- 8
    q <- 3

    set.seed(1)
    test_x <- matrix(
      rnorm(n * p), n, p,
      dimnames = list(NULL, paste0("f", seq_len(p)))
    )
    test_y <- matrix(rnorm(n * q), n, q)

    mixomics_pls <- mixOmics::pls(test_x, test_y, ncomp = 3, scale = FALSE)
    mixomics_fitted <- predict(mixomics_pls, test_x)$predict[, , 3]

    # the response transformed (as a tensor) and untransformed (as a matrix)
    for (y_input in list(array(test_y, dim = c(n, q, 1)), test_y)) {
      tpls_obj <- tpls(array(test_x, dim = c(n, p, 1)), y_input, ncomp = 3)

      # the X coefficients change sign with the component
      signs <- sign(colSums(tpls_obj$x_loadings * mixomics_pls$loadings$X))
      expect_equal(
        unname(sweep(tpls_obj$x_reg_coefs, 2, signs, "*")),
        unname(mixomics_pls$mat.c)
      )

      # the fitted response, from the variates and the Y coefficients
      fitted <- sweep(
        tcrossprod(tpls_obj$x_projected, tpls_obj$y_reg_coefs),
        2, STATS = tpls_obj$mean_slice_y[, 1], FUN = "+"
      )
      expect_equal(unname(fitted), unname(mixomics_fitted))
    }
  }
)
