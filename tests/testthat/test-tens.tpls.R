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
