context("tpca")

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
  "basic tpca sense checks",
  code = {
    n <- 4
    p <- 5
    t <- 3
    ncomp_input <- 2
    test_tensor <- array(1:(n * p * t), dim = c(n, p, t))
    tpca_obj <- tpca(test_tensor, ncomp = ncomp_input)
    # removed explained_variance output pending theory development to ensure
    # it's legit first
    # expect_equal(length(tpca_obj$explained_variance), ncomp_input)
    expect_equal(dim(tpca_obj$variates), c(n, ncomp_input))
    expect_equal(dim(tpca_obj$loadings), c(p, ncomp_input))
  }
)

test_that(
  "tpca agrees with mixOmics pca",
  code = {
    n <- 3
    p <- 10
    t <- 1
    k <- min(n, p)
    ncomp_input <- 2

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    mixomics_pca <- mixOmics::pca(
      test_x[, , 1],
      ncomp = ncomp_input
    )

    tensor_pca <- tpca(
      test_x,
      ncomp = ncomp_input
    )

    expect_equal(
      .make_signs_consistent(mixomics_pca$loadings$X),
      .make_signs_consistent(tensor_pca$loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_pca$variates$X),
      .make_signs_consistent(tensor_pca$variates)
    )
  }
)

test_that(
  "tpca propagates names appropriately",
  code = {
    n <- 3
    p <- 4
    t <- 5

    test_tensor <- array(1:(n * p * t), dim = c(n, p, t))
    dimnames(test_tensor) <- list(
      paste0("r", seq_len(n)),
      paste0("c", seq_len(p)),
      paste0("t", seq_len(t))
    )

    tpca_mat <- tpca(test_tensor, matrix_output = TRUE)
    tpca_tens <- tpca(test_tensor, matrix_output = FALSE)

    # rows of loadings output named after the features of the original data
    expect_equal(
      dimnames(tpca_mat$loadings),
      list(
        paste0("c", seq_len(p)),
        NULL
      )
    )
    expect_equal(
      dimnames(tpca_tens$loadings),
      list(
        paste0("c", seq_len(p)),
        NULL,
        paste0("t", seq_len(t))
      )
    )

    # rows of the variates (projections) output named after the samples of the
    # original data
    expect_equal(
      dimnames(tpca_mat$variates),
      list(
        paste0("r", seq_len(n)),
        NULL
      )
    )
    expect_equal(
      dimnames(tpca_tens$variates),
      list(
        paste0("r", seq_len(n)),
        NULL,
        paste0("t", seq_len(t))
      )
    )
  }
)
