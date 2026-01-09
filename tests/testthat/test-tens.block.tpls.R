context("block tpls")

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
  "block_tpls obtains the same results using the primal and dual formulations",
  code = {
    n <- 7
    ps <- c(10, 5, 2)
    q <- 12
    t <- 4
    ncomp_input <- 3
    j_total <- length(ps) + 1

    set.seed(1)

    # generate 'blocks' of tensor data and output variable y
    test_a <- lapply(
      c(ps, q),
      function(p) array(rnorm(n * p * t, mean = 0, sd = p), dim = c(n, p, t))
    )
    test_y_ind <- length(ps) + 1

    supported_pairwise_scoring <- c(
      "top_singular",
      "sum_singular",
      "frobenius"
    )

    for (pairwise_scoring_test in supported_pairwise_scoring) {
      block_tpls_obj_dual <- block_tpls(
        test_a,
        test_y_ind,
        ncomp = ncomp_input,
        pairwise_scoring = pairwise_scoring_test,
        solve_dual = TRUE
      )
      block_tpls_obj_primal_only <- block_tpls(
        test_a,
        test_y_ind,
        ncomp = ncomp_input,
        pairwise_scoring = pairwise_scoring_test,
        solve_dual = FALSE
      )

      for (j in seq_len(j_total)) {
        expect_equal(
          block_tpls_obj_dual$loadings[[j]],
          block_tpls_obj_primal_only$loadings[[j]]
        )
        expect_equal(
          block_tpls_obj_dual$projected[[j]],
          block_tpls_obj_primal_only$projected[[j]]
        )
      }
    }
  }
)

test_that(
  "regression: block_tpls agrees with MixOmics block.pls in 2D case",
  code = {
    n <- 7
    ps <- c(10, 5, 4)
    q <- 12
    # BLTODO: need one for regression and separate one for canonical with
    # `ncomp_input` > 1 to test deflation
    ncomp_input <- 3

    j_total <- length(ps) + 1
    ps_merged <- c(ps, q)
    y_ind <- length(ps_merged)
    block_names <- paste0("block", seq_len(length(ps_merged)))

    set.seed(1)
    test_a_mats <- lapply(
      ps_merged,
      function(p) matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_a_mats) <- block_names

    test_a_tens <- lapply(
      test_a_mats,
      function(mat) array(mat, dim = c(nrow(mat), ncol(mat), 1))
    )
    names(test_a_tens) <- block_names

    mixomics_block_pls <- mixOmics::block.pls(
      test_a_mats,
      indY = y_ind,
      ncomp = ncomp_input,
      scale = FALSE
    )
    block_tpls_obj <- block_tpls(
      test_a_tens,
      y_ind,
      ncomp = ncomp_input
    )

    for (j in seq_len(j_total)) {
      expect_equal(
        .make_signs_consistent(mixomics_block_pls$loadings[[j]]),
        .make_signs_consistent(block_tpls_obj$loadings[[j]])
      )
      expect_equal(
        .make_signs_consistent(mixomics_block_pls$variates[[j]]),
        .make_signs_consistent(block_tpls_obj$projected[[j]])
      )
    }
  }
)

test_that(
  "canonical: block_tpls agrees with MixOmics block.pls in 2D case",
  code = {
    n <- 7
    ps <- c(10, 5, 6)
    q <- 12
    # BLTODO: need one for regression and separate one for canonical with
    # `ncomp_input` > 1 to test deflation
    ncomp_input <- 3

    j_total <- length(ps) + 1
    ps_merged <- c(ps, q)
    y_ind <- length(ps_merged)
    block_names <- paste0("block", seq_len(length(ps_merged)))

    set.seed(1)
    test_a_mats <- lapply(
      ps_merged,
      function(p) matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_a_mats) <- block_names

    test_a_tens <- lapply(
      test_a_mats,
      function(mat) array(mat, dim = c(nrow(mat), ncol(mat), 1))
    )
    names(test_a_tens) <- block_names

    mixomics_block_pls <- mixOmics::block.pls(
      test_a_mats,
      indY = y_ind,
      ncomp = ncomp_input,
      mode = "canonical",
      scale = FALSE
    )
    block_tpls_obj <- block_tpls(
      test_a_tens,
      y_ind,
      ncomp = ncomp_input,
      mode = "canonical"
    )

    for (j in seq_len(j_total)) {
      expect_equal(
        .make_signs_consistent(mixomics_block_pls$loadings[[j]]),
        .make_signs_consistent(block_tpls_obj$loadings[[j]])
      )
      expect_equal(
        .make_signs_consistent(mixomics_block_pls$variates[[j]]),
        .make_signs_consistent(block_tpls_obj$projected[[j]])
      )
    }
  }
)

test_that(
  "block_tpls error for non-supported pairwise_scoring input",
  code = {
    n <- 3
    ps <- c(3, 3, 2)
    q <- 1
    t <- 2

    set.seed(1)

    # generate 'blocks' of tensor data and output variable y
    test_x <- lapply(
      c(ps, q),
      function(p) array(rnorm(n * p * t, mean = 0, sd = p), dim = c(n, p, t))
    )
    test_y_ind <- length(ps) + 1

    unsupported_pairwise_scoring <- "crap"
    expect_error(
      block_tpls(
        test_x,
        test_y_ind,
        pairwise_scoring = unsupported_pairwise_scoring
      ),
      # code produces one-liner error message
      paste0(
        "Please ensure facewise_scoring is one of: top_singular, sum_singular,",
        " frobenius"
      )
    )
  }
)

test_that(
  "block_tpls propagates names appropriately",
  code = {
    n <- 7
    ps <- c(10, 5, 2)
    t <- 4
    ncomp_input <- 3

    set.seed(1)
    # generate 'blocks' of tensor data and output variable y
    test_a <- lapply(
      ps,
      function(p) array(
        rnorm(n * p * t, mean = 0, sd = p),
        dim = c(n, p, t),
        dimnames = list(
          paste0("r", seq_len(n)),
          paste0(p, "_c", seq_len(p)),
          paste0("t", seq_len(t))
        )
      )
    )

    block_tpls_obj <- block_tpls(
      test_a,
      1, # doesn't really matter what the "y" block is for this test
      ncomp = ncomp_input
    )

    for (i in seq_along(test_a)) {
      expect_equal(
        dimnames(block_tpls_obj$loadings[[i]]),
        list(paste0(ps[i], "_c", seq_len(ps[i])), NULL)
      )
      expect_equal(
        dimnames(block_tpls_obj$projected[[i]]),
        list(paste0("r", seq_len(n)), NULL)
      )
    }
  }
)
