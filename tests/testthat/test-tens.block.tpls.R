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
        design = "full",
        pairwise_scoring = pairwise_scoring_test,
        solve_dual = TRUE
      )
      block_tpls_obj_primal_only <- block_tpls(
        test_a,
        test_y_ind,
        ncomp = ncomp_input,
        design = "full",
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

    for (curr_design in c("full", "null")) {
      mixomics_block_pls <- mixOmics::block.pls(
        test_a_mats,
        indY = y_ind,
        ncomp = ncomp_input,
        design = curr_design,
        scale = FALSE
      )
      block_tpls_obj <- block_tpls(
        test_a_tens,
        y_ind,
        ncomp = ncomp_input,
        design = curr_design
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

    for (curr_design in c("full", "null")) {
      mixomics_block_pls <- mixOmics::block.pls(
        test_a_mats,
        indY = y_ind,
        ncomp = ncomp_input,
        mode = "canonical",
        design = curr_design,
        scale = FALSE
      )
      block_tpls_obj <- block_tpls(
        test_a_tens,
        y_ind,
        ncomp = ncomp_input,
        mode = "canonical",
        design = curr_design
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
        design = "full",
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
      ncomp = ncomp_input,
      design = "full"
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

test_that(
  "block_tpls face selection weights each pair of blocks by its design entry",
  code = {
    n <- 20
    t <- 4
    ncomp_input <- 2

    set.seed(1)
    # the two x blocks share a pattern that alternates over time, so their
    # pairwise score is large on a face other than the first; y is the same at
    # every time point, so after the transform it is non-zero on face 1 only
    shared <- outer(matrix(rnorm(n * 6), n, 6), 10 * c(1, -1, 1, -1))
    test_a <- list(
      shared + array(rnorm(n * 6 * t), dim = c(n, 6, t)),
      shared + array(rnorm(n * 6 * t), dim = c(n, 6, t)),
      array(rep(rnorm(n * 2), t), dim = c(n, 2, t))
    )

    # with design = "null" the x-x pair carries no weight, so only face 1,
    # where y is non-zero, can be selected
    block_tpls_obj <- block_tpls(
      test_a, 3,
      ncomp = ncomp_input,
      design = "null"
    )
    expect_equal(as.vector(block_tpls_obj$faces), rep(1, ncomp_input))
  }
)

test_that(
  "block_tpls and block_tplsda require a design",
  code = {
    n <- 7
    ps <- c(10, 5, 2)
    t <- 4

    set.seed(1)
    test_a <- lapply(
      ps,
      function(p) array(rnorm(n * p * t), dim = c(n, p, t))
    )

    expect_error(block_tpls(test_a, 3), "design must be explicitly set")
    expect_error(
      block_tplsda(test_a[1:2], rep(c("a", "b"), length.out = n)),
      "design must be explicitly set"
    )
  }
)

test_that(
  "block_tpls returns what is needed to project new samples",
  code = {
    n <- 15
    q <- 3
    t <- 4

    set.seed(1)
    # the second block has more features than samples, so it is solved in
    # its dual form
    test_x <- list(
      array(rnorm(n * 12 * t), dim = c(n, 12, t)),
      array(rnorm(n * 20 * t), dim = c(n, 20, t))
    )
    test_y_tens <- array(rnorm(n * q * t), dim = c(n, q, t))
    test_y_mat <- matrix(rnorm(n * q), n, q)

    for (curr_mode in c("regression", "canonical")) {
      for (test_y in list(test_y_tens, test_y_mat)) {
        block_obj <- block_tpls(
          c(test_x, list(test_y)), 3,
          ncomp = 4, mode = curr_mode, design = 0.5
        )
        # components that share a face test the deflation coefficients, not
        # just the loadings
        expect_true(anyDuplicated(block_obj$faces) > 0)
        expect_true(block_obj$center)
        # the transform pair, minv undoing m
        expect_equal(block_obj$minv(block_obj$m(test_x[[1]])), test_x[[1]])
        expect_equal(
          block_obj$mean_slices[[3]],
          if (is.matrix(test_y)) matrix(colMeans(test_y), q, t) else
            apply(test_y, c(2, 3), mean)
        )
        expect_null(block_obj$x_reg_coefs[[3]])
        if (curr_mode == "regression") {
          expect_null(block_obj$y_reg_coefs[[3]])
          for (b in 1:2) {
            expect_equal(dim(block_obj$y_reg_coefs[[b]]), c(q, 4))
          }
        } else {
          expect_null(block_obj$y_reg_coefs)
        }

        # from the raw data and the stored quantities alone, the components of
        # each face project with that block's rotation W (t(P) W)^-1
        for (b in 1:2) {
          expect_equal(
            block_obj$mean_slices[[b]],
            apply(test_x[[b]], c(2, 3), mean)
          )
          xhat <- block_obj$m(sweep(
            test_x[[b]], c(2, 3),
            STATS = block_obj$mean_slices[[b]], FUN = "-"
          ))
          projected <- matrix(0, n, block_obj$ncomp)
          for (f in unique(block_obj$faces)) {
            idx <- which(block_obj$faces == f)
            w <- block_obj$loadings[[b]][, idx, drop = FALSE]
            p_coefs <- block_obj$x_reg_coefs[[b]][, idx, drop = FALSE]
            projected[, idx] <- matrix(xhat[, , f], nrow = n) %*%
              w %*% solve(crossprod(p_coefs, w))
          }
          expect_equal(projected, unname(block_obj$projected[[b]]))
        }
      }
    }
  }
)

test_that(
  "block_tpls with one X block has the deflation coefficients of tpls",
  code = {
    n <- 15
    p <- 10
    q <- 2
    t <- 4

    set.seed(1)
    test_x <- array(rnorm(n * p * t), dim = c(n, p, t))

    # with two blocks the RGCCA step is the PLS step; a strict tolerance lets
    # it converge to the same weights
    for (test_y in list(
      array(rnorm(n * q * t), dim = c(n, q, t)),
      matrix(rnorm(n * q), n, q)
    )) {
      tpls_obj <- tpls(test_x, test_y, ncomp = 3)
      block_obj <- block_tpls(
        list(test_x, test_y), 2,
        ncomp = 3, design = "full", tol = 1e-14, max_iter = 1e5
      )
      expect_equal(as.vector(block_obj$faces), as.vector(tpls_obj$faces))

      # the coefficients change sign with the component
      signs <- sign(colSums(block_obj$loadings[[1]] * tpls_obj$x_loadings))
      expect_equal(
        unname(sweep(block_obj$x_reg_coefs[[1]], 2, signs, "*")),
        unname(tpls_obj$x_reg_coefs),
        tolerance = 1e-6
      )
      expect_equal(
        unname(sweep(block_obj$y_reg_coefs[[1]], 2, signs, "*")),
        unname(tpls_obj$y_reg_coefs),
        tolerance = 1e-6
      )
    }
  }
)

test_that(
  "block_tpls names its coefficients after the blocks and their features",
  code = {
    n <- 10
    t <- 3

    set.seed(1)
    test_a <- list(
      A = array(rnorm(n * 6 * t), dim = c(n, 6, t),
                dimnames = list(NULL, paste0("a", 1:6), NULL)),
      B = array(rnorm(n * 4 * t), dim = c(n, 4, t),
                dimnames = list(NULL, paste0("b", 1:4), NULL)),
      Y = array(rnorm(n * 2 * t), dim = c(n, 2, t),
                dimnames = list(NULL, c("y1", "y2"), NULL))
    )
    block_obj <- block_tpls(test_a, 3, ncomp = 2, design = "full")

    expect_equal(names(block_obj$x_reg_coefs), c("A", "B", "Y"))
    expect_equal(names(block_obj$y_reg_coefs), c("A", "B", "Y"))
    expect_equal(rownames(block_obj$x_reg_coefs$A), paste0("a", 1:6))
    expect_equal(rownames(block_obj$y_reg_coefs$B), c("y1", "y2"))
  }
)
