context("tsvdm")
# bltodo: add tests for permutations of different configurations and input sizes
# like the `tred` library

test_that(
  "reconstructing original matrix from tsvdm",
  code = {
    test_tensor <- array(1:24, dim = c(2, 4, 3))
    tsvdm_decomposition <- tsvdm(test_tensor)
    u <- tsvdm_decomposition$u
    s <- tsvdm_decomposition$s
    v <- tsvdm_decomposition$v
    vt <- ft(v)

    expect_equal(
      test_tensor,
      m_product(u, s, vt)
    )
  }
)

test_that(
  "reconstructing original matrix from tsvdm when transform=FALSE",
  code = {
    test_tensor <- array(1:24, dim = c(2, 4, 3))
    tsvdm_decomposition <- tsvdm(test_tensor, transform = FALSE)
    u <- tsvdm_decomposition$u
    s <- tsvdm_decomposition$s
    v <- tsvdm_decomposition$v
    vt <- ft(v)

    expect_equal(
      test_tensor,
      facewise_product(u, s, vt)
    )
  }
)

test_that(
  "`svals_matrix_form = TRUE` gives consistent output",
  code = {
    test_tensor <- array(1:24, dim = c(2, 4, 3))

    tsvdm_decomposition_tensor <- tsvdm(test_tensor)
    s_tensor <- tsvdm_decomposition_tensor$s
    tsvdm_decomposition_matrix <- tsvdm(test_tensor, svals_matrix_form = TRUE)
    s_matrix <- tsvdm_decomposition_matrix$s

    s_tensor_reconstruct <- array(0, dim = dim(s_tensor))
    for (i in seq_len(dim(s_tensor)[3])) {
      s_tensor_reconstruct[, , i] <- diag(s_matrix[, i])
    }
    expect_equal(s_tensor, s_tensor_reconstruct)
  }
)
