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
