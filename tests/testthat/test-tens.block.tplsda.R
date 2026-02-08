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

context("block tplsda")

test_that(
  "block tplsda single time point y imputation is consistent",
  code = {
    n <- 7
    ps <- c(10, 5, 2)
    t <- 4

    y_ind <- length(ps) + 1

    set.seed(1)
    test_x <- lapply(
      ps,
      function(p) array(rnorm(n * p * t, mean = 0, sd = p), dim = c(n, p, t))
    )
    test_y_vec <- letters[1:n]
    test_y_mat <- array(test_y_vec, dim = c(n, 1))
    test_y_tens <- array(test_y_vec, dim = c(n, 1, 1))

    block_tplsda_vec <- block_tplsda(test_x, test_y_vec)
    block_tplsda_mat <- block_tplsda(test_x, test_y_mat)
    block_tplsda_tens <- block_tplsda(test_x, test_y_tens)

    expect_equal(block_tplsda_vec$a[[y_ind]], block_tplsda_mat$a[[y_ind]])
    expect_equal(block_tplsda_mat$a[[y_ind]], block_tplsda_tens$a[[y_ind]])
  }
)

test_that(
  "block tplsda multi_label y transformation is appropriate",
  code = {
    # there's only 26 letters in the alphabet, stored in `letters` variable,
    # so this test will automatically fail if n * t > 26
    n <- 7
    ps <- c(10, 5, 2)
    t <- 3

    y_ind <- length(ps) + 1

    set.seed(1)
    test_x <- lapply(
      ps,
      function(p) array(rnorm(n * p * t, mean = 0, sd = p), dim = c(n, p, t))
    )
    test_y_mat <- array(letters[1:(n * t)], dim = c(n, t))
    test_y_tens <- array(letters[1:(n * t)], dim = c(n, 1, t))

    block_tplsda_mat <- block_tplsda(test_x, test_y_mat, multi_label = TRUE)
    block_tplsda_tens <- block_tplsda(test_x, test_y_tens, multi_label = TRUE)

    # both multi_label y inputs should lead to the same underlying tpls call
    expect_equal(block_tplsda_mat$a[[y_ind]], block_tplsda_tens$a[[y_ind]])

    # the y tensor that goes into the tpls call should have the appropriate
    # number of columns corresponding to all unique categories found across
    # all timepoints specified in y
    expect_equal(dim(block_tplsda_mat$a[[y_ind]]), c(n, n * t, t))
  }
)
