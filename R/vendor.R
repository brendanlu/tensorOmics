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

# ==============================================================================
# Various internal utility files vendored from external packages
# Credit to
# > mixOmics:
#  - tau.estimate
#  - unmap
# > RGCCA:
#  - cov2
# ==============================================================================

#' Vendored from https://github.com/mixOmicsTeam/mixOmics on 16/10/2025
#'
#' Estimation of tau accoring to Strimmer formula
#' @author mixOmics Team 2018
#' @keywords internal
#' @noRd
tau.estimate <- function (x)
{
    if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
        stop("The data matrix must be numeric!")
    
    p = NCOL(x)
    n = NROW(x)
    #covm = cov(x)
    corm = stats::cor(x)
    xs = scale(x, center = TRUE, scale = TRUE)
    xs2 = xs^2
    v = (n/((n - 1)^3)) * (crossprod(xs2) - 1/n * (crossprod(xs))^2)
    diag(v) = 0
    m = matrix(rep(apply(xs2, 2, mean), p), p, p)
    I = diag(NCOL(x))
    d = (corm - I)^2
    tau = (sum(v))/sum(d)
    tau = max(min(tau, 1), 0)
    return(tau)
}

#' Vendored from https://github.com/mixOmicsTeam/mixOmics on 16/10/2025
#'
#' Dummy matrix for an outcome factor; converts a class or group vector or
#' factor into a matrix of indicator variables.
#' @author Ignacio Gonzalez, Kim-Anh Le Cao, Pierre Monget, AL J Abadi
#' @references
#' C. Fraley and A. E. Raftery (2002). Model-based
#' clustering, discriminant analysis, and density estimation. \emph{Journal of
#' the American Statistical Association 97:611-631}.
#'
#' C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012). mclust Version
#' 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification,
#' and Density Estimation. Technical Report No. 597, Department of Statistics,
#' University of Washington.
#' @keywords internal
#' @noRd
unmap <- function(classification, groups = NULL, noise = NULL)
    {
        n = length(classification)
        u = sort(unique(classification))
        levels =  levels(classification)### Add levels
        
        if (is.null(groups))
        {
            groups = u
        } else {
            if (any(match(u, groups, nomatch = 0) == 0))
                stop("groups incompatible with classification")
            miss = match(groups, u, nomatch = 0) == 0
        }
        
        cgroups = as.character(groups)
        if (!is.null(noise))
        {
            noiz = match(noise, groups, nomatch = 0)
            if (any(noiz == 0))
                stop("noise incompatible with classification")
            
            groups = c(groups[groups != noise], groups[groups == noise])
            noise = as.numeric(factor(as.character(noise), levels = unique(groups)))
        }
        
        groups = as.numeric(factor(cgroups, levels = unique(cgroups)))
        classification = as.numeric(factor(as.character(classification), levels = unique(cgroups)))
        k = length(groups) - length(noise)
        nam = levels(groups)
        
        if (!is.null(noise))
        {
            k = k + 1
            nam = nam[1:k]
            nam[k] = "noise"
        }
        
        z = matrix(0, n, k, dimnames = c(names(classification), nam))
        for (j in 1:k) z[classification == groups[j], j] = 1
        attr(z, "levels") = levels
        z
    }

#' Vendored from https://github.com/rgcca-factory/RGCCA on 16/10/2025
#'
#' Covariance; cov2() is similar to stats::cov() but has an additional argument.
#' The denominator \eqn{n} (bias = TRUE) can be used (instead of \eqn{n-1})
#' to give a biased estimator of the (co)variance.
#' Note that if values are missing the "pairwise.complete.obs"
#' option is used resulting in (co)-variance matrix that are not necessarily
#' positive definite.
#' @keywords internal
#' @noRd
cov2 <- function(x, y = NULL, bias = TRUE) {
  n <- NROW(x)

  if (is.null(y)) {
    x <- as.matrix(x)
    if (bias) {
      C <- ((n - 1) / n) * stats::cov(x, use = "pairwise.complete.obs")
    } else {
      C <- stats::cov(x, use = "pairwise.complete.obs")
    }
  } else {
    if (bias) {
      C <- ((n - 1) / n) * stats::cov(x, y, use = "pairwise.complete.obs")
    } else {
      C <- stats::cov(x, y, use = "pairwise.complete.obs")
    }
  }
  return(C)
}
