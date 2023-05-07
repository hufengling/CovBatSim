#' Generate batch effects in the covariance structure
#'
#' This function generates a matrix of PCA residuals with specified batch effects. These batch effects will be rotated to hide them in a covariance matrix.
#'
#' @param n_subjects The number of subjects in the study
#' @param n_features The number of features in the study
#' @param covariates A matrix of size n_subjects x (1 + n_covariates) containing the covariates used to generate the batch effects. The first column should be binary and represent the batch, while the remaining columns should be continuous and represent other covariates that should be controlled for.
#' @param pca_effects A matrix of size (1 + n_covariates) x n_features containing the batch effects and biological effects. Not all effects will be used, since batch effects are only corrected in the first few eigendimensions.
#' @param pca_delta A scalar parameter controlling the magnitude of the batch effects in the PCA residuals. A value of 1 means there is no scale batch effects.
#' @param pca_cutoff A scalar parameter controlling the number of principal components to retain. The default value of 0.95 means that batch effects will only be simulated in the principal components explaining up to 95% of the total variance.
#'
#' @return A matrix of size n_subjects x n_features containing the residuals with the specified batch effects in PC space.
#' @export
#'
#' @importFrom funData eVal
#' @import matrixStats
#' @import MASS
#'
#' @examples
#' \dontrun{
#' pca_effects = generate_effect_sizes(300, 3, c(3, 2, 2, -2), rep(.5, 4))
#' generate_pca_residuals(n_subjects = 1000, n_features = 300, covariates,
#' pca_effects, pca_delta = 1.5, pca_cutoff = 0.95)
#' }
generate_pca_residuals <- function(n_subjects, n_features, covariates,
                                   pca_effects, pca_delta, pca_cutoff = 0.95) {
  eigen_vals <- eVal(n_features, "wiener")
  eigen_vals <- eigen_vals / sum(eigen_vals)
  batch_dim <- cumsum(eigen_vals) < 0.95
  eigen_vals[batch_dim] <- eigen_vals[batch_dim] / 2 # Make sure the final covariance matrix has correct eigenvalues

  pca_means <- covariates %*% pca_effects[, batch_dim]
  pca_means_0 <- (pca_means - matrix(rep(colMeans(pca_means), nrow(pca_means)),
                                     ncol = ncol(pca_means), byrow = TRUE)) / matrix(rep(colSds(pca_means) / (eigen_vals[batch_dim]), nrow(pca_means)),
                                                                                     ncol = ncol(pca_means), byrow = TRUE) # Center and scale the PCs where there are batch effects to mean 0, variance 0.5.
  pca_means_all <- cbind(pca_means_0, matrix(numeric(n_subjects * (n_features - ncol(pca_means))),
                                             nrow = n_subjects))

  tmp_resids <- matrix(numeric(n_subjects * n_features), ncol = n_features)
  for (i in 0:1) {
    is_i <- covariates[, 1] == i
    tmp_resids[is_i, ] <- mvrnorm(sum(is_i), numeric(n_features), diag(eigen_vals))
    tmp_resids[is_i, batch_dim] <- tmp_resids[is_i, batch_dim] * (pca_delta^i)
  }

  return(pca_means_all + tmp_resids)
}
