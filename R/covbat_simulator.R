#' Generate simulations under the CovBat model
#'
#' This function generates simulated data with batch effects assuming the CovBat model.
#'
#' @param n_sims Number of simulated datasets to generate
#' @param covariates A matrix with the covariates to be used to adjust for the batch effects. The first column should be a binary variable indicating the batch (0 or 1).
#' @param effect_sizes A matrix with the effect sizes for each feature for batch and biological covariates.
#' @param pca_effect_sizes A matrix with the effect sizes for each feature for batch and biological covariates in the covariance structure.
#' @param delta A scalar value that determines the strength of the scale batch effect.
#' @param pca_delta A scalar value that determines the strength of the scale batch effect in the covariance.
#' @param cores Number of cores to use. Default is NULL which uses 1 core.
#' @param pca_cutoff A scalar parameter controlling the number of principal components to retain. The default value of 0.95 means that batch effects will only be simulated in the principal components explaining up to 95% of the total variance.
#'
#' @return A list of matrices with the CovBat simulated data.
#' @export
#'
#' @import pracma
#' @import parallel
#'
#' @examples
#' \dontrun{
#' covariates = generate_covariate_matrix(n_subject = 1000, p_batch = 0.6,
#' n_covariates = 3,
#' batchwise_offset = c(0, 2, 3),
#' correlation_strength = "moderate")
#' effect_sizes = generate_effect_sizes(300, 3, c(3, 2, 2, -2), rep(.5, 4))
#' pca_effect_sizes = generate_effect_sizes(300, 3, c(1, -1, 0, 1), rep(.3, 4))
#' generate_covbat_simulations(covariates, effect_sizes, pca_residuals,
#' delta = 1.2)
#' }
generate_covbat_simulations <- function(n_sims, covariates,
                                        effect_sizes, pca_effect_sizes,
                                        delta, pca_delta,
                                        cores = NULL, pca_cutoff = 0.95) {
  if (!is.null(cores)) {
    covbat_sim_list <- mclapply(1:n_sims, function(x) {
      pca_residuals = generate_pca_residuals(n_subjects = nrow(covariates), n_features = ncol(effect_sizes),
                                             covariates, pca_effect_sizes, pca_delta = pca_delta, pca_cutoff = pca_cutoff)
      eigenvector_mat <- randortho(ncol(pca_residuals), type = "orthonormal")
      residuals <- pca_residuals %*% eigenvector_mat
      for (i in 0:1) {
        is_i <- covariates[, 1] == i
        residuals[is_i, ] <- residuals[is_i, ] * (delta^i)
      }

      return(covariates %*% effect_sizes + residuals)
    }, mc.cores = cores)
  } else {
    covbat_sim_list <- lapply(1:n_sims, function(x) {
      pca_residuals = generate_pca_residuals(n_subjects = nrow(covariates), n_features = ncol(effect_sizes),
                                             covariates, pca_effect_sizes, pca_delta = pca_delta, pca_cutoff = pca_cutoff)
      eigenvector_mat <- randortho(ncol(pca_residuals), type = "orthonormal")
      residuals <- pca_residuals %*% eigenvector_mat
      for (i in 0:1) {
        is_i <- covariates[, 1] == i
        residuals[is_i, ] <- residuals[is_i, ] * (delta^i)
      }

      return(covariates %*% effect_sizes + residuals)
    })
  }

  return(covbat_sim_list)
}
