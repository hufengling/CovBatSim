#' Generate effect sizes
#'
#' This function generates a matrix of effect sizes.
#'
#' @param n_features The number of features in the study
#' @param n_covariates The number of covariates in the study (i.e., the number of columns in the design matrix)
#' @param effect_size_means A vector of length n_features + 1 containing the means of the normal distributions from which the effect sizes will be drawn. The first element corresponds to the batch effect, and the remaining elements correspond to the effect sizes for the covariates.
#' @param effect_size_variances A vector of length n_features + 1 containing the variances of the normal distributions from which the effect sizes will be drawn. The first element corresponds to the batch effect, and the remaining elements correspond to the effect sizes for the covariates.
#'
#' @return A matrix of size (n_covariates + 1) x n_features, where the first row contains the intercepts and the remaining rows contain the effect sizes for the predictors.
#' @export
#'
#' @import stats
#'
#' @examples
#' \dontrun{
#' generate_effect_sizes(300, 3, c(3, 2, 2, -2), rep(.5, 4))
#' }
generate_effect_sizes <- function(n_features, n_covariates,
                                  effect_size_means,
                                  effect_size_variances) {
  effect_sizes <- matrix(numeric(n_features * (1 + n_covariates)),
                         ncol = n_features)
  for (i in 1:nrow(effect_sizes)) {
    effect_sizes[i, ] <- rnorm(n_features, effect_size_means[i], effect_size_variances[i])
  }

  return(effect_sizes)
}
