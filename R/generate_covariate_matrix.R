#' Generate simulated batch and covariates-
#'
#' Generates an n x (n_covariates+1) matrix of covariates and batch. Allows for specification of how batch is correlated with covariates through batchwise_offset and how much covariates are correlated with each other.
#'
#' @param n_subject Number of subjects
#' @param p_batch Proportion of subjects in batch 1 ()
#' @param n_covariates Number of biological covariates to preserve (p)
#' @param batchwise_offset Vector describing how much each covariate's mean should be offset in batch 1 compared to batch 0. All covariates are generated under normal distribution with variance of 1.
#' @param correlation_strength "none", "weak", "moderate", or "strong"
#'
#' @return n x (n_covariates+1) matrix of covariates and batch. Batch is the first column.
#' @export
#'
#' @import MASS
#' @import stats
#'
#' @examples
#' \dontrun{
#' generate_covariate_matrix(n_subject = 1000, p_batch = 0.6,
#' n_covariates = 3,
#' batchwise_offset = c(0, 2, 3),
#' correlation_strength = "moderate")
#' }
generate_covariate_matrix <- function(n_subject, p_batch,
                                      n_covariates,
                                      batchwise_offset = NULL,
                                      correlation_strength = "none") {
  match.arg(correlation_strength, c("none", "weak", "moderate", "strong"))
  stopifnot("p_batch must be between 0 and 1"=
              p_batch > 0 & p_batch < 1)
  stopifnot("n_subject and n_covariates must be integers greater than 0"=
              n_subject %% 1 == 0 & n_covariates %% 1 == 0)
  stopifnot("batchwise_offset must be vector of length n_covariates"=
              is.null(batchwise_offset) | length(batchwise_offset) == n_covariates)
  if (is.null(batchwise_offset)) {
    batchwise_offset <- numeric(n_covariates)
  }
  if (correlation_strength == "none") {
    covariance <- diag(n_covariates)
  }
  if (correlation_strength == "weak") {
    covariance <- matrix(c(1, .2, .2, .2, 1, .2, .2, .2, 1), ncol = 3)
  }
  if (correlation_strength == "moderate") {
    covariance <- matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), ncol = 3)
  }
  if (correlation_strength == "strong") {
    covariance <- matrix(c(1, .8, .8, .8, 1, .8, .8, .8, 1), ncol = 3)
  }

  batch <- rbinom(n_subject, 1, p_batch)
  covariates <- matrix(numeric(n_subject * n_covariates), ncol = n_covariates)
  for (i in 0:1) {
    is_i <- batch == i
    covariates[is_i, ] <- mvrnorm(sum(is_i), batchwise_offset * i, covariance)
  }

  return(cbind(batch, covariates))
}
