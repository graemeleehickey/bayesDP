#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @title Bayesian Discount Prior: Comparison Between Current and Historical Data
#' @description An educational calculation of the current-historical comparison
#'   for normal summary data. This helper is not used internally. The fixed
#'   method returns twice the smaller directional posterior probability; the
#'   Monte Carlo method returns a normal-tail comparison for each draw.
#' @param mu scalar. Mean of the current data.
#' @param sigma scalar. Standard deviation of the current data.
#' @param N scalar. Number of observations of the current data.
#' @param mu0 scalar. Mean of the historical data.
#' @param sigma0 scalar. Standard deviation of the historical data.
#' @param N0 scalar. Number of observations of the historical data.
#' @param number_mcmc scalar. Number of Monte Carlo simulations. Default is
#'   10000.
#' @param method character. Default \code{"fixed"} estimates a directional
#'   probability q and returns \code{2 * min(q, 1 - q)}. With \code{"mc"},
#'   each pair of mean draws is standardized by its conditional standard error
#'   and converted to a two-sided normal-tail comparison. See the normal
#'   vignette for details.
#' @details This is a borrowing comparison, not a posterior probability of
#'   equality or equivalence. Weight variation under \code{"mc"} is part of
#'   that procedure; more draws improve the precision of its distributional
#'   summaries. The baseline prior for each source is proportional to
#'   \eqn{1 / \sigma^2}. The fixed method uses direct Monte Carlo simulation.
#' @return A numeric comparison in [0, 1]: one value under \code{"fixed"}, or
#'   a vector of length \code{number_mcmc} under \code{"mc"}.
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' probability_discount(
#'   mu = 0, sigma = 1, N = 100,
#'   mu0 = 0.1, sigma0 = 1, N0 = 100
#' )
#' @rdname probability_discount
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset
#'   model.response pweibull quantile rbeta rgamma rnorm var vcov pnorm
#' @aliases probability_discount,ANY-method
#' @export probability_discount
probability_discount <- setClass("probability_discount")

setGeneric(
  "probability_discount",
  function(mu = NULL,
           sigma = NULL,
           N = NULL,
           mu0 = NULL,
           sigma0 = NULL,
           N0 = NULL,
           number_mcmc = 10000,
           method = "fixed") {
    standardGeneric("probability_discount")
  }
)

setMethod(
  "probability_discount",
  signature(),
  function(mu = NULL,
           sigma = NULL,
           N = NULL,
           mu0 = NULL,
           sigma0 = NULL,
           N0 = NULL,
           number_mcmc = 10000,
           method = "fixed") {
    if (!(method %in% c("fixed", "mc"))) {
      stop("method entered incorrectly. Must be one of 'fixed' or 'mc'.")
    }
    
    ### Preposterior of current mu using flat prior
    posterior_flat_sigma2 <- 1 / rgamma(number_mcmc, (N - 1) / 2, ((N - 1) * sigma^2) / 2)
    s <- (posterior_flat_sigma2 / ((N - 1) + 1))^0.5
    posterior_flat_mu <- rnorm(number_mcmc, mu, s)
    
    ### Posterior of historical data parameters using flat prior
    prior_sigma2 <- 1 / rgamma(number_mcmc, (N0 - 1) / 2, ((N0 - 1) * sigma0^2) / 2)
    s0 <- (prior_sigma2 / ((N0 - 1) + 1))^0.5
    prior_mu <- rnorm(number_mcmc, mu0, s0)
    
    ### Test of model vs real
    p_hat <- mean(posterior_flat_mu < prior_mu)
    
    ### Transform probability to an intuitive value
    p_hat <- 2 * ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
    
    if (method == "mc") {
      Z <- abs(posterior_flat_mu - prior_mu) / sqrt(s^2 + s0^2)
      p_hat <- 2 * (1 - pnorm(Z))
    }
    
    return(p_hat)
  }
)
