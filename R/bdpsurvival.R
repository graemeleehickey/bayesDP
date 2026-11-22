#' @title Bayesian Discount Prior: Survival Analysis
#'
#' @description \code{bdpsurvival} is used to estimate the survival probability
#'   (single arm trial; OPC) or hazard ratio (two-arm trial; RCT) for
#'   right-censored data using the survival analysis implementation of the
#'   Bayesian discount prior. Current treatment data are required; control data
#'   in a two-arm analysis can be current, historical, or both. This code is modeled after the methodologies
#'   developed in Haddad et al. (2017).
#'
#' @param formula an object of class "formula." Must have a survival object on
#'   the left side and at most one input on the right side, treatment. See
#'   "Details" for more information. Response expressions such as
#'   \code{Surv(followup / 12, status == 1)} and \code{survival::Surv()}
#'   are supported. Only right-censored responses are supported.
#' @param data a data frame containing the current data variables in the model.
#'   Supply the variables used in the formula. Evaluated times must be finite
#'   and non-negative, with no missing outcomes. If treatment is specified,
#'   use numeric 0/1 or logical values, with no missing values; factors and
#'   character codes are not converted automatically.
#' @param data0 optional. A data frame containing the historical data variables
#'   in the model, subject to the same response and treatment checks as current
#'   data. It must supply the variables used in the same formula.
#' @param breaks numeric vector. Positive, finite internal cutpoints of the
#'   piecewise exponential model, without zero or infinity. Values are sorted
#'   and duplicates removed. Use \code{numeric(0)} for one constant hazard.
#'   Default cutpoints are pooled observed-time quantiles.
#' @param a0 positive, finite scalar. Prior value for the gamma shape of the piecewise
#'   exponential hazards. Default is 0.1.
#' @param b0 positive, finite scalar. Prior value for the gamma rate of the piecewise exponential
#'   hazards. Default is 0.1.
#' @param surv_time positive, finite scalar. Time of interest for computing the
#'   probability of survival for a single arm (OPC) trial. Default is overall,
#'   i.e., current+historical, median observed event-or-censoring time.
#' @param discount_function character. Specify the discount function to use.
#'   Currently supports \code{weibull}, \code{scaledweibull}, and
#'   \code{identity}. The discount function \code{scaledweibull} scales the
#'   output of the Weibull CDF to have a max value of 1. The \code{identity}
#'   discount function uses the transformed comparison directly as the discount
#'   weight. Default value is "\code{identity}".
#' @param alpha_max scalar in [0, 1]. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param number_mcmc integer of at least 2. Number of Monte Carlo simulations.
#'   Default is 10000.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. For a two-arm trial, users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group.
#' @param method character. Analysis method with respect to estimation of the
#'   weight parameter alpha. Default method "\code{mc}" estimates alpha for each
#'   Monte Carlo iteration, giving a distribution of generated weights.
#'   Increasing the number of draws improves numerical precision but does not
#'   remove the spread of that distribution. These weights are not draws from
#'   a separately specified prior/posterior model for alpha.
#'   Alternate value "\code{fixed}" estimates alpha
#'   once and holds it fixed throughout the analysis. See the
#'   \code{bdpsurvival} vignette \cr
#'   \code{vignette("bdpsurvival-vignette", package = "bayesDP")} for more
#'   details.
#' @param compare logical. Should a comparison object be included in the fit?
#'   For a one-arm analysis, the comparison object is simply the posterior draws
#'   of the treatment group parameter. For a two-arm analysis, the comparison
#'   object is the posterior draws of the treatment effect that compares
#'   treatment and control. If \code{compare = TRUE}, the comparison object is
#'   accessible in the \code{final} slot, else the \code{final} slot is
#'   \code{NULL}. Default is \code{TRUE}.
#'
#' @details Historical borrowing is assessed separately within each arm by
#'   comparing survival outcomes. With \code{method = "fixed"}, estimate
#'   a directional probability q comparing current and historical outcomes within an arm.
#'   The discount input is \code{p = 2 * min(q, 1 - q)} and the weight is
#'   \code{alpha_max * W(p)} for the selected discount function W. This
#'   comparison is not a posterior probability of equality or equivalence.
#'   Choosing \code{fix_alpha = TRUE} prespecifies the weight instead.
#'
#'   For one arm, q compares survival at \code{surv_time}; for two arms it
#'   compares current and historical hazards separately within treatment and
#'   control using analytic inverse-variance weights for the interval log
#'   hazard differences. With \code{method = "mc"}, each draw's difference
#'   is divided by its posterior standard deviation and mapped to
#'   \code{2 * pnorm(abs(Z), lower.tail = FALSE)}. The variance is calculated
#'   from exact Gamma moments: hazard variances for log survival differences,
#'   and trigamma log-hazard variances for two-arm comparisons. The normal-tail
#'   mapping remains an approximation, especially with sparse events.
#'   The returned \code{p_hat} uses the selected method's transformation even
#'   with \code{fix_alpha = TRUE}; the prescribed weight is unaffected.
#'
#'   In the second stage, the historical likelihood is raised to the weight
#'   \code{alpha}, and the baseline Gamma prior is included once. For an
#'   interval with current event count D and exposure T, and historical
#'   event count D0 and exposure T0, the augmented Gamma shape and rate are
#'   \code{a0 + D + alpha * D0} and \code{b0 + T + alpha * T0}.
#'   The prior shape \code{a0} is distinct from the historical weight alpha.
#'   At alpha = 0 this gives the current-only posterior; at alpha = 1 it
#'   gives the pooled posterior using the same interval boundaries and baseline
#'   prior. With \code{method = "mc"}, each draw uses its generated weight;
#'   the resulting marginal distribution is a mixture of Gamma distributions.
#'
#'   To carry out a single arm (OPC) analysis, data for the current and
#'   historical treatments are specified in separate data frames, data and
#'   data0, respectively. The formula is evaluated separately in each data
#'   frame to obtain the survival (censor) time and event indicator. Missing
#'   outcomes or treatment values cause an error rather than dropping rows.
#'   Counting-process, interval-censored, and multistate responses are not
#'   supported. The results are then
#'   based on the posterior probability of survival at \code{surv_time} for the
#'   current data augmented by the historical data.
#'
#'   Two-arm (RCT) analyses use the same response specification, with an
#'   additional column named 'treatment' to denote treatment and control data.
#'   Use 1 (or TRUE) for treatment and 0 (or FALSE) for control in each data
#'   frame. The current data are augmented by historical data (if present) and
#'   the results are then based on the posterior distribution of the hazard
#'   ratio between the treatment and control groups.
#'
#'   For more details, see the \code{bdpsurvival} vignette: \cr
#'   \code{vignette("bdpsurvival-vignette", package = "bayesDP")}
#'
#' @return \code{bdpsurvival} returns an object of class "bdpsurvival". The
#'   functions \code{\link[=summary,bdpsurvival-method]{summary}} and
#'   \code{\link[=print,bdpsurvival-method]{print}} are used to obtain and print
#'   a summary of the results, including user inputs. The
#'   \code{\link[=plot,bdpsurvival-method]{plot}} function displays visual
#'   outputs as well.
#'
#' An object of class "\code{bdpsurvival}" is a list containing at least
#' the following components:
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the treatment group:
#'    \itemize{
#'      \item \code{alpha_discount}
#'        numeric. Alpha value, the weighting parameter of the historical data.
#'      \item \code{p_hat}
#'        numeric. The two-sided current-historical comparison in [0, 1]:
#'        one value for \code{method = "fixed"}, or \code{number_mcmc} values
#'        for \code{method = "mc"}. See Details for the transformation.
#'      \item \code{posterior_survival}
#'        matrix. If one-arm trial, one column with \code{number_mcmc} rows
#'        containing the posterior probability draws of survival at
#'        \code{surv_time}.
#'      \item \code{posterior_flat_survival}
#'        matrix. If one-arm trial, one column with \code{number_mcmc} rows
#'        containing the probability draws of survival at \code{surv_time}
#'        for the current treatment not augmented by historical treatment.
#'      \item \code{prior_survival}
#'        matrix. If one-arm trial, one column with \code{number_mcmc} rows
#'        containing the probability draws of survival at \code{surv_time}
#'        for the historical treatment.
#'      \item \code{posterior_hazard}
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks) + 1}
#'        columns containing the posterior draws of the piecewise hazards
#'        for each interval break point.
#'      \item \code{posterior_flat_hazard}
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks) + 1}
#'        columns containing the draws of piecewise hazards for each interval
#'        break point for current treatment not augmented by historical treatment.
#'      \item \code{prior_hazard}
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks) + 1}
#'        columns containing the draws of piecewise hazards for each interval break point
#'        for historical treatment.
#'      \item \code{posterior_log_hazard}, \code{posterior_flat_log_hazard},
#'        \code{prior_log_hazard}
#'        matrices. Log hazards corresponding to the three hazard matrices.
#'        These retain finite log draws when very small hazards round to zero
#'        and are used for borrowing comparisons and the two-arm effect.
#'    }
#'  }
#'  \item{\code{posterior_control}}{
#'    list. If two-arm trial, contains values related to the control group
#'    analogous to the \code{posterior_treatment} output.}
#'  \item{\code{final}}{
#'    list. Contains the final comparison object, dependent on the analysis type:
#'    \itemize{
#'      \item One-arm analysis:
#'        vector. Posterior chain of survival probability at requested time.
#'      \item Two-arm analysis:
#'        vector. Posterior chain of log-hazard rate comparing treatment and control groups.
#'    }
#'  }
#'  \item{\code{args1}}{
#'    list. Entries contain user inputs. In addition, the following elements
#'    are output:
#'    \itemize{
#'      \item \code{S_t}, \code{S_c}, \code{S0_t}, \code{S0_c}
#'        survival objects. Used internally to pass survival data between
#'        functions.
#'      \item \code{arm2}
#'        logical. Used internally to indicate one-arm or two-arm analysis.
#'      \item \code{participants}
#'        data frame. One record per participant, with observed \code{time},
#'        event \code{status}, \code{treatment}, and \code{historical} source.
#'        Used for reported sample sizes and risk sets; \code{data} contains
#'        the split interval records used for hazard estimation.
#'    }
#'  }
#' }
#'
#' @seealso \code{\link[=summary,bdpsurvival-method]{summary}},
#'   \code{\link[=print,bdpsurvival-method]{print}},
#'   and \code{\link[=plot,bdpsurvival-method]{plot}} for details of each of the
#'   supported methods.
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' # One-arm trial (OPC) example - survival probability at 5 years
#'
#' # Collect data into data frames
#' df_ <- data.frame(
#'   status = rexp(50, rate = 1 / 30),
#'   time = rexp(50, rate = 1 / 20)
#' )
#' df_$status <- ifelse(df_$time < df_$status, 1, 0)
#'
#' df0 <- data.frame(
#'   status = rexp(50, rate = 1 / 30),
#'   time = rexp(50, rate = 1 / 10)
#' )
#' df0$status <- ifelse(df0$time < df0$status, 1, 0)
#'
#'
#' fit1 <- bdpsurvival(Surv(time, status) ~ 1,
#'   data = df_,
#'   data0 = df0,
#'   surv_time = 5,
#'   method = "fixed"
#' )
#'
#' print(fit1)
#' \dontrun{
#' plot(fit1)
#' }
#'
#' # Two-arm trial example
#' # Collect data into data frames
#' df_ <- data.frame(
#'   time = c(
#'     rexp(50, rate = 1 / 20), # Current treatment
#'     rexp(50, rate = 1 / 10)
#'   ), # Current control
#'   status = rexp(100, rate = 1 / 40),
#'   treatment = c(rep(1, 50), rep(0, 50))
#' )
#' df_$status <- ifelse(df_$time < df_$status, 1, 0)
#'
#' df0 <- data.frame(
#'   time = c(
#'     rexp(50, rate = 1 / 30), # Historical treatment
#'     rexp(50, rate = 1 / 5)
#'   ), # Historical control
#'   status = rexp(100, rate = 1 / 40),
#'   treatment = c(rep(1, 50), rep(0, 50))
#' )
#' df0$status <- ifelse(df0$time < df0$status, 1, 0)
#'
#' fit2 <- bdpsurvival(Surv(time, status) ~ treatment,
#'   data = df_,
#'   data0 = df0,
#'   method = "fixed"
#' )
#' summary(fit2)
#'
#' ### Fix alpha at 1
#' fit2_1 <- bdpsurvival(Surv(time, status) ~ treatment,
#'   data = df_,
#'   data0 = df0,
#'   fix_alpha = TRUE,
#'   method = "fixed"
#' )
#' summary(fit2_1)
#' @rdname bdpsurvival
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset
#'   model.response pweibull quantile rbeta rgamma rnorm var vcov update
#' @importFrom survival Surv survSplit is.Surv
#' @aliases bdpsurvival-method
#' @aliases bdpsurvival,ANY-method
#' @export bdpsurvival
bdpsurvival <- setClass("bdpsurvival", slots = c(
  posterior_treatment = "list",
  posterior_control = "list",
  final = "list",
  args1 = "list"
))
setGeneric(
  "bdpsurvival",
  function(formula = formula,
           data = data,
           data0 = NULL,
           breaks = NULL,
           a0 = 0.1,
           b0 = 0.1,
           surv_time = NULL,
           discount_function = "identity",
           alpha_max = 1,
           fix_alpha = FALSE,
           number_mcmc = 10000,
           weibull_scale = 0.135,
           weibull_shape = 3,
           method = "mc",
           compare = TRUE) {
    standardGeneric("bdpsurvival")
  }
)

setMethod(
  "bdpsurvival",
  signature(),
  function(formula = formula,
           data = data,
           data0 = NULL,
           breaks = NULL,
           a0 = 0.1,
           b0 = 0.1,
           surv_time = NULL,
           discount_function = "identity",
           alpha_max = 1,
           fix_alpha = FALSE,
           number_mcmc = 10000,
           weibull_scale = 0.135,
           weibull_shape = 3,
           method = "mc",
           compare = TRUE) {


    ### Check validity of data input
    call <- match.call()
    if (missing(data)) {
      stop("Current data not input correctly.")
    }
    for (parameter in c("a0", "b0")) {
      value <- get(parameter)
      if (!is.numeric(value) || length(value) != 1L ||
          !is.finite(value) || value <= 0) {
        stop(parameter, " must be a positive, finite scalar.", call. = FALSE)
      }
    }
    if (length(method) != 1L || is.na(method) || !method %in% c("fixed", "mc")) {
      stop("Unrecognized method. Use one of 'fixed' or 'mc'.", call. = FALSE)
    }
    if (!is.numeric(number_mcmc) || length(number_mcmc) != 1L ||
        !is.finite(number_mcmc) || number_mcmc < 2 || number_mcmc != floor(number_mcmc)) {
      stop("number_mcmc must be an integer of at least 2.", call. = FALSE)
    }
    if (!is.numeric(alpha_max) || !length(alpha_max) %in% c(1L, 2L) ||
        any(!is.finite(alpha_max)) || any(alpha_max < 0 | alpha_max > 1)) {
      stop("alpha_max must contain one or two finite values in [0, 1].", call. = FALSE)
    }
    if (!is.logical(fix_alpha) || length(fix_alpha) != 1L || is.na(fix_alpha)) {
      stop("fix_alpha must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(surv_time) && (!is.numeric(surv_time) || length(surv_time) != 1L ||
        !is.finite(surv_time) || surv_time <= 0)) {
      stop("surv_time must be a positive, finite scalar.", call. = FALSE)
    }

    ##############################################################################
    ### Parse current data
    ##############################################################################

    dataCurrent <- bdp_survival_data(formula, data, "Current", 0)
    dataHistorical <- if (!is.null(data0)) {
      bdp_survival_data(formula, data0, "Historical", 1)
    } else {
      NULL
    }

    # Check that discount_function is input correctly
    all_functions <- c("weibull", "scaledweibull", "identity")
    if (length(discount_function) != 1L || is.na(discount_function) ||
        !discount_function %in% all_functions) {
      stop("discount_function input incorrectly.")
    }
    if (discount_function != "identity" && !fix_alpha) {
      for (parameter in c("weibull_shape", "weibull_scale")) {
        value <- get(parameter)
        if (!is.numeric(value) || !length(value) %in% c(1L, 2L) ||
            any(!is.finite(value)) || any(value <= 0)) {
          stop(parameter, " must contain one or two positive, finite values.", call. = FALSE)
        }
      }
    }

    historical <- NULL
    treatment <- NULL

    ##############################################################################
    # Quick check, if alpha_max, weibull_scale, or weibull_shape have length 1,
    # repeat input twice
    ##############################################################################

    if (length(alpha_max) == 1) {
      alpha_max <- rep(alpha_max, 2)
    }

    if (length(weibull_scale) == 1) {
      weibull_scale <- rep(weibull_scale, 2)
    }

    if (length(weibull_shape) == 1) {
      weibull_shape <- rep(weibull_shape, 2)
    }

    ##############################################################################
    # Format input data
    ##############################################################################
    ### If no breaks input, create intervals along quantiles
    if (is.null(breaks)) {
      breaks <- quantile(c(dataCurrent$time, dataHistorical$time), probs = c(0.2, 0.4, 0.6, 0.8))
    }
    if (!is.numeric(breaks) || any(!is.finite(breaks)) || any(breaks < 0)) {
      stop("breaks must contain finite, non-negative cutpoints.", call. = FALSE)
    }
    breaks <- sort(unique(breaks))

    ### If zero is present in breaks, remove and give warning
    if (any(breaks == 0)) {
      breaks <- breaks[!(breaks == 0)]
      warning("Breaks vector included 0. The zero value was removed.")
    }

    ### Combine evaluated responses and retain the original participant records.
    dataALL <- rbind(dataCurrent, dataHistorical)
    participants <- dataALL
    rownames(participants) <- NULL

    # A cut at the last observation also handles a constant-hazard model.
    split_cuts <- if (length(breaks)) breaks else max(dataALL$time)
    dataSplit <- survSplit(Surv(time, status) ~ treatment + historical,
      cut = split_cuts, start = "start", end = "time", event = "status",
      episode = "interval", data = dataALL, na.action = stats::na.fail)

    ### Compute exposure time within each interval
    dataSplit$exposure <- dataSplit$time - dataSplit$start

    ### Retain all intervals, including those with no observations in one source.
    dataSplit$interval <- factor(dataSplit$interval,
                                 levels = seq_len(length(breaks) + 1L),
                                 labels = c(breaks, Inf)
    )

    ### Parse out the historical and current data
    S_t <- subset(dataSplit, historical == 0 & treatment == 1)
    S_c <- subset(dataSplit, historical == 0 & treatment == 0)
    S0_t <- subset(dataSplit, historical == 1 & treatment == 1)
    S0_c <- subset(dataSplit, historical == 1 & treatment == 0)

    if (nrow(S0_t) == 0) S0_t <- NULL
    if (nrow(S_c) == 0) S_c <- NULL
    if (nrow(S0_c) == 0) S0_c <- NULL

    ### Compute arm2, internal indicator of a two-arm trial
    if (is.null(S_c) & is.null(S0_c)) {
      arm2 <- FALSE
    } else {
      arm2 <- TRUE
    }

    # if(arm2) stop("Two arm trials are not currently supported.")

    ### If surv_time is null, replace with median of combined
    ### current + historical times
    if (is.null(surv_time) & !arm2) {
      all_times <- dataALL$time
      surv_time <- median(all_times)
    }
    if (!arm2 && (!is.finite(surv_time) || surv_time <= 0)) {
      stop("surv_time must be a positive, finite scalar.", call. = FALSE)
    }

    ### Check inputs
    if (!arm2) {
      if (nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
      if (is.null(S0_t)) warning("Historical treatment data missing or input incorrectly.")
    } else if (arm2) {
      if (nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
      if (is.null(S_c)) warning("Current control data missing or input incorrectly.")
      if (is.null(S0_t) & is.null(S0_c)) warning("Historical data input incorrectly.")
    }

    posterior_treatment <- posterior_survival(
      S = S_t,
      S0 = S0_t,
      surv_time = surv_time,
      discount_function = discount_function,
      alpha_max = alpha_max[1],
      fix_alpha = fix_alpha,
      a0 = a0,
      b0 = b0,
      number_mcmc = number_mcmc,
      weibull_shape = weibull_shape[1],
      weibull_scale = weibull_scale[1],
      breaks = breaks,
      arm2 = arm2,
      method = method
    )

    if (arm2) {
      posterior_control <- posterior_survival(
        S = S_c,
        S0 = S0_c,
        discount_function = discount_function,
        alpha_max = alpha_max[2],
        fix_alpha = fix_alpha,
        a0 = a0,
        b0 = b0,
        surv_time = surv_time,
        number_mcmc = number_mcmc,
        weibull_shape = weibull_shape[2],
        weibull_scale = weibull_scale[2],
        breaks = breaks,
        arm2 = arm2,
        method = method
      )
    } else {
      posterior_control <- NULL
    }

    args1 <- list(
      S_t = S_t,
      S_c = S_c,
      S0_t = S0_t,
      S0_c = S0_c,
      discount_function = discount_function,
      alpha_max = alpha_max,
      fix_alpha = fix_alpha,
      a0 = a0,
      b0 = b0,
      surv_time = surv_time,
      number_mcmc = number_mcmc,
      weibull_scale = weibull_scale,
      weibull_shape = weibull_shape,
      method = method,
      arm2 = arm2,
      breaks = breaks,
      data = dataSplit,
      participants = participants,
      data_current = data
    )


    ##############################################################################
    ### Create final (comparison) object
    ##############################################################################

    if (!compare) {
      final <- NULL
    } else {
      if (arm2) {
        final <- list()
        final$posterior_loghazard <- bdp_survival_effect(posterior_treatment, posterior_control)
      } else {
        final <- list()
        final$posterior_survival <- posterior_treatment$posterior_survival
      }
    }


    me <- list(
      posterior_treatment = posterior_treatment,
      posterior_control = posterior_control,
      final = final,
      args1 = args1
    )

    class(me) <- "bdpsurvival"

    return(me)
  }
)


################################################################################
# Survival posterior estimation
# 1) Estimate the discount function (if current+historical data both present)
# 2) Estimate the posterior of the augmented data
################################################################################
### Combine  loss function and posterior estimation into one function
posterior_survival <- function(S, S0, surv_time, discount_function,
                               alpha_max, fix_alpha, a0, b0,
                               number_mcmc, weibull_shape, weibull_scale,
                               breaks, arm2, method) {
  nInt <- length(breaks) + 1L
  current <- bdp_survival_source(S, nInt, a0, b0, number_mcmc)
  historical <- bdp_survival_source(S0, nInt, a0, b0, number_mcmc)
  p_hat <- weight <- NULL

  if (!is.null(current) && !is.null(historical)) {
    p_hat <- bdp_survival_comparison(current, historical, surv_time, breaks, arm2, method)
    weight <- if (fix_alpha) alpha_max else alpha_discount(p_hat,
      discount_function = discount_function, alpha_max = alpha_max,
      weibull_shape = weibull_shape, weibull_scale = weibull_scale)

    posterior_log_hazard <- matrix(NA_real_, number_mcmc, nInt)
    for (j in seq_len(nInt)) {
      # Discount historical events and exposure; include the baseline prior once.
      posterior_log_hazard[, j] <- bdp_log_gamma(number_mcmc,
        a0 + current$events[j] + weight * historical$events[j],
        b0 + current$exposure[j] + weight * historical$exposure[j])
    }
  } else {
    source <- if (is.null(current)) historical else current
    posterior_log_hazard <- source$log_hazard
  }

  survival_draws <- function(log_hazard) {
    if (arm2 || is.null(log_hazard)) return(NULL)
    matrix(exp(-exp(bdp_log_cumulative_hazard(log_hazard, surv_time, breaks))), ncol = 1)
  }
  list(
    alpha_discount = weight,
    p_hat = p_hat,
    posterior_survival = survival_draws(posterior_log_hazard),
    posterior_flat_survival = survival_draws(current$log_hazard),
    prior_survival = survival_draws(historical$log_hazard),
    posterior_hazard = exp(posterior_log_hazard),
    posterior_flat_hazard = if (is.null(current)) NULL else exp(current$log_hazard),
    prior_hazard = if (is.null(historical)) NULL else exp(historical$log_hazard),
    posterior_log_hazard = posterior_log_hazard,
    posterior_flat_log_hazard = current$log_hazard,
    prior_log_hazard = historical$log_hazard
  )
}

# Separate Gamma posteriors for each source and interval.
bdp_survival_source <- function(data, nInt, a0, b0, n) {
  if (is.null(data)) return(NULL)
  events <- exposure <- shape <- rate <- numeric(nInt)
  log_hazard <- matrix(NA_real_, n, nInt)
  for (j in seq_len(nInt)) {
    rows <- as.integer(data$interval) == j
    events[j] <- sum(data$status[rows])
    exposure[j] <- sum(data$exposure[rows])
    shape[j] <- a0 + events[j]
    rate[j] <- b0 + exposure[j]
    log_hazard[, j] <- bdp_log_gamma(n, shape[j], rate[j])
  }
  list(events = events, exposure = exposure, shape = shape, rate = rate,
    log_hazard = log_hazard)
}

# For a < 1, Gamma(a, 1) = Gamma(a + 1, 1) * U^(1/a), independently.
# Work on the log scale so small hazards need no artificial positive offset.
bdp_log_gamma <- function(n, shape, rate) {
  if (any(!is.finite(shape) | shape <= 0) || any(!is.finite(rate) | rate <= 0)) {
    stop("Gamma shapes and rates must be positive and finite.", call. = FALSE)
  }
  shape <- rep_len(shape, n)
  rate <- rep_len(rate, n)
  small <- shape < 1
  draws <- log(rgamma(n, shape + small, rate = 1)) - log(rate)
  draws[small] <- draws[small] + log(stats::runif(sum(small))) / shape[small]
  if (any(!is.finite(draws))) {
    stop("Log-Gamma draws exceed numerical precision; use less extreme prior parameters.", call. = FALSE)
  }
  draws
}

# Sum positive quantities while retaining their logarithms.
bdp_log_sum_exp <- function(x) {
  if (!is.matrix(x)) x <- matrix(x, nrow = 1)
  largest <- apply(x, 1, max)
  largest + log(rowSums(exp(sweep(x, 1, largest, "-"))))
}

bdp_interval_lengths <- function(time, breaks) {
  pmax(0, pmin(time, c(breaks, Inf)) - c(0, breaks))
}

bdp_log_cumulative_hazard <- function(log_hazard, time, breaks) {
  lengths <- bdp_interval_lengths(time, breaks)
  active <- lengths > 0
  bdp_log_sum_exp(sweep(log_hazard[, active, drop = FALSE],
    2, log(lengths[active]), "+"))
}

# Comparisons use separate current and historical posteriors, within one arm.
bdp_survival_comparison <- function(current, historical, time, breaks, arm2, method) {
  if (arm2) {
    variance <- trigamma(current$shape) + trigamma(historical$shape)
    if (any(!is.finite(variance) | variance <= 0)) {
      stop("Log-hazard variances exceed numerical precision; use less extreme prior shapes.", call. = FALSE)
    }
    # Rescaling inverse variances avoids overflow without changing the weights.
    weights <- min(variance) / variance
    weights <- weights / sum(weights)
    contrast <- drop((historical$log_hazard - current$log_hazard) %*% weights)
    standard_error <- sqrt(min(variance) / sum(min(variance) / variance))
    Z <- abs(contrast) / standard_error
  } else {
    logH <- bdp_log_cumulative_hazard(current$log_hazard, time, breaks)
    logH0 <- bdp_log_cumulative_hazard(historical$log_hazard, time, breaks)
    # Survival is higher when cumulative hazard is smaller. Comparing log H
    # also distinguishes draws whose survival probabilities round to 0 or 1.
    contrast <- logH0 - logH
    lengths <- bdp_interval_lengths(time, breaks)
    active <- lengths > 0
    log_variance <- bdp_log_sum_exp(c(
      log(current$shape[active]) + 2 * (log(lengths[active]) - log(current$rate[active])),
      log(historical$shape[active]) + 2 * (log(lengths[active]) - log(historical$rate[active]))))
    Z <- abs(exp(logH - log_variance / 2) - exp(logH0 - log_variance / 2))
  }
  if (method == "fixed") {
    # Split numerical ties equally between the two directions.
    q <- mean(contrast > 0) + 0.5 * mean(contrast == 0)
    return(2 * min(q, 1 - q))
  }
  2 * pnorm(Z, lower.tail = FALSE)
}

# The augmented treatment effect still uses empirical variances across draws,
# since draw-specific borrowing produces a mixture rather than one Gamma law.
bdp_survival_effect <- function(treatment, control) {
  log_hazard <- function(arm) {
    if (!is.null(arm$posterior_log_hazard)) return(arm$posterior_log_hazard)
    log(arm$posterior_hazard) # Fits saved before log hazards were retained.
  }
  contrast <- log_hazard(treatment) - log_hazard(control)
  variance <- apply(contrast, 2, var)
  if (any(!is.finite(variance) | variance <= 0)) {
    stop("Cannot estimate the log hazard ratio: interval variances must be positive and finite.", call. = FALSE)
  }
  weights <- min(variance) / variance
  contrast %*% (weights / sum(weights))
}
