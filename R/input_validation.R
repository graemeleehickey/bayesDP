# Common checks for supplied summaries and borrowing settings.
bdp_numeric_vector <- function(x) {
  is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x))
}

bdp_summary_group <- function(values, group, normal = FALSE, required = FALSE) {
  missing <- vapply(values, is.null, logical(1))
  if (!required && all(missing)) return(FALSE)
  if (any(missing)) {
    for (name in names(values)[missing]) cat(name, " missing\n", sep = "")
    if (group == "Historical treatment") stop("Historical treatment incomplete.", call. = FALSE)
    stop(group, " not provided/incomplete.", call. = FALSE)
  }
  for (name in names(values)) {
    value <- values[[name]]
    if (!bdp_numeric_vector(value) || length(value) != 1L) {
      stop(name, " must be a finite numeric scalar.", call. = FALSE)
    }
  }
  n <- values[[length(values)]]
  minimum <- if (normal) 2 else 1
  if (n < minimum || n != floor(n)) {
    stop(names(values)[length(values)], " must be an integer of at least ", minimum, ".", call. = FALSE)
  }
  if (normal) {
    if (values[[2L]] <= 0) stop(names(values)[2L], " must be positive.", call. = FALSE)
  } else {
    events <- values[[1L]]
    if (events < 0 || events > n || events != floor(events)) {
      stop(names(values)[1L], " must be an integer between 0 and ", names(values)[2L], ".", call. = FALSE)
    }
  }
  TRUE
}

bdp_summary_settings <- function(method, discount_function, alpha_max, fix_alpha,
                                 compare, number_mcmc, weibull_shape, weibull_scale) {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("fixed", "mc")) {
    stop("method must be 'fixed' or 'mc'.", call. = FALSE)
  }
  if (!is.character(discount_function) || length(discount_function) != 1L ||
      is.na(discount_function) ||
      !discount_function %in% c("identity", "weibull", "scaledweibull")) {
    stop("discount_function input incorrectly.", call. = FALSE)
  }
  for (name in c("fix_alpha", "compare")) {
    value <- get(name)
    if (!is.logical(value) || !is.null(dim(value)) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!bdp_numeric_vector(number_mcmc) || length(number_mcmc) != 1L ||
      number_mcmc < 2 || number_mcmc != floor(number_mcmc)) {
    stop("number_mcmc must be an integer of at least 2.", call. = FALSE)
  }
  if (!bdp_numeric_vector(alpha_max) || !length(alpha_max) %in% c(1L, 2L) ||
      any(alpha_max < 0 | alpha_max > 1)) {
    stop("alpha_max must contain one or two finite values in [0, 1].", call. = FALSE)
  }
  if (discount_function != "identity" && !fix_alpha) {
    for (name in c("weibull_shape", "weibull_scale")) {
      value <- get(name)
      if (!bdp_numeric_vector(value) || !length(value) %in% c(1L, 2L) || any(value <= 0)) {
        stop(name, " must contain one or two positive, finite values.", call. = FALSE)
      }
    }
  }
  invisible(NULL)
}

# Evaluate the response once and keep each participant aligned with treatment.
bdp_survival_data <- function(formula, data, source, historical) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  response <- stats::model.response(mf)
  if (!survival::is.Surv(response)) {
    stop("Data input incorrectly. ", source, " response must be a Surv object.", call. = FALSE)
  }
  if (!identical(attr(response, "type"), "right") || ncol(response) != 2L) {
    stop(source, " data must use a right-censored Surv(time, event) response.", call. = FALSE)
  }
  if (!nrow(mf) || any(!is.finite(response[, 1L])) || any(response[, 1L] < 0) ||
      anyNA(response[, 2L]) || any(!response[, 2L] %in% c(0, 1))) {
    stop(source, " survival response must have finite, non-negative times and non-missing binary event indicators.", call. = FALSE)
  }
  treatment <- if ("treatment" %in% names(mf)) mf$treatment else rep.int(1, nrow(mf))
  if (!(is.numeric(treatment) || is.logical(treatment)) || is.complex(treatment) ||
      !is.null(dim(treatment)) || anyNA(treatment) || any(!treatment %in% c(0, 1))) {
    stop(source, " treatment must contain only numeric 0/1 or logical values, with no missing values.", call. = FALSE)
  }
  data.frame(time = response[, 1L], status = response[, 2L],
    treatment = as.integer(treatment), historical = historical)
}
