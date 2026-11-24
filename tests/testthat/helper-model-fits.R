# Small actual fits exercise method dispatch and all supported result structures.
make_tidier_fits <- function() {
  withr::local_seed(2981)
  current <- expand.grid(treatment = 0:1, site = factor(c("A", "B")), replicate = 1:12)
  current$x <- rnorm(nrow(current), 5, 1)
  current$Y <- 2 + current$treatment + 0.3 * current$x + rnorm(nrow(current))
  historical <- current
  historical$x <- rnorm(nrow(historical), 6, 1)
  historical$Y <- 2 + historical$treatment + 0.3 * historical$x + rnorm(nrow(historical))
  continuous <- lapply(c("fixed", "mc"), function(method) {
    bdplm(Y ~ treatment + x + site, data = current, data0 = historical,
      method = method, number_mcmc_alpha = 300, number_mcmc_sigmagrid = 100,
      number_mcmc_sigma = 50, number_mcmc_beta = 301)
  })
  # Both outcomes occur in every arm/site group, avoiding separation.
  current$Y <- as.numeric(current$replicate %% 3 == 0)
  historical$Y <- as.numeric(historical$replicate %% 4 == 0)
  logistic <- bdplogit(Y ~ treatment + x + site, data = current, data0 = historical,
    method = "fixed", number_mcmc_alpha = 300, number_mcmc_beta = 400)
  survival_data <- data.frame(time = rep(c(1, 2, 2, 4, 6, 8), 2),
    status = rep(c(1, 0, 1, 0, 1, 1), 2), treatment = rep(0:1, each = 6))
  survival_historical <- transform(survival_data, time = time + 0.5)
  list(lm_fixed = continuous[[1]], lm_mc = continuous[[2]], logit = logistic,
    survival_one = bdpsurvival(Surv(time, status) ~ 1,
      data = subset(survival_data, treatment == 1),
      data0 = subset(survival_historical, treatment == 1), breaks = c(2, 5),
      surv_time = 3, method = "fixed", number_mcmc = 300),
    survival_two = bdpsurvival(Surv(time, status) ~ treatment,
      data = survival_data, data0 = survival_historical, breaks = c(2, 5),
      method = "mc", compare = FALSE, number_mcmc = 300))
}
