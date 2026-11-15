#' @title bdpnormal Object Plot
#' @description \code{plot} method for class \code{bdpnormal}.
#'
#' @param x object of class \code{bdpnormal}. The result of a call to the
#'   \code{\link{bdpnormal}} function.
#' @param type character. Optional. Select plot type to print. Supports the
#'   following: "discount" gives the discount function; "posteriors" gives the
#'   posterior distributions by information source; "density" gives the
#'   augmented arm distributions; and "effect" gives the treatment-minus-control
#'   difference for a two-arm analysis. Leave NULL to display the first three
#'   plots in sequence.
#' @param print logical. Optional. Produce a plot (\code{print = TRUE}; default)
#'   or return a ggplot object (\code{print = FALSE}). When \code{print =
#'   FALSE}, it is possible to use \code{ggplot2} syntax to modify the plot
#'   appearance.
#' @details The \code{plot} method for \code{bdpnormal} returns up to three
#'   plots: (1) posterior density curves by information source; (2) augmented
#'   arm densities; and (3) the discount function. A call to \code{plot} that omits
#'   the \code{type} input returns one plot at a time and prompts the user to
#'   click the plot or press return to proceed to the next plot. Otherwise, the
#'   user can specify a plot \code{type} to display the requested plot.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes ggtitle ylim guides guide_legend theme
#' @importFrom ggplot2 .data
#'   element_blank
#' @importFrom graphics par
#' @importFrom stats sd density is.empty.model median model.offset
#'   model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("plot", signature(x = "bdpnormal"), function(x, type = NULL, print = TRUE) {
  if (identical(type, "effect")) {
    if (!x$args1$arm2) stop("An effect plot requires treatment and control data.")
    p <- bdp_effect_plot(x$posterior_treatment$posterior_mu - x$posterior_control$posterior_mu, "Mean difference (treatment - control)")
    if (print) return(base::print(p))
    return(p)
  }
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  discount_function <- x$args1$discount_function
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  N_t <- x$args1$N_t
  N_c <- x$args1$N_c
  arm2 <- x$args1$arm2
  method <- x$args1$method

  information_sources <- NULL

  ### Create data frames for plotting via ggplot2
  D1 <- D2 <- D3 <- D5 <- D6 <- NULL

  if (arm2) {
    dens1 <- density(posterior_control$posterior_mu,
                     adjust = 0.5
    )
    D1 <- data.frame(
      information_sources = "Augmented posterior",
      group = "Control",
      x = dens1$x,
      y = dens1$y
    )

    if (!is.null(posterior_control$posterior_flat_mu)) {
      dens2 <- density(posterior_control$posterior_flat_mu,
                       adjust = 0.5
      )
      D2 <- data.frame(
        information_sources = "Current-only posterior",
        group = "Control",
        x = dens2$x,
        y = dens2$y
      )
    }

    if (!is.null(posterior_control$prior_mu)) {
      dens3 <- density(posterior_control$prior_mu,
                       adjust = 0.5
      )
      D3 <- data.frame(
        information_sources = "Historical-only posterior",
        group = "Control",
        x = dens3$x,
        y = dens3$y
      )
    }
  }

  dens4 <- density(posterior_treatment$posterior_mu,
                   adjust = 0.5
  )
  D4 <- data.frame(
    information_sources = "Augmented posterior",
    group = "Treatment",
    x = dens4$x,
    y = dens4$y
  )

  if (!is.null(posterior_treatment$posterior_flat_mu)) {
    dens5 <- density(posterior_treatment$posterior_flat_mu,
                     adjust = 0.5
    )
    D5 <- data.frame(
      information_sources = "Current-only posterior",
      group = "Treatment",
      x = dens5$x,
      y = dens5$y
    )
  }

  if (!is.null(posterior_treatment$prior_mu)) {
    dens6 <- density(posterior_treatment$prior_mu,
                     adjust = 0.5
    )
    D6 <- data.frame(
      information_sources = "Historical-only posterior",
      group = "Treatment",
      x = dens6$x,
      y = dens6$y
    )
  }

  D <- rbind(D1, D2, D3, D4, D5, D6)

  D$information_sources <- factor(D$information_sources,
                                  levels = c("Augmented posterior", "Current-only posterior", "Historical-only posterior")
  )
  D$group <- factor(D$group, levels = c("Treatment", "Control"))

  ##############################################################################
  ### Posterior Type Plots
  ##############################################################################

  post_typeplot <- ggplot(D, aes(x = .data[["x"]], y = .data[["y"]])) +
    bdp_source_scales() +
    geom_line(
      linewidth = 2,
      aes(
        colour = .data[["information_sources"]],
        lty = .data[["information_sources"]]
      )
    ) +
    facet_wrap(~group, ncol = 1, scales = "free") +
    ylab("Density (PDF)") +
    xlab("Population mean") +
    theme_bw() +
    ggtitle("Posterior distributions by information source") +
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.title = element_blank())

  ##############################################################################
  ### Density Plots
  ##############################################################################

  densityplot <- ggplot(
    subset(D, information_sources == "Augmented posterior"),
    aes(x = .data[["x"]], y = .data[["y"]])
  ) +
    geom_line(linewidth = 2, aes(colour = .data[["group"]])) +
    ylab("Density (PDF)") +
    xlab("Population mean") +
    theme_bw() +
    ggtitle("Augmented posterior distributions") +
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.title = element_blank())

  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################

  discountfun_plot <- bdp_discount_plot(x, !is.null(N_t) && !is.null(N0_t),
                                        !is.null(N_c) && !is.null(N0_c))

  if (print) {
    if (is.null(type)) {
      op <- par(ask = interactive())
      plot(post_typeplot)
      plot(densityplot)
      if (!is.null(discountfun_plot)) {
        plot(discountfun_plot)
      }
      par(op)
    } else if (type == "discount") {
      if (!is.null(discountfun_plot)) {
        plot(discountfun_plot)
      }
    } else if (type == "posteriors") {
      plot(post_typeplot)
    } else if (type == "density") {
      plot(densityplot)
    }
  } else {
    if (is.null(type)) {
      out <- list(
        posteriors = post_typeplot,
        density = densityplot
      )
      if (!is.null(discountfun_plot)) {
        out <- c(out, discount = list(discountfun_plot))
        out
      } else {
        out
      }
    } else if (type == "discount") {
      if (!is.null(discountfun_plot)) {
        discountfun_plot
      }
    } else if (type == "posteriors") {
      post_typeplot
    } else if (type == "density") {
      densityplot
    }
  }
})


#' @title bdpbinomial Object Plot
#' @description \code{plot} method for class \code{bdpbinomial}.
#'
#' @param x object of class \code{bdpbinomial}. The result of a call to the
#'   \code{\link{bdpbinomial}} function.
#' @param type character. Optional. Select plot type to print. Supports the
#'   following: "discount" gives the discount function; "posteriors" gives the
#'   posterior distributions by information source; "density" gives the
#'   augmented arm distributions; and "effect" gives the treatment-minus-control
#'   difference for a two-arm analysis. Leave NULL to display the first three
#'   plots in sequence.
#' @param print logical. Optional. Produce a plot (\code{print = TRUE}; default)
#'   or return a ggplot object (\code{print = FALSE}). When \code{print =
#'   FALSE}, it is possible to use \code{ggplot2} syntax to modify the plot
#'   appearance.
#'
#' @details The \code{plot} method for \code{bdpbinomial} returns up to three
#'   plots: (1) posterior density curves by information source; (2) augmented
#'   arm densities; and (3) the discount function. A call to \code{plot} that omits
#'   the \code{type} input returns one plot at a time and prompts the user to
#'   click the plot or press return to proceed to the next plot. Otherwise, the
#'   user can specify a plot \code{type} to display the requested plot.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes ggtitle ylim guides guide_legend theme
#' @importFrom ggplot2 .data
#'   element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("plot", signature(x = "bdpbinomial"), function(x, type = NULL, print = TRUE) {
  if (identical(type, "effect")) {
    if (!x$args1$arm2) stop("An effect plot requires treatment and control data.")
    p <- bdp_effect_plot(x$posterior_treatment$posterior - x$posterior_control$posterior, "Risk difference (treatment - control)")
    if (print) return(base::print(p))
    return(p)
  }
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  discount_function <- x$args1$discount_function
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  N_t <- x$args1$N_t
  N_c <- x$args1$N_c
  arm2 <- x$args1$arm2
  method <- x$args1$method

  information_sources <- NULL

  ### Create data frames for plotting via ggplot2
  D1 <- D2 <- D3 <- D5 <- D6 <- NULL

  if (arm2) {
    dens1 <- density(posterior_control$posterior,
                     adjust = 0.5
    )
    D1 <- data.frame(
      information_sources = "Augmented posterior",
      group = "Control",
      x = dens1$x,
      y = dens1$y
    )

    if (!is.null(posterior_control$posterior_flat)) {
      dens2 <- density(posterior_control$posterior_flat, adjust = 0.5)
      D2 <- data.frame(
        information_sources = "Current-only posterior",
        group = "Control",
        x = dens2$x,
        y = dens2$y
      )
    }

    if (!is.null(posterior_control$prior)) {
      dens3 <- density(posterior_control$prior, adjust = 0.5)
      D3 <- data.frame(
        information_sources = "Historical-only posterior",
        group = "Control",
        x = dens3$x,
        y = dens3$y
      )
    }
  }

  dens4 <- density(posterior_treatment$posterior, adjust = 0.5)
  D4 <- data.frame(
    information_sources = "Augmented posterior",
    group = "Treatment",
    x = dens4$x,
    y = dens4$y
  )

  if (!is.null(posterior_treatment$posterior_flat)) {
    dens5 <- density(posterior_treatment$posterior_flat, adjust = 0.5)
    D5 <- data.frame(
      information_sources = "Current-only posterior",
      group = "Treatment",
      x = dens5$x,
      y = dens5$y
    )
  }

  if (!is.null(posterior_treatment$prior)) {
    dens6 <- density(posterior_treatment$prior, adjust = 0.5)
    D6 <- data.frame(
      information_sources = "Historical-only posterior",
      group = "Treatment",
      x = dens6$x,
      y = dens6$y
    )
  }

  D <- rbind(D1, D2, D3, D4, D5, D6)

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Augmented posterior", "Current-only posterior", "Historical-only posterior"))
  )
  D$group <- factor(D$group, levels = c("Treatment", "Control"))

  ##############################################################################
  ### Posterior Type Plots
  ##############################################################################

  post_typeplot <- ggplot(D, aes(x = .data[["x"]], y = .data[["y"]])) +
    bdp_source_scales() +
    geom_line(linewidth = 2, aes(color = .data[["information_sources"]], lty = .data[["information_sources"]])) +
    theme_bw() +
    facet_wrap(~group, ncol = 1, scales = "free") +
    ylab("Density (PDF)") +
    xlab("Event probability") +
    ggtitle("Posterior distributions by information source") +
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.title = element_blank())

  ##############################################################################
  ### Density Plots
  ##############################################################################

  densityplot <- ggplot(
    subset(D, information_sources == "Augmented posterior"),
    aes(x = .data[["x"]], y = .data[["y"]])
  ) +
    geom_line(linewidth = 2, aes(color = .data[["group"]])) +
    ylab("Density (PDF)") +
    xlab("Event probability") +
    theme_bw() +
    ggtitle("Augmented posterior distributions") +
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.title = element_blank())

  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################

  discountfun_plot <- bdp_discount_plot(x, !is.null(N_t) && !is.null(N0_t),
                                        !is.null(N_c) && !is.null(N0_c))

  if (print) {
    if (is.null(type)) {
      op <- par(ask = interactive())
      plot(post_typeplot)
      plot(densityplot)
      if (!is.null(discountfun_plot)) {
        plot(discountfun_plot)
      }
      par(op)
    } else if (type == "discount") {
      if (!is.null(discountfun_plot)) {
        plot(discountfun_plot)
      }
    } else if (type == "posteriors") {
      plot(post_typeplot)
    } else if (type == "density") {
      plot(densityplot)
    }
  } else {
    if (is.null(type)) {
      out <- list(
        posteriors = post_typeplot,
        density = densityplot
      )
      if (!is.null(discountfun_plot)) {
        out <- c(out, discount = list(discountfun_plot))
        out
      } else {
        out
      }
    } else if (type == "discount") {
      if (!is.null(discountfun_plot)) {
        discountfun_plot
      }
    } else if (type == "posteriors") {
      post_typeplot
    } else if (type == "density") {
      densityplot
    }
  }
})


#' @title bdpsurvival Object Plot
#' @description \code{plot} method for class \code{bdpsurvival}.
#'
#' @param x object of class \code{bdpsurvival}. The result of a call to the
#'   \code{\link{bdpsurvival}} function.
#' @param type character. Optional. Select plot type to print. Supports the
#'   following: "discount" gives the discount function and "survival" gives the
#'   survival curves. Leave NULL to display all plots in sequence.
#' @param print logical. Optional. Produce a plot (\code{print = TRUE}; default)
#'   or return a ggplot object (\code{print = FALSE}). When \code{print =
#'   FALSE}, it is possible to use \code{ggplot2} syntax to modify the plot
#'   appearance.
#'
#' @details The \code{plot} method for \code{bdpsurvival} returns up to two
#'   plots: (1) posterior survival curves and (2) the discount function. A call
#'   to \code{plot} that omits the \code{type} input returns one plot at a time
#'   and prompts the user to click the plot or press return to proceed to the
#'   next plot. Otherwise, the user can specify a plot \code{type} to display
#'   the requested plot.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes ggtitle ylim guides guide_legend theme
#' @importFrom ggplot2 .data
#'   element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("plot", signature(x = "bdpsurvival"), function(x, type = NULL, print = TRUE) {
  args1 <- x$args1
  discount_function <- x$args1$discount_function
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  data <- args1$data
  breaks <- args1$breaks
  arm2 <- args1$arm2
  method <- args1$method

  ##############################################################################
  ### Survival curve(s)
  ### - Only computed for one-arm trial
  ##############################################################################

  ### Organize data for current treatment
  time_t <- sort(unique(args1$S_t$time))
  survival_times_posterior_flat <- ppexp_times(time_t,
                                          posterior_treatment$posterior_flat_hazard,
                                          cuts = c(0, breaks)
  )
  survival_median_posterior_flat <- 1 - sapply(survival_times_posterior_flat, median)

  D1 <- data.frame(
    source = "Current-only posterior",
    group = "Treatment",
    x = time_t,
    y = survival_median_posterior_flat
  )

  ### Organize data for historical treatment
  if (!is.null(args1$S0_t)) {
    time0_t <- sort(unique(args1$S0_t$time))
    survival_times_prior <- ppexp_times(time0_t,
                                   posterior_treatment$prior_hazard,
                                   cuts = c(0, breaks)
    )
    survival_median_prior <- 1 - sapply(survival_times_prior, median)

    D2 <- data.frame(
      source = "Historical-only posterior",
      group = "Treatment",
      x = time0_t,
      y = survival_median_prior
    )
  } else {
    D2 <- NULL
  }

  ### Organize data for treatment posterior
  survival_times_posterior <- ppexp_times(time_t, posterior_treatment$posterior_hazard, cuts = c(0, breaks))
  survival_median_posterior <- 1 - sapply(survival_times_posterior, median)

  D3 <- data.frame(
    source = "Augmented posterior",
    group = "Treatment",
    x = time_t,
    y = survival_median_posterior
  )

  ### Organize data for current control
  if (!is.null(args1$S_c)) {
    time_c <- sort(unique(args1$S_c$time))
    survival_times_posterior_flat <- ppexp_times(time_c,
                                            posterior_control$posterior_flat_hazard,
                                            cuts = c(0, breaks)
    )
    survival_median_posterior_flat <- 1 - sapply(survival_times_posterior_flat, median)

    D4 <- data.frame(
      source = "Current-only posterior",
      group = "Control",
      x = time_c,
      y = survival_median_posterior_flat
    )
  } else {
    D4 <- NULL
  }

  ### Organize data for historical control
  if (!is.null(args1$S0_c)) {
    time0_c <- sort(unique(args1$S0_c$time))
    survival_times_prior <- ppexp_times(time0_c,
                                   posterior_control$prior_hazard,
                                   cuts = c(0, breaks)
    )
    survival_median_prior <- 1 - sapply(survival_times_prior, median)

    D5 <- data.frame(
      source = "Historical-only posterior",
      group = "Control",
      x = time0_c,
      y = survival_median_prior
    )
  } else {
    D5 <- NULL
  }

  ### Organize data for control posterior
  if (!is.null(args1$S_c) & !is.null(args1$S0_c)) {
    survival_times_posterior <- ppexp_times(time_c, posterior_control$posterior_hazard, cuts = c(0, breaks))
    survival_median_posterior <- 1 - sapply(survival_times_posterior, median)

    D6 <- data.frame(
      source = "Augmented posterior",
      group = "Control",
      x = time_c,
      y = survival_median_posterior
    )
  } else if (is.null(args1$S_c) & !is.null(args1$S0_c)) {
    survival_times_posterior <- ppexp_times(time0_c, posterior_control$posterior_hazard, cuts = c(0, breaks))
    survival_median_posterior <- 1 - sapply(survival_times_posterior, median)

    D6 <- data.frame(
      source = "Augmented posterior",
      group = "Control",
      x = time0_c,
      y = survival_median_posterior
    )
  } else if (!is.null(args1$S_c) & is.null(args1$S0_c)) {
    survival_times_posterior <- ppexp_times(time_c, posterior_control$posterior_hazard, cuts = c(0, breaks))
    survival_median_posterior <- 1 - sapply(survival_times_posterior, median)

    D6 <- data.frame(
      source = "Augmented posterior",
      group = "Control",
      x = time_c,
      y = survival_median_posterior
    )
  } else {
    D6 <- NULL
  }

  D <- rbind(D4, D5, D6, D1, D2, D3)
  D$group <- factor(D$group, levels = c("Treatment", "Control"))

  ### Plot survival curve
  survival_curves <- ggplot(D, aes(x = .data[["x"]], y = .data[["y"]])) +
    bdp_source_scales() +
    geom_line(linewidth = 1.4, aes(color = .data[["source"]], lty = .data[["source"]])) +
    facet_wrap(~group, ncol = 1, scales = "free") +
    ylab("Survival probability") +
    xlab("Time") +
    theme_bw() +
    ggtitle("Posterior median survival curves") +
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.title = element_blank())

  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################

  discountfun_plot <- bdp_discount_plot(x, !is.null(args1$S_t) && !is.null(args1$S0_t),
                                        !is.null(args1$S_c) && !is.null(args1$S0_c))

  if (print) {
    if (is.null(type)) {
      op <- par(ask = interactive())
      plot(survival_curves)
      if (!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))) {
        plot(discountfun_plot)
      }
      par(op)
    } else if (type == "discount") {
      if (!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))) {
        plot(discountfun_plot)
      }
    } else if (type == "survival") {
      plot(survival_curves)
    }
  } else {
    if (is.null(type)) {
      out <- list(survival = survival_curves)
      if (!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))) {
        out <- c(out, discount = list(discountfun_plot))
        out
      } else {
        out
      }
    } else if (type == "discount") {
      if (!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))) {
        discountfun_plot
      }
    } else if (type == "survival") {
      survival_curves
    }
  }
})

# Consistent colours and line types identify information sources across outcomes.
bdp_source_scales <- function() {
  labels <- c("Augmented posterior", "Current-only posterior", "Historical-only posterior")
  list(
    ggplot2::scale_colour_manual(values = stats::setNames(c("#0072B2", "#666666", "#D55E00"), labels)),
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed", "dotdash"), labels))
  )
}

# Plot the weight actually applied, including alpha_max and prescribed weights.
bdp_discount_plot <- function(x, treatment_present, control_present) {
  present <- c(treatment_present, control_present)
  if (!any(present)) return(NULL)
  curves <- comparisons <- weights <- list()
  grid <- seq(0, 1, length.out = 101)
  for (i in which(present)) {
    arm <- c("Treatment", "Control")[i]
    posterior <- list(x$posterior_treatment, x$posterior_control)[[i]]
    maximum <- rep(x$args1$alpha_max, length.out = 2)[i]
    values <- if (x$args1$fix_alpha) {
      rep(maximum, length(grid))
    } else {
      alpha_discount(grid,
        discount_function = x$args1$discount_function,
        alpha_max = maximum,
        weibull_shape = rep(x$args1$weibull_shape, length.out = 2)[i],
        weibull_scale = rep(x$args1$weibull_scale, length.out = 2)[i]
      )
    }
    curves[[i]] <- data.frame(group = arm, p = grid, alpha = values)
    comparisons[[i]] <- data.frame(group = arm, p = median(posterior$p_hat))
    weights[[i]] <- data.frame(group = arm, alpha = median(posterior$alpha_discount))
  }
  ggplot(do.call(rbind, curves), aes(x = .data[["p"]], y = .data[["alpha"]])) +
    geom_line(colour = "#0072B2", linewidth = 1) +
    geom_vline(data = do.call(rbind, comparisons), aes(xintercept = .data[["p"]]), lty = 2) +
    geom_hline(data = do.call(rbind, weights), aes(yintercept = .data[["alpha"]]), lty = 2) +
    facet_wrap(~group, ncol = 1) +
    theme_bw() + ylim(0, 1) +
    ggplot2::labs(
      x = "Current-historical comparison, p", y = "Historical weight, alpha",
      title = if (x$args1$fix_alpha) "Prespecified historical weight" else "Historical borrowing",
      subtitle = if (identical(x$args1$method, "mc") && !x$args1$fix_alpha) "Dashed lines show medians across draws" else NULL
    )
}

bdp_effect_plot <- function(draws, label) {
  d <- density(draws, adjust = 0.5)
  ggplot(data.frame(effect = d$x, density = d$y), aes(x = .data[["effect"]], y = .data[["density"]])) +
    geom_line(colour = "#0072B2", linewidth = 1) +
    geom_vline(xintercept = 0, lty = 2, colour = "#666666") +
    theme_bw() + ggplot2::labs(x = label, y = "Density", title = "Posterior treatment effect")
}
