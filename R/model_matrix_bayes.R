#' @importFrom stats contrasts<- model.frame model.matrix.default terms
#'   terms.formula
#' @keywords internal
#' @noRd
model.matrixBayes <- function(object, data = environment(object), contrasts.arg = NULL,
                              xlev = NULL, keep.order = FALSE, drop.baseline = FALSE, ...) {
  t <- if (missing(data)) {
    terms(object)
  } else {
    terms.formula(object, data = data, keep.order = keep.order)
  }

  attr(t, "intercept") <- attr(object, "intercept")
  if (is.null(attr(data, "terms"))) {
    data <- model.frame(object, data, xlev = xlev)
  } else {
    reorder <- match(sapply(attr(t, "variables"), deparse, width.cutoff = 500)[-1], names(data))
    if (anyNA(reorder)) {
      stop("model frame and formula mismatch in model.matrix()")
    }
    if (!identical(reorder, seq_len(ncol(data)))) {
      data <- data[, reorder, drop = FALSE]
    }
  }

  int <- attr(t, "response")
  if (length(data)) { # otherwise no rhs terms, so skip all this
    if (drop.baseline) {
      contr.funs <- as.character(getOption("contrasts"))
    } else {
      contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
    }

    namD <- names(data)

    ## turn  character columns into factors
    for (i in namD) {
      if (is.character(data[[i]])) {
        data[[i]] <- factor(data[[i]])
        warning(gettextf("variable '%s' converted to a factor", i), domain = NA)
      }
    }
    isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
    isF[int] <- FALSE
    isOF <- vapply(data, is.ordered, NA)
    for (nn in namD[isF]) { # drop response
      if (is.null(attr(data[[nn]], "contrasts"))) {
        contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
      }
    }

    ## it might be safer to have numerical contrasts:
    ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
    if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
      if (is.null(namC <- names(contrasts.arg))) {
        stop("invalid 'contrasts.arg' argument")
      }
      for (nn in namC) {
        if (is.na(ni <- match(nn, namD))) {
          warning(gettextf("variable '%s' is absent, its contrast will be ignored", nn), domain = NA)
        }
        else {
          ca <- contrasts.arg[[nn]]
          if (is.matrix(ca)) {
            contrasts(data[[ni]], ncol(ca)) <- ca
          }
          else {
            contrasts(data[[ni]]) <- contrasts.arg[[nn]]
          }
        }
      }
    }
  } else {
    isF <- FALSE
    data <- data.frame(x = rep(0, nrow(data)))
  }
  ans <- model.matrix.default(object = t, data = data)
  cons <- if (any(isF)) {
    lapply(data[isF], function(x) attr(x, "contrasts"))
  } else {
    NULL
  }
  attr(ans, "contrasts") <- cons
  ans
}
