# Check the generated equations using the same KaTeX version as the website.
# Run from the package directory: Rscript .github/scripts/validate-katex.R docs

validate_katex <- function(site) {
  pages <- list.files(site, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
  if (!length(pages)) {
    stop("No HTML pages found in ", site, ". Build the pkgdown site first.", call. = FALSE)
  }

  renderers <- new.env(parent = emptyenv())
  javascript <- tempfile(fileext = ".js")
  on.exit(unlink(javascript), add = TRUE)
  total <- 0L
  failures <- character()

  for (page in pages) {
    html <- xml2::read_html(page)
    equations <- xml2::xml_find_all(
      html, '//*[contains(concat(" ", normalize-space(@class), " "), " math ")]'
    )
    if (!length(equations)) next

    scripts <- xml2::xml_attr(xml2::xml_find_all(html, "//script[@src]"), "src")
    katex_url <- unique(grep(
      "^https://cdn[.]jsdelivr[.]net/npm/katex@[0-9]+[.][0-9]+[.][0-9]+/dist/katex[.]min[.]js$",
      scripts, value = TRUE
    ))
    if (length(katex_url) != 1L) {
      stop("Expected one versioned KaTeX script in ", page, ".", call. = FALSE)
    }

    if (!exists(katex_url, envir = renderers, inherits = FALSE)) {
      utils::download.file(katex_url, javascript, mode = "wb", quiet = TRUE)
      ctx <- V8::v8()
      invisible(ctx$source(javascript))
      invisible(ctx$eval('
        function validateMath(equations) {
          var macros = {};
          var failures = [];
          equations.forEach(function(equation, i) {
            try {
              katex.renderToString(equation.tex, {
                displayMode: equation.display,
                throwOnError: true,
                strict: "error",
                macros: macros
              });
            } catch (error) {
              failures.push("Equation " + (i + 1) + ": " + error.message +
                "\\nLaTeX: " + equation.tex);
            }
          });
          return failures;
        }
      '))
      assign(katex_url, ctx, envir = renderers)
      message("Using KaTeX ", ctx$get("katex.version"), " from ", katex_url)
    }

    ctx <- get(katex_url, envir = renderers, inherits = FALSE)
    formulas <- lapply(equations, function(equation) {
      list(
        tex = xml2::xml_text(equation),
        display = grepl("(^|[[:space:]])display([[:space:]]|$)", xml2::xml_attr(equation, "class"))
      )
    })
    errors <- ctx$call("validateMath", formulas)
    total <- total + length(equations)
    if (length(errors)) {
      failures <- c(failures, paste0(page, ": ", errors))
    }
    message(page, ": ", length(equations), " equations, ", length(errors), " failures")
  }

  if (total == 0L) {
    stop("No equations found. Check the site output and equation markup.", call. = FALSE)
  }
  if (length(failures)) {
    stop(paste(c("KaTeX validation failed:", failures), collapse = "\n"), call. = FALSE)
  }
  message("Strict KaTeX validation passed for all ", total, " equations.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript .github/scripts/validate-katex.R [docs]")
validate_katex(if (length(args)) args[[1L]] else "docs")
