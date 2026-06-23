# Open a null graphics device for the duration of the test run so that plot
# tests do not write a stray Rplots.pdf into the test directory.
pdf(NULL)
withr::defer(
  {
    # Close any devices opened during testing (the null pdf above, and any
    # left open by a failing test).
    while (!is.null(dev.list())) dev.off()
  },
  teardown_env()
)
