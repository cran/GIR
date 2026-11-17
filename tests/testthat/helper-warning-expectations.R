# Check the entire ordered warning set, not just the presence of one warning.
# With no expected messages this also verifies that a successful call is quiet.
expect_gir_warnings <- function(expr, expected = character()) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(length(messages) == length(expected),
              info = paste(messages, collapse = "\n"))
  for (i in seq_len(min(length(messages), length(expected)))) {
    expect_match(messages[i], expected[i], fixed = TRUE)
  }
  invisible(value)
}
