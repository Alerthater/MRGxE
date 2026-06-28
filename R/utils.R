# ---------------------------------------------------------------------------
# Common utility operators and helpers
# ---------------------------------------------------------------------------

#' Default value operator
#'
#' Returns `b` if `a` is NULL, otherwise returns `a`.
#' Analogous to `??` in other languages.
#'
#' @param a Value to check
#' @param b Default value
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
