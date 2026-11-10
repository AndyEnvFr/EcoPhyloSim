#' @title Toroidal matrix expansion (helper)
#' @param mat Input matrix
#' @param r Radius to expand
#' @return Matrix with \code{r} rows and columns added on each side, filled
#'   with the values from the opposite side (periodic boundaries)
#' @keywords internal
getTorusMatrix <- function(mat, r) {
  rows <- (seq(-r, nrow(mat) + r - 1) %% nrow(mat)) + 1
  cols <- (seq(-r, ncol(mat) + r - 1) %% ncol(mat)) + 1
  mat[rows, cols, drop = FALSE]
}
