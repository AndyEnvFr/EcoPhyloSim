#' Expand matrix with toroidal border
#'
#' @title Internal torus extension
#' @description Adds toroidal border of radius \code{r} to a 2D matrix.
#' The resulting matrix has \code{dim + 2*r} rows and columns.
#'
#' @param matrix A 2D numeric or integer matrix
#' @param max_neighborhood_radius Radius of the toroidal wrap
#' @return A matrix enlarged by 2*r in both dimensions
#' @keywords internal
torus_in <- function(matrix, max_neighborhood_radius) {
  r <- max_neighborhood_radius
  
  mrow <- nrow(matrix)
  mcol <- ncol(matrix)
  mbig <- matrix(NA, nrow = mrow + 2 * r, ncol = mcol + 2 * r)
  
  str_r <- r + 1
  end_r <- nrow(mbig) - r
  str_c <- r + 1
  end_c <- ncol(mbig) - r
  
  mbig[str_r:end_r, str_c:end_c] <- matrix
  
  # wrap edges
  mbig[1:r, str_c:end_c] <- matrix[(mrow - r + 1):mrow, ]            # upper
  mbig[(end_r + 1):(end_r + r), str_c:end_c] <- matrix[1:r, ]        # lower
  mbig[str_r:end_r, 1:r] <- matrix[, (mcol - r + 1):mcol]            # left
  mbig[str_r:end_r, (end_c + 1):(end_c + r)] <- matrix[, 1:r]        # right
  
  # corners
  mbig[1:r, 1:r] <- matrix[(mrow - r + 1):mrow, (mcol - r + 1):mcol]             # upper-left
  mbig[1:r, (end_c + 1):(end_c + r)] <- matrix[(mrow - r + 1):mrow, 1:r]         # upper-right
  mbig[(end_r + 1):(end_r + r), 1:r] <- matrix[1:r, (mcol - r + 1):mcol]         # lower-left
  mbig[(end_r + 1):(end_r + r), (end_c + 1):(end_c + r)] <- matrix[1:r, 1:r]     # lower-right
  
  return(mbig)
}
