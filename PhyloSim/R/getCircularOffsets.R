#' Get circular neighborhood coordinate offsets
#'
#' @title Compute coordinate offsets for circular neighborhoods
#' @param radius Numeric value specifying the radius of the circular neighborhood. 
#'   Must be a non-negative integer. The actual neighborhood will include all integer 
#'   coordinate pairs (dx, dy) where sqrt(dx² + dy²) ≤ radius.
#' @return A data frame with two integer columns:
#'   \describe{
#'     \item{dx}{X-coordinate offsets relative to center cell}
#'     \item{dy}{Y-coordinate offsets relative to center cell}
#'   }
#'   The center coordinate (0, 0) is excluded. Rows are ordered by x-offset, 
#'   then by y-offset within each x-value.
#' @description
#' Calculates coordinate offsets for all cells within a circular neighborhood around a central cell,
#' with a given radius matching the simulation's density cutoff distance.
#'
#' The simulation's densitycut parameter limits density-dependent effects to cells within this circular radius.
#' This function identifies exactly those neighboring cells that would be affected by density effects
#' in the simulation, generating all integer coordinate pairs within the specified radius (excluding
#' the center cell itself).
#'
#' The radius definition is consistent with the simulation's C++ implementation.
#' 
#' @details
#' The algorithm iterates through all possible x-offsets from -radius to +radius,
#' then for each x-offset, calculates the maximum y-offset that still falls within
#' the circular boundary using the Pythagorean theorem: y_max = sqrt(radius² - x²).
#' 
#' @example /inst/examples/getCircularOffsets-help.R
#' 
#' @seealso 
#' \code{\link{expand.grid}} for generating rectangular neighborhoods
#' 
#' @export
getCircularOffsets <- function(radius) {
  if (length(radius) != 1 || is.na(radius) || radius < 0 || radius != round(radius)) stop("radius must be a single non-negative integer")
  offsets <- data.frame()
  
  for (dx in -radius:radius) {
    y_lims <- floor(sqrt(radius^2 - dx^2)) # based on cpp code
    for (dy in -y_lims:y_lims) {
      if (!(dx == 0 && dy == 0)) {
        offsets <- rbind(offsets, data.frame(dx = dx, dy = dy))
      }
    }
  }
  
  return(offsets)
}
