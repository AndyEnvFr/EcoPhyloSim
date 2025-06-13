#' Get circular neighborhood coordinate offsets
#'
#' @title Compute toroidal neighbor shifts
#' @param radius Integer radius defining the circular neighborhood
#' @return A data frame with columns \code{dx} and \code{dy} indicating x and y offsets
#' @description Computes relative coordinate offsets that define a circular neighborhood
#' of a given radius around a central cell, excluding the center (0, 0). Used in conspecific
#' neighborhood analysis on torus-expanded matrices.
#' 
#' @examples
#' getCircularOffsets(1)
#' getCircularOffsets(3)
#' 
#' @export
getCircularOffsets <- function(radius) {
  offsets <- data.frame()
  
  for (dx in -radius:radius) {
    y_lims <- floor(sqrt(radius^2 - dx^2))
    for (dy in -y_lims:y_lims) {
      if (!(dx == 0 && dy == 0)) {
        offsets <- rbind(offsets, data.frame(dx = dx, dy = dy))
      }
    }
  }
  
  return(offsets)
}
