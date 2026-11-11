# Offsets (dx, dy) of all cells with dx^2 + dy^2 <= radius^2 around a focal cell, without the
# focal cell; the same neighborhood as for density dependence in Landscape::densityUpdate (Grid.cpp).
# Used by getMortalityData().
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
