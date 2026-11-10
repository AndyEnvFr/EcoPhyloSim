# Small neighborhood (4 direct neighbours, no diagonals)
getCircularOffsets(1)
#>   dx dy
#> 1 -1  0
#> 2  0 -1
#> 3  0  1
#> 4  1  0

# Larger circular neighborhood
offsets_r3 <- getCircularOffsets(3)
nrow(offsets_r3)  # Number of cells in neighborhood
#> [1] 28

# Visualize the circular pattern
plot(offsets_r3$dx, offsets_r3$dy, asp = 1, pch = 16,
     main = "Circular neighborhood (radius = 3)",
     xlab = "dx", ylab = "dy")
points(0, 0, pch = 4, col = "red", cex = 2)  # Center point
