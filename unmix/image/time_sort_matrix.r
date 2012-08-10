# Sort matrix to integrate timeseries data.

# cell volume (and weight)
source("git/unmix/image/cell_volume.r")

time.points = c(0,30,60,90,120,140,180,240,270,300,330,360,390,
  420,450,480,540,570,600,630,660,690,720)

# Computes a row of the sorting matrix for a given time range,
# weighted by time and volume.
# Args:
#   cell.weights - the cell weights matrix
#   a, b - start and end time
# Returns: proportion of cells in that "time slice"
time.sort.matrix.row = function(cell.weights, a, b) {

  # compute proportion of each cell's lifetime that is within
  # the range [a,b]
  t1 = pmax(cell.weights$on, a)
  t2 = pmin(cell.weights$off, b)
  t.overlap = ifelse(t2 > t1,
    (t2 - t1) / (cell.weights$off - cell.weights$on),
    0)

  # the weights of cells (by overlap, not volume)
  w.time.weighted = t.overlap    # cell.weights$w * t.overlap
  names(w.time.weighted) = rownames(cell.weights)

  r = w.time.weighted[lin.node.names]
  names(r) = lin.node.names
  r[is.na(r)] = 0
  r
}

time.sort.matrix.unweighted = NULL
for(i in c(2:18)) {
  time.sort.matrix.unweighted = rbind(time.sort.matrix.unweighted,
    time.sort.matrix.row(cell.weights, time.points[i-1], time.points[i]))
}
rownames(time.sort.matrix) = paste("t.", time.points[2:18], sep="")

time.sort.matrix = t( t(time.sort.matrix.unweighted) * cell.weights[lin.node.names,"w"] )
time.sort.matrix[ is.na(time.sort.matrix) ] = 0

save(time.sort.matrix, time.sort.matrix.unweighted, file="git/unmix/image/time_sort_matrix.Rdata")

