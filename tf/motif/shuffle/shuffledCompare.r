# Compares motifs with shuffled versions of themselves.

library("MotIV")

source("git/utils.r")

load("git/tf/motif/meme.format.pwm.Rdata")

# get motif names to run this on

# XXX
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}
hughes.motif = sub(".Rdata", "",
  list.files("git/tf/motif/count/upstreamMotifCount/hughes_20141202/"))
motif.list = c(orig.motif.list, hughes.motif)

# Shuffles a motif many times, and computes the
# distance between the original and shuffled motif.
# Args:
#   m - the motif to shuffle
#   num.shuffles - number of random shuffles to do
# Returns: vector of distances
shuffled.motif.dist.1 = function(m, num.shuffles=1000) {
  r = rep(NA, num.shuffles)
  motif.length = ncol(m)

  for(iter in 1:num.shuffles) {
    m1 = m[ , sample(motif.length) ]
    a = motifDistances(list(m=m, m1=m1))
    r[iter] = a[1]
    a = NULL    # trying to prevent a memory leak
  }

  r
}

# Gets shuffled motif distributions for many motifs.
# Args:
#   motif.list - list of motifs to get distributions for
#   num.shuffles - number of times to shuffle these
shuffled.motif.dist = function(motif.names, num.shuffles=1000) {
  r = matrix(NA, nrow=length(motif.names), ncol=num.shuffles)
  rownames(r) = motif.names

  for(m in motif.names) {
    write.status(m)
    if (m %in% names(meme.format.pwm))
      r[m,] = shuffled.motif.dist.1(meme.format.pwm[[m]], num.shuffles)
  }

  r
}

shuffled.compare = shuffled.motif.dist(motif.list)

save(shuffled.compare, file="git/tf/motif/shuffle/shuffledCompare.Rdata")

