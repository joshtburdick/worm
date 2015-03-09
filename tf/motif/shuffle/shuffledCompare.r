# Compares motifs with shuffled versions of themselves.

library("MotIV")

source("git/utils.r")

load("git/tf/motif/meme.format.pwm.Rdata")

output.dir = "git/tf/motif/shuffle/shuffledCompare/"

# get motif names to run this on
# XXX
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}
hughes.motif = sub(".Rdata", "",
  list.files("git/tf/motif/count/upstreamMotifCount/hughes_20141202/"))
motif.list = c(orig.motif.list, hughes.motif)

system(paste0("mkdir -p ", output.dir))

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
#    a = NULL    # trying to prevent a memory leak
  }

  r
}

# Gets shuffled motif distributions for many motifs.
# Args:
#   motif.list - list of motifs to get distributions for
#   num.shuffles - number of times to shuffle these
write.shuffled.motif.dist = function(motif.names, num.shuffles=1000) {
  for(m in motif.names) {
    write.status(m)
    if (m %in% names(meme.format.pwm)) {
      shuffled.motif.dist = shuffled.motif.dist.1(meme.format.pwm[[m]], num.shuffles)
      save(shuffled.motif.dist, file=paste0(output.dir, "/", m, ".Rdata"))
    }
  }
}

# Runs this on one "chunk" of motifs.
run.chunk = function(low, hi) {

  write.shuffled.motif.dist(motif.list[low:hi])
}

# XXX amazingly lame parallelism: run each of these "by hand"
# run.chunk(1, 499)
# run.chunk(500, 799)
# run.chunk(800, 1099)
# run.chunk(1100, 1299)
# run.chunk(1300, 1497)

