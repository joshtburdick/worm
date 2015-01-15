# Somewhat like approx_region, but models posterior as a Dirichlet.
# Also potentially unmixes more than one gene at a time.
# For now, only includes exact constraints.

source("git/unmix/ept/approx_region.r")

# for printing updates
backspace = paste(rep("\b", 70), collapse="")

# Conditions a Dirichlet on some proportions.
# Args:
#   x - the count parameters of the Dirichlet (_not_ the natural parameters)
#   A, b - these give the constraint
# Returns: x counts, conditional on Ax \propto b
# FIXME: check that this calculation is right
# ??? eventually vectorize this?
dirichlet.cond = function(x, A, b) {
  A = as.vector(A)
  b = as.vector(b)

  # the current counts of each gene in this sample
  s = as.vector( x %*% A )

  # the current proportions
  b1 = s / sum(s)

  # what the counts should be (if the original counts
  # are reapportioned to match the observed proportions)
  s1 = as.vector( b * sum(s) )

  w = t( t(x) * A )
  w = w / apply(w %*% A, 1, sum)

  x1 = x + w * (s1 - s)

  # print how far off this is
  b.new = as.vector( x1 %*% A )
  b.new = b.new / sum(b.new)
# cat(range(b.new - b), "")    # how far off this is; should be tiny

  x1
}

# Approximates a region constrained such that x >= 0.
# Args:
#   A, b - these give the constraint that
#     forall j. A x[,j] \propto b[,j]
# Returns: list with components
#   x - the posterior counts
#   FIXME update.stats - matrix with mean and variance of update sizes
approx.region.dirichlet = function(A, b,
  converge.tolerance = 1e-3, max.iters=100) {

  # for now, we're assuming each row of b sums to 1
  # (as we're not including any variance in the constraints)
  b = b / apply(b, 1, sum)

  # the messages from each factor; initially flat
  m = array(0, dim=c(ncol(b), ncol(A), nrow(b)),
    dimnames=list(colnames(b), colnames(A), rownames(b)))

  prior = 0 * m[,,1]
  # XXX hack to deal with case when many genes are aggregated with "_"
  prior["_",] = 2* (19723 - ncol(b)) - 1
  prior = prior + 1

  # the posterior
  x = prior + apply(m, c(1,2), sum)

  # messages to and from each factor
  m.to = 0 * m
  m.from = 0 * m

  # for tracking convergence
  update.stats = NULL

  for(iter in 1:max.iters) {

    # compute messages to each factor
    # (note that x and m are different sizes, so this computes
    # x - m[,,i], for each i)
    m.to = as.vector(x) - m

    # compute what each factor would "like" the result to be
    for(i in 1:nrow(A)) {
cat(backspace, i, "")
      m.from[,,i] =
        dirichlet.cond(m.to[,,i] + 1, A[i,], b[i,]) - 1
    }

    # update messages from each factor
    # (again, m.from and x are different sizes; see above about m.to)
    m.change = m.from - as.vector(x)
    m = 0.05 * m.change + m       # includes damping

    # update posterior (as sum of messages from each factor)
    x = prior + apply(m, c(1,2), sum)

    # for tracking convergence
    update.stats = rbind(update.stats, c(mean=mean(m.change), sd=sd(m.change)))
    cat(backspace, iter, update.stats[nrow(update.stats),], "\n")

    # stop if we've converged
    if (max(m.change) < converge.tolerance)
      break
  }
cat("\n")

  x1 = x + 1
  b1 = A %*% t(x1)
  b1 = b1 / apply(b1, 1, sum)
  list(x = x1, b1 = b1, update.stats = update.stats)
}

if (TRUE) {
# test data sets
t1 = list(
  A = rbind(c(1,1,1,1,0),
    c(0,0,0,1,1)),
  b = rbind(c(0.5,0.5),
    c(0.8,0.2)))
t1$A = t1$A / apply(t1$A, 1, sum)
rownames(t1$A) = c("a", "b")

x = matrix(sample(10), nrow=2)

A = t1$A
b = t1$b
max.iters=2

# x1 = dirichlet.cond(x, A, b)

r = approx.region.dirichlet(t1$A, t1$b, max.iters=10)
x1 = r$x
x1 = x1 / apply(x1, 1, sum)
}

