# The pseudoinverse (on fractions, rather than enrichments.)

source("git/sort_paper/unmix/otherMethods/crossValAccuracyRpm.r")

# Unmixes using the pseudoinverse.
unmix.pseudoinverse1 = function(m, r) {
  mn = m / apply(m, 1, sum)
  x = t( pseudoinverse(mn) %*% t(r))
cat("mean x < 0 =", mean(x < 0), "\n")
  x[ x < 0 ] = 0
  list(x = x)
}

write.crossval.graphs(unmix.pseudoinverse1, "pseudoinverse1")

