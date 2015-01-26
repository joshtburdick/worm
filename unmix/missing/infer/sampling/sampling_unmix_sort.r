# Attempt at simultaneously estimating expression,
# and the sort matrix.
# (Probably omitting double-sorted fractions for the moment.)

source("git/utils.r")

source("git/unmix/missing/infer/sampling/mixture_prior.r")

# the sorting matrix to use
source("git/sort_paper/unmix/sortMatrix.r")
# pull these numbers somewhat away from 0 and 1
m = 0.8 * m.unnormalized + 0.1

# for now, omitting double-sorted fractions, and arbitrarily
# picking one of the fractions with replicates
m = m[ 1:14 , ]
rownames(m) = sub("pha-4", "pha-4 5/9", rownames(m))
rownames(m) = sub("cnd-1", "cnd-1 8/19", rownames(m))

# the expression data, in reads per million (transposed)
x = t(as.matrix(read.tsv("git/cluster/readsPerMillion.tsv")))

# Construct a simple prior. For now, we assume a flat
# prior on the volume.
sort.prior = sort.prior.dirichlet.1(rep(1, 1341), m, 10)

x0 = sort.init(sort.prior)   # XXX rename this
m0 = sort.matrix.normalize(x0)

sort.fractions = sort(intersect(rownames(m0), rownames(x)))
x = x[ sort.fractions , ]

# for prototyping (and to reduce noise), only keeping
# most-highly-expressed things
x = x[ , apply(x, 2, max) >= 100 ]
# normalize the rows
x = x / apply(x, 1, sum)

m0 = m0[sort.fractions,]

# foo = expr.init(m0, x)
# baz = m0 %*% foo$x - x
# print(sum(baz^2))

foo = find.init(x0, x)

