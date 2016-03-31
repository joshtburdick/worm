# Tests of Message.

library("Matrix")

source("git/unmix/ept/practice/Message.R")

a = new("Message", x = list(
  a = c(1:10),
  r = matrix(1:3, nrow=2, ncol=3),
  y = Matrix(1:6, nrow=2)))

# XXX it might be nice to not have to put
# Message around the lists here.
b = new("Message", x = list(a = a, z = 1))


z = new("Message", x=list(a = 3 * c(1:10)))


print(1+a)
print(a/10)

print(z + z)
print(a + z)

# test of whether "recycling" works with different-dimension
# arrays (for better or worse)
y = new("Message", x=list(r = c(0:2)))
print(b + y)
# XXX doesn't seem to be working

