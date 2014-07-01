// Utilities for vector math.
// Here, what's called a "vector" is just an array of numbers.

// XXX this probably isn't how people usually do
// namespaces in JavaScript.
vectorMath = {};


// The sum of a vector.
vectorMath.sum = function(x) {
  s = 0;
  for(i=0; i<x.length; i++)
    s = s + x[i];
  return s;
};

// Some stats about a vector.
vectorMath.meanAndStandardDeviation = function (x) {
  n = x.length;
  s1 = vectorMath.sum( x );
  s2 = vectorMath.sum( x.map(
    function(a) { return a * a; } ));
  m1 = s1 / n;

  return { mean: m1,
    stdDev:
      Math.sqrt( (n/(n-1)) * (s2/n - m1*m1) ) };
};

// Standardizes a vector (to have mean 0 and s.d. 1.)
vectorMath.standardize = function (x) {
  s = vectorMath.meanAndStandardDeviation(x);
  return x.map(
    function(a) { return (a - s.mean) / s.stdDev; });
}

// Dot product of two vectors.
vectorMath.dot = function(x, y) {
  s = 0;
  for(i=0; i<x.length; i++)
    s = s + x[i] * y[i];
  return s;
}


