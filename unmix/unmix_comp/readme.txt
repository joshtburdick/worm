
This is data and software for the paper

Deconvolution of gene expression in single cells from population data

If you have questions, please e-mail

Josh Burdick (jburdick@mail.med.upenn.edu) or
John Murray (jmurr@mail.med.upenn.edu).

The software is written in R, and requires the following R packages:

corpcor
limSolve (v.1.5.3)
Matrix
Rcpp

It has mostly been tested using R 2.15.1 running
on x86 64-bit Linux (with those package versions),
but should run with other versions of these.
(Rcpp is probably the most complex of these dependencies,
but is only used for the sampling.)

To run the tests, start R in the "src/" directory, and at
the prompt, type

source("run_all.r")

This will create simulated data, pick reporters, and simulate
several methods of unmixing.

Note that for the timing results of the truncated pseudoinverse
(with and without correlation), we used lsei() with the
"type=1" option, as that was faster, and worked for the datasets
being benchmarked there. However, it didn't work for all datasets,
so we used the "type=2" option.

Josh Burdick
20130118

