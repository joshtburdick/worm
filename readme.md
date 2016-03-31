
This repository contains various code I wrote in John Murray's
lab at UPenn.

Unfortunately, there is no unified architecture to re-run all the analysis.
Most of the code is written in R (see www.r-project.org), and can be run
by going to an empty directory, and:

- at the UNIX prompt, typing

mkdir worm
cd worm
git clone -b dev https://github.com/joshtburdick/worm.git git
R

- at the R prompt, typing

source("git/plot/seqLogoSVG.r")
test1()

(Admittedly, "git" shouldn't be in all the pathnames.)

For instance, much of the code used in the paper describing gene
expression in FACS-sorted C. elegans cells is in the directory
git/sort_paper .

Much of the code used in the earlier unmixing simulation paper
is in the directory git/unmix/unmix_comp (on the dev branch), and
has to be run from that directory.

Josh Burdick, 20160331

