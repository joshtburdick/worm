# practice using gmake

%.sampling :
	Rscript unmix_sampling.r --args $@

genes = $(shell cut -f1 ~/data/image/exprCell.tsv  | tail --lines=+2 | tr -d ['\"'])

s = $(genes:%=%.sampling)

all : $(s)

