# practice using gmake

%.sampling.Rdata :
	Rscript unmix_sampling_1.r --args $@

genes = $(shell echo a b)


s = $(genes:%=%.sampling.Rdata)



all : $(s)

