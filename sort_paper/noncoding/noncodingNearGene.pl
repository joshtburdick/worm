#!/usr/bin/perl -w
# Finds noncoding genes near genes.
# Uses bedtools.

# my $chr_sizes = "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.fai";

# the gene bounds, sorted (really this should be sorted
# in the first place)
system("bedtools sort -i ../../unmix/seq/quant/geneBounds.tsv > geneBounds.tsv");
system("egrep '(anr|linc)-' geneBounds.tsv | bedtools sort -i - > nc.bed");

# system("bedtools closest -d -a geneBounds_WS220.tsv -b tmp.bed | bedtools sort -i - > noncodingNearCoding.bed");

# first, create the gene bounds, widened
# system("bedtools slop -i geneBounds.tsv -g $chr_sizes -b 10000 > geneRegions.bed");

# then, find the intersection
# system("bedtools intersect -wa -wb -a nc.bed -b geneRegions.bed > noncodingNearGene.tsv");

# also, find the nearest gene (a la Nam and Bartel)
# first, find list of "coding" (or at least neither anc- nor linc-)
system("egrep -v '(anr|linc)-' geneBounds.tsv | bedtools sort -i - > coding.bed");
system("bedtools closest -d -a nc.bed -b coding.bed > closestToNoncoding.tsv");


# clean up
# unlink("coding.bed");
# unlink("nc.bed");
# unlink("geneRegions.bed");

