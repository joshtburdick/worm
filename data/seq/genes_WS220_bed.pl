#!/usr/bin/perl -w
# Writes out the genes as a .bed file.
# XXX possibly not used, as this loses the exon structure.

use strict;

# for SAM->BAM conversion
my $genome_dict = "/var/tmp/data/tophat2/WS220/genome_with_MtDNA.dict";

# first, convert to SAM format
system("gtf_to_sam genes_tophat_WS220.gtf genes_WS220.sam");

# convert that to bam
system("cat $genome_dict genes_WS220.sam | samtools view -bS - > genes_WS220.bam");

system("bedtools bamtobed -color 0,0,0 -bed12 -i genes_WS220.bam > genes_WS220.bed");



