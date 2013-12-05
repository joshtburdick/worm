#!/usr/bin/perl -w
# Constructs a .gff file for Tophat, which includes the
# shipped annotation, plus lincRNAs and anrRNAs.
# (This is intended for mapping, not so much for quantification.)

use strict;

my $ncrna_dir = "~/gcb/data/expression/lincRNA";

# first, include the given annotation file
system("cp ~/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Annotation/Genes/genes.gtf genes_tophat_WS220.gtf");

# tack on the non-coding RNAs
system("sed -e 's/^chr//' $ncrna_dir/ancRNAs_W3PSeq3_ce10.gtf " .
  " >> genes_tophat_WS220.gtf");
system("sed -e 's/^chr//' $ncrna_dir/lincRNAs_W3PSeq3_ce10.gtf " .
  " >> genes_tophat_WS220.gtf");

# sort it
system("bedtools sort -i genes_tophat_WS220.gtf > tmp.gtf");
system("mv tmp.gtf genes_tophat_WS220.gtf");

