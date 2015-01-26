#!/usr/bin/perl -w
# Makes a file containing the merged upstream regions.
# Deprecated.

open IN, "sort -k1,1 -k2,2n upstreamRegionsWS220_1kb.bed | bedtools merge -i |" || die;
open OUT, ">upstreamRegionsWS220_1kb_merged.bed" || die;

while (<IN>) {
  chomp;
  my($chr, $a, $b) = split /\t/;
  print OUT (join "\t", ($chr, $a, $b, "upstreamRegion", "0", "+"));
  print OUT "\n";
}

close IN;
close OUT;

