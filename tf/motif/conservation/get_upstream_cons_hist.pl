#!/usr/bin/perl -w
# Computes upstream conservation for several sizes of
# upstream region.

use strict;

foreach my $s (1..5) {
  print("processing $s kb upstream\n");
  my $regions_file =
    "../motifCount/regions/upstream_" . $s . "kb_WS220.bed";
  my $output_file = "cons_hist_WS220_" . $s . "kb_upstream.tsv.gz";
  system("../../../perl/motif/get_conservation_hist_in_region.pl " .
    "< $regions_file | gzip -c > $output_file");
}

