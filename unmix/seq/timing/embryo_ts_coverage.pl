#!/usr/bin/perl -w
# Gets coverage of a .wig file at some locations.

use strict;

my $gene_bed_file = "../quant/geneBounds.tsv";
my $bw_dir = "/murrlab/seq/igv/expression/embryo_ts/";
my $output_dir = "embryo_ts_coverage";

system("mkdir -p $output_dir");

# rename to remove the "chr"
system("sed -e 's/^chr//;' < $gene_bed_file > geneBoundsNoChr.bed");
system("bedSort geneBoundsNoChr.bed geneBoundsNoChr.bed");

foreach my $f (<$bw_dir/*.bw>) {
  die if not ($f =~ /\/([^\/]+)\.bw/);
  my $output_file = "$output_dir/$1.tsv";

  print "[quantifying $f]\n";

  system("bigWigAverageOverBed $f geneBoundsNoChr.bed $output_file");

  system("bzip2 --best $output_file");
}


