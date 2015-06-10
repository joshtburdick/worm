#!/usr/bin/perl -w
# Runs FIMO on other similar genomes.

use strict;

# my $output_base = "/home/jburdick/tmp/fimo/";
my $output_base = "/media/jburdick/disk2/jburdick/fimo/";

# relevant scripts
my $run_fimo = "../../../../perl/motif/run_fimo_bam_output.pl";

# Runs FIMO on one set of motifs.
sub run_fimo {
  my($name, $fa) = @_;

  system("$run_fimo ../../../../data/tf/hughes/Ce_1.02.meme $output_base/$name" .
    "_" . "Ce_1.02 $fa");
  foreach my $org ("Hs_1.02", "Mm_1.02", "Dm_1.02") {
    system("$run_fimo $org.meme $output_base/$name" . "_" . "$org $fa");
  }
}

# run this on various genome(s)
run_fimo("cb3",
  "/media/jburdick/disk2/data/ftp/hgdownload.soe.ucsc.edu/goldenPath/cb3/bigZips/cb3.fa");

