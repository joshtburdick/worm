#!/usr/bin/perl -w
# Utility to reheader a directory of .bam files, for
# submitting them to the SRA.

use strict;

my $new_header = $eARGV[0];
my $input_dir = $ARGV[1];
my $output_dir = $ARGV[2];

# loop through the files
foreach (<$input_dir/*.bam>) {
  chomp;
  /\/([^\/]+)\.bam$/;
  my $a = $1;

  print "$a\n";

  system("samtools reheader $new_header $input_dir/$a.bam > $output_dir/$a.bam");

  system("samtools index $output_dir/$a.bam");
}


