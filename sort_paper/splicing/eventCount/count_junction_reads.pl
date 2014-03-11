#!/usr/bin/perl -w
# Counts reads across splice junctions, using sjcount.

use strict;

my $read_dir = "/murrlab/seq/tophat2/WS220_20140111/";

my $output_dir = "spliceJunctionCounts/";

# Counts all the junctions in the files in a directory.
# Args:
#   read_dir - a directory full of .bam files
#   output_dir - directory in which to write files
# Side effects: writes counts to output_dir.
sub count_junctions {
  my($read_dir, $output_dir) = @_;

  print "[counting reads in $read_dir, writing to $output_dir]\n";

  system("mkdir -p $output_dir");

  foreach my $f (<$read_dir/*.bam>) {
    print "[counting $f]\n";

    die if not ($f =~ /\/([^\/]+)\.bam/);
    my $output_name = $1;

    #   Note that we reverse-complement both reads, for
    # convenience, as that matches the gene direction.
    #   For now, omitting the "splice spanning counts" that
    # would be printed with "-ssc".
    system("sjcount -bam $f -quiet -read1 1 -read2 1 " .
      " -ssj $output_dir/$output_name.tsv");
  }
}

foreach my $a (`ls $read_dir`) {
  chomp $a;
  print "$a\n";

  count_junctions("$read_dir/$a", "$output_dir/$a");
}

