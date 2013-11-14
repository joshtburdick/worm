#!/usr/bin/perl -w
# Writes out a histogram of the amount of conservation in each
# of a set of regions.
#
# Note that Bio::DB::BigWig can be somewhat painful to install.
#
# Input (from stdin): regions to get amount of conservation of
#   (in BED4 format: chromosome, start, end, name of region)
# Output (to stdout): histogram of amount of conservation in
#   each region.

use strict;

use POSIX;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::DB::Fasta;

# the conservation track to use
my $wig = Bio::DB::BigFile->bigWigFileOpen(
  "/murrlab/seq/igv/conservation/ce10.phastCons7way.bw");

# number of bins in the histogram
my $num_bins = 20;

# maximum value of histogram (we assume min is 0)
my $max = 1;

# print header
print "name";
foreach my $i (0..($num_bins-1)) {
  print "\t" . ($i / $num_bins);
}
print "\n";

my $scale = $num_bins / $max;

while (<>) {
  chomp;
  my($chr, $a, $b, $name) = split /\t/;

  $chr =~ s/chr//;

  write_conservation_hist($chr, $a, $b, $name);
}

# Writes out a histogram of the number of bases with a
# given amount of conservation in each region.
sub write_conservation_hist {
  my($chr, $a, $b, $name) = @_;

  my @hist;
  foreach my $i (0..($num_bins-1)) {
    $hist[$i] = 0;
  }

  my $n = $b - $a;

  # get bins of statistical summary data
  # (this should be per-base, so we just use max)
  # we compute the regions on the "+" strand
  my $r = $wig->bigWigSummaryArray("chr$chr",$a=>$b,bbiSumMax,$n);

  # if there were no results, don't print anything
  # (??? and print a warning?)
  if (not defined $r) {
    print "$name\t$n";
    foreach my $i (1..($num_bins-1)) {
      print "\t0";
    }
    print "\n";
    return;
  }

  my @conservation = @$r;

  # loop through each base, adding to the counts
  for(my $i=0; $i<$n; $i++) {

    my $c = $conservation[$i];

    # note that "no conservation", a common case, is fast
    if (!defined($c)) {
      $hist[0]++;
    }
    else {
      my $j = floor($scale * $c);
      if ($j < 0) {
        $j = 0;
      }
      if ($j >= $num_bins) {
        $j = $num_bins - 1;
      }
      $hist[$j]++;
    }
  }

  # write counts
  print join "\t", ($name, @hist);
  print "\n";
}

