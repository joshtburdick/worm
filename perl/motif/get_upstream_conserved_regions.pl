#!/usr/bin/perl -w
# Writes out upstream regions with some degree of conservation
# (the actual regions, not the DNA.)
#
# Note that Bio::DB::BigWig can be somewhat painful to install.
#
# Input (from stdin): regions to get conserved portions of
# Output (to stdout): just the portions with at least that
#   amount of conservation.
# Note that if input regions overlap, then output regions may as well.

# use strict;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::DB::Fasta;

my $cutoff = 0.5;

my $min_size = 20;

# the conservation track to use
my $wig = Bio::DB::BigFile->bigWigFileOpen(
  "/murrlab/seq/igv/conservation/ce10.phastCons7way.bw");

while (<>) {
  chomp;
  my($chr, $a, $b, $name, $score, $strand) = split /\t/;

  $chr =~ s/chr//;

  write_conserved_regions(0.5, $chr, $a, $b, $name, $score, $strand);
}

# Prints FASTA-formatted DNA at a region, masked at some
# amount of conservation.
# Args:
#   the fields from a .bed-formatted file
# Side effects: prints the regions
sub write_conserved_regions {
  my($cutoff, $chr, $a, $b, $name, $score, $strand) = @_;
  die if not ($strand eq "+" || $strand eq "-");
  my $n = $b - $a;

  # get bins of statistical summary data
  # (this should be per-base, so we just use max)
  my $r = $wig->bigWigSummaryArray("chr$chr",$a=>$b,bbiSumMax,$n);

  # if there were no results, don't print anything
  # (??? and print a warning?)
  if (not defined $r) {
    return;
  }

  my @conservation = @$r;

  # the regions (as a list of regions to print)
  my @regions = ();

  my $region_start = undef;

  for(my $i=0; $i<$n; $i++) {

    # if it's undef, conservation is presumed to be zero
    my $c = (defined $conservation[$i] ? $conservation[$i] : 0);

    # if conservation is now above the cutoff,
    # and a region hasn't started yet, start one
    if ($c >= $cutoff && !defined($region_start)) {
      $region_start = $i;
    }

    # if conservation has dropped below the cutoff,
    # write out the region
    if ($c < $cutoff && defined($region_start)) {
      if ($i - $region_start >= $min_size) {
        push @regions, (join "\t",
          ($chr, $a + $region_start, $a + $i, $name, $score, $strand)) . "\n";
      }
      $region_start = undef;
    }
  }

  # possibly tack on the final region
  if (defined $region_start) {
    if ($n - $region_start >= $min_size) {
      push @regions, (join "\t",
        ($chr, $a + $region_start, $a + $n, $name, $score, $strand)) . "\n";
    }
  }

  # print the regions, if there were any
  if (@regions > 0) {
    if ($strand eq "+") {
      print join "", @regions;
    }
    else {
      print join "", (reverse @regions);
    }
  }
}




