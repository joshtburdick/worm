#!/usr/bin/perl -w
# Writes out upstream regions with some degree of conservation
# (the actual regions, not the DNA.)
#
# Note that Bio::DB::BigWig can be somewhat painful to install.
#
# Args:
#   cutoff - the conservation cutoff to use
#   regions - BED file of regions to compute conservation of
# Output (to stdout): just the portions with at least that
#   amount of conservation.
# Note that if input regions overlap, then output regions will also.

# use strict;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::DB::Fasta;

my $cutoff = $ARGV[0];

my $min_size = 10;

# the conservation track to use
my $wig = Bio::DB::BigFile->bigWigFileOpen(
  "/murrlab/seq/igv/conservation/ce10.phastCons7way.bw");

open IN, "<" . $ARGV[1] || die;
while (<IN>) {
  chomp;
  my($chr, $a, $b, $name, $score, $strand) = split /\t/;

  $chr =~ s/chr//;

  write_conserved_regions($cutoff, $chr, $a, $b, $name, $score, $strand);
}
close IN;

# Prints regions which have at least some amount of conservation.
# Args:
#   the fields from a .bed-formatted file
# Side effects: prints the regions
sub write_conserved_regions {
  my($cutoff, $chr, $a, $b, $name, $score, $strand) = @_;
  die if not ($strand eq "+" || $strand eq "-");
  my $n = $b - $a;

  # get bins of statistical summary data
  # (this should be per-base, so we just use max)
  # we compute the regions on the "+" strand
  my $r = $wig->bigWigSummaryArray("chr$chr",$a=>$b,bbiSumMax,$n);

  # if there were no results, don't print anything
  # (??? and print a warning?)
  if (not defined $r) {
    return;
  }

  my @conservation = @$r;

  # the regions (as a list of regions to print)
  my @block_sizes = ();
  my @block_starts = ();

  my $region_start = undef;

  my $count = 0;
  for(my $i=0; $i<$n; $i++) {

    # if it's undef, conservation is presumed to be zero
    my $c = (defined $conservation[$i] ? $conservation[$i] : 0);

    # if conservation is now above the cutoff,
    # and a region hasn't started yet, start one
    if ($c >= $cutoff && !defined($region_start)) {
      $region_start = $i;
    }

    # if conservation has dropped below the cutoff,
    # write out the region's start and end
    if ($c < $cutoff && defined($region_start)) {
      my $size = $i - $region_start;
      if ($size >= $min_size) {
        $count++;
        push @block_sizes, $size;
        push @block_starts, $region_start;
      }
      $region_start = undef;
    }
  }

  # possibly tack on the final region
  if (defined $region_start) {
    my $size = $n - $region_start;
    if ($size >= $min_size) {
      $count++;
      push @block_sizes, $size;
      push @block_starts, $region_start;
    }
  }

  # print the regions as "exons" (if there were any)
  if (@block_sizes > 0) {
    my $n = @block_sizes + 0;
    my $thick_start = $a + $block_starts[0];
    my $thick_end = $a + $block_starts[$n-1] + $block_sizes[$n-1];

    print join "\t", ($chr, $a, $b, $name, $score, $strand,
      $thick_start, $thick_end, "0,0,0", $n,
      join (",", @block_sizes), join (",", @block_starts));
    print "\n";
  }
}

