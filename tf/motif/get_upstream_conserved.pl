#!/usr/bin/perl -w
# Writes out upstream regions with some degree of conservation.
#
# Note that Bio::DB::BigWig can be somewhat painful to install.
# Reads from standard input, writes to standard output.

use strict;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::DB::Fasta;

# command line args

# the file of regions to search
my $region_file = $ARGV[0];

# conservation cutoffs to use (between 0 and 1)
my $cutoff = $ARGV[1];       #  0.50;

# configuration

# set this to wherever the .fa files are
my $fasta = Bio::DB::Fasta->new(
  "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa");

# the conservation track to use
my $wig = Bio::DB::BigFile->bigWigFileOpen(
  "/murrlab/seq/igv/conservation/ce10.phastCons7way.bw");

### end configuration

open IN, "<$region_file" || die;
while (<IN>) {
  chomp;
  my($chr, $a, $b, $name, $score, $strand) = split /\t/;

  $chr =~ s/chr//;

  write_conservation_masked_fasta("$name upstream region",
    $chr, $a, $b, $strand);
}
close IN;

# Prints FASTA-formatted DNA at a region, masked at some
# amount of conservation.
# Args:
#   name - descriptive name
#   chr, a, b, strand - the location
# Side effects: prints FASTA-formatted sequence to standard output
sub write_conservation_masked_fasta {
  my($name, $chr, $a, $b, $strand) = @_;
  die if not ($strand eq "+" || $strand eq "-");
  my $n = $b - $a;

  # get sequence
  my $dna = ($strand eq "+" ?
    $fasta->seq($chr, $a=>$b-1) : $fasta->seq($chr, $b-1=>$a));

  # get bins of statistical summary data
  # (this should be per-base, so we just use max)
  my $r = $wig->bigWigSummaryArray("chr$chr",$a=>$b,bbiSumMax,$n);

  # if there were no results, only print sequence if
  # conservation cutoff was 0 (and then return)
  if (not defined $r) {
#    print "N" x $n;
#    print "\n";
    if ($cutoff == 0) {
      print ">$name, $chr:$a-$b ($strand), phastCons7way >= $cutoff\n";
      print "$dna\n";
    }
    return;
  }

  my @conservation = @$r;  
  if ($strand eq "-") {
    @conservation = reverse(@conservation);
  }

  # get the sequence
  my $s = "";
  my $num_conserved = 0;
  for (my $i=0;$i<$n;$i++) {
    my $c = $conservation[$i];

    # print this base if either the cutoff is 0, or the cutoff
    # is nonzero, and this is at least that much conserved
    if (($cutoff == 0) || (defined $c && $c >= $cutoff)) {
      $s = $s . substr($dna,$i,1);
      $num_conserved++;
    }
    else {
      $s = $s . "N";
    }
  }

  # if sequence is long enough, print it (with a header)
  if ($num_conserved >= 10) {
    print ">$name, $chr:$a-$b ($strand), phastCons7way >= $cutoff\n";
    print "$s\n";
  }
}

