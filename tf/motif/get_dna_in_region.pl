#!/usr/bin/perl -w
# Gets some chunks of sequence.
# Input: a BED file of regions.

use strict;

use Bio::DB::Fasta;

# set this to wherever the .fa files are
my $fasta = Bio::DB::Fasta->new(
  "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa");

while (<>) {
  chomp;
  my($chr, $a, $b, $name, $score, $strand) = split /\t/;

  $chr =~ s/chr//;

  print ">$name\n";
  if ($strand eq "+") {
    print $fasta->seq($chr, $a => $b) . "\n";
  }
  elsif ($strand eq "-") {
    print $fasta->seq($chr, $b => $a) . "\n";
  }
  else {
    die "unknown strand $strand\n";
  }
}

