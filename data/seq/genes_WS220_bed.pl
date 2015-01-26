#!/usr/bin/perl -w
# Writes out the genes as a .bed file.
# Most of the work is done by gtf2bed.pl; this also renames
# them (when possible), and colors them blue.
#
# Possibly not used.

use strict;

# Reads in mapping from transcript names to human-readable
# WormBase names (when available.)
sub get_wormbase_names {
  open IN, "gunzip -c ../../../data/wormbase/c_elegans.PRJNA13758.WS240.xrefs.txt.gz |" || die;
  my %h = ();
  while (<IN>) {
    chomp;
    my @a = split /\t/;
    if (!($a[2] eq ".")) {
      $h{$a[0]} = $a[2];
      $h{$a[3]} = $a[2];
    }
  }

  return \%h;
}

my $wormbase_name = get_wormbase_names();

open IN, "gtf2bed.pl genes_tophat_WS220.gtf |" || die;
open OUT, ">genes_WS220.bed" || die;
while (<IN>) {
  my @a = split /\t/;
  my $name = $a[3];
  if (defined $wormbase_name->{$name}) {
    $name = $wormbase_name->{$name};
  }
  next if $name =~ /^21ur-/;
  $a[8] = "0,0,178";
  print OUT (join "\t", @a);
}
close IN;
close OUT;

