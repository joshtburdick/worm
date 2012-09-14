#!/usr/bin/perl -w
# Gets upstream sequence for a set of genes.
# Usage:
#   get_upstream.pl gene_list rel_upstream_start rel_upstream_end
# where
#   gene_list - a list of gene identifiers
#   rel_upstream_start, rel_upstream_end - region relative to
#     gene starts to include. For "3kb upstream", use
#     -3000 and 0, respectively.
# Output: FASTA file of upstream sequences

use strict;

my $gene_list = $ARGV[0];
my $rel_upstream_start = $ARGV[1];
my $rel_upstream_end = $ARGV[2];

# Location of BED file of genes.
my $gene_bed_file = "refGene_WS220.bed";

# temporary BED file
my $bed_tmp = "tmp.bed";

# read in list of genes
my %genes = ();
open GENE, "<$gene_list" || die;
while (<GENE>) {
  chomp;
  $genes{$_} = 1;
}
close GENE;

open OUT, ">$bed_tmp" || die;

open IN, "<$gene_bed_file" || die;
while (<IN>) {
  chomp;
  my($chr, $a, $b, $name, $score, $strand) = split /\t/;

  # only include genes in the list
  next if not ($genes{$name});

  # start and end coordinates, relative to start
  my $a1 = $strand eq '+' ?
    $a + $rel_upstream_start :
    $b + $rel_upstream_end;
  my $b1 = $strand eq '+' ?
    $a + $rel_upstream_end :
    $b + $rel_upstream_start;

  # skipping cases which would involve truncation
  # FIXME should also check for "off the end of the chromosome"
  next if ($a1 < 1 || $b1 < 1);

  print OUT (join "\t", ($chr, $a1, $b1, $name, 0, $strand));
  print OUT "\n";
}

close IN;
close OUT;

system "./get_dna_in_region.pl $bed_tmp";
unlink($bed_tmp);

