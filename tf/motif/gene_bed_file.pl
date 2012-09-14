#!/usr/bin/perl -w
# Gets list of genes from a _-format file to .BED format.

open IN, "</home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Annotation/Genes/refGene.txt" || die;

my %h = ();

while (<IN>) {
  chomp;
  my @a = split /\t/;
  $h{$a[12]} = join "\t", ($a[2], $a[4], $a[5], $a[12], 0, $a[3]);
}

close IN;

foreach (values %h) {
  print "$_\n";
}

