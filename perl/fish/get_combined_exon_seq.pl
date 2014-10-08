#!/usr/bin/perl -w
# Script to get sequence for FISH probe design.
# Requires BEDtools.
#
# Usage:
#   get_combined_exon_seq.pl genes.txt output.fa
# where
#   genes.txt - a list of gene names, one per line
#   output.fa - where to write output

die if not @ARGV == 2;
my $genes = $ARGV[0];
my $output_file = $ARGV[1];

# Genomic sequence, in FASTA format (this is from the Tophat2
# downloads.)
my $genome_fasta = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/WholeGenomeFasta/genome.fa";

# BED file of exonic regions for each gene. (This is not the BED
# file of mRNAs, but rather the exonic regions merged together.)
my $exon_bed = "/home/jburdick/gcb/git/data/seq/combined_exons_WS220.bed";

# my $tmp_bed = "tmp.bed";     # FIXME use tmpnam()

# put the gene names to get in a hash
open IN, "<$genes" || die;
my %g = ();
while (<IN>) {
  chomp;
  $g{$_} = 1;
}
close IN;

# get the relevant lines from the file of genes
open IN, "<$exon_bed" || die;
open OUT, "|bedtools getfasta -s -name -split -fi $genome_fasta -bed - -fo $output_file";
while (<IN>) {
  my $line = $_;
  chomp $line;
  my @a = split /\t/, $line;
  if ($g{$a[3]}) {
    print OUT "$line\n";
  }
}
close IN;
close OUT;

