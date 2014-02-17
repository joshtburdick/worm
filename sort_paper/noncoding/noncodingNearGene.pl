#!/usr/bin/perl -w
# Finds noncoding genes near genes.
# Uses bedtools.

# my $chr_sizes = "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.fai";

# gene bounds, and class (e.g. coding or non-coding)
my $gene_bounds = "../../data/seq/merged_genes_WS220.bed";
my $gene_class = "../../data/seq/gene_class.tsv.gz";

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

# Gets genes of a particular class.
sub get_genes_by_class {
  my($class_pattern, $output_file) = @_;

  # get list of genes of that class
  my %g = ();
  open IN, "gunzip -c $gene_class |" || die;
  while (<IN>) {
    chomp;
    my($gene, $class) = split /\t/;

    # if it's in this class...
    if ($class =~ $class_pattern) {

      # possibly translate the name
      if (defined $wormbase_name->{$gene}) {
        $gene = $wormbase_name->{$gene};
      }
      $g{$gene} = 1;
    }
  }
  close IN;

  open IN, $gene_bounds || die;
  open OUT, ">$output_file" || die;
  while (<IN>) {
    my @a = split /\t/;
    if (defined $g{$a[3]}) {
      print OUT (join "\t", @a);
    }
  }
  close IN;
  close OUT;
}

# bounds of (some) non-coding genes
get_genes_by_class("^(anc|linc)RNA\$", "nc.bed");

# bounds of other genes
get_genes_by_class("^protein_coding\$", "coding.bed");

# find the nearest gene (a la Nam and Bartel)
system("bedtools closest -d -a nc.bed -b coding.bed > closestToNoncoding.tsv");

# find pairs of nearby genes
system("bedtools closest -io -a coding.bed -b coding.bed > nearbyGenes.tsv");

