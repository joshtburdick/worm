#!/usr/bin/perl -w
# Converts a list of gene annotations downloaded from Wormbase
# into a file with genes with the same name merged.
# Requires Bedtools.

use strict;

# use Graph::UnionFind;

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

# Creates a BED6 file containing merged exons, and the names of all
# genes overlapping each exon.
sub get_exons {
  open OUT, "| bedtools sort -i - > exons_WS220.bed" || die;

  my $infile = "genes_tophat_WS220.gtf";

  # first, the majority of genes (from WS220; technically from the
  # iGenomes build, but appears to be the same)
  open IN, "bedtools sort -i $infile |" || die;
  while (<IN>) {
    my @a = split /\t/;
    $a[0] =~ s/^CHROMOSOME_//;

    if (!($a[1] eq "transcript")) {
      die @a if not $a[8] =~ /gene_id "(\S+)"/;
      my $name = $1;

      # skipping genes whose name starts with "21ur-"
      if (defined $wormbase_name->{$name}) {
        $name = $wormbase_name->{$name};
      }

      if (!($name =~ /^21ur-/)) {
        print OUT (join "\t", $a[0], $a[3], $a[4], $name, 0, $a[6]);
        print OUT "\n";
      }
    }
  }
  close IN;

  close OUT;

  system("bedtools merge -s -nms -i exons_WS220.bed " .
    "| bedtools sort -i - > merged_exons_WS220.bed");
}

# Given a BED6 file, writes a BED6 file with one record per name
# (covering the range of all records having that name.)
# Args:
#   in_file - the file to read
#   out_file - the file to write
sub get_gene_bounds {
  my($in_file, $out_file) = @_;

  # record of bounds, each in BED6 format
  my %b = ();

  # group records by name
  open IN, "<$in_file" || die;
  while (<IN>) {
    chomp;
    my @a = split /\t/;

    # have we seen this name before?
    if (!defined($b{$a[3]})) {
      # add initial record
      $b{ $a[3] } = \@a;
    }
    else {
      # update bounds
      my @r = @{ $b{ $a[3] } };
      $r[1] = ($a[1] < $r[1] ? $a[1] : $r[1]);
      $r[2] = ($a[2] > $r[2] ? $a[2] : $r[2]);
      die if not ($a[5] eq $r[5]);
      $b{ $a[3] } = \@r;
    }
  }
  close IN;

  # write the output (and sort it)
  open OUT, "| bedtools sort -i - > $out_file" || die;
  foreach (values(%b)) {
    print OUT (join "\t", @{$_});
    print OUT "\n";
  }
  close OUT;
}

get_exons();
get_gene_bounds("exons_WS220.bed", "merged_genes_WS220.bed");

