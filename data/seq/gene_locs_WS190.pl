#!/usr/bin/perl -w
# Converts a list of gene annotations downloaded from Wormbase
# into a file with genes with the same name merged.
# Requires Bedtools.

# XXX this is printing the name of an arbitrary transcript,
# rather than the locus. This may be fixable.

use strict;

use Graph::UnionFind;

my $infile = "~/data/wormbase/c_elegans.WS190.annotations.gff2.gz";

# Creates a BED6 file containing merged exons, and the names of all
# genes overlapping each exon.
sub get_exons {
  open OUT, "| bedtools sort -i - > exons_WS190.bed" || die;

  # first, the majority of genes (from Wormbase)
  open IN, "gunzip -c $infile | egrep \"(Coding_transcript|Non_coding_transcript)\" | bedtools sort -i - |" || die;
  while (<IN>) {
    my @a = split /\t/;
    $a[0] =~ s/^CHROMOSOME_//;

    if (($a[1] eq "Coding_transcript" ||
        $a[1] eq "Non_coding_transcript")
        && $a[2] eq "exon") {
      die @a if not $a[8] =~ /Transcript "(\S+)"/;
      my $name = $1;
      print OUT (join "\t", $a[0], $a[3], $a[4], $name, 0, $a[6]);
      print OUT "\n";
    }
  }
  close IN;

  # then, the non-coding RNAs
  my $dir = "../../../data/expression/lincRNA/";
  open IN,
    "cat $dir/ancRNAs_W3PSeq3_ce6.gtf $dir/lincRNAs_W3PSeq3_ce6.gtf |"
    || die;
  while (<IN>) {
    my @a = split /\t/;
    $a[0] =~ s/^chr//;

    if ($a[2] eq "exon") {
      die @a if not $a[8] =~ /gene_id "(\S+)"/;
      my $name = $1;
      print OUT (join "\t", $a[0], $a[3], $a[4], $name, 0, $a[6]);
      print OUT "\n";
    }
  }
  close IN;

  close OUT;

  system("bedtools merge -s -nms -i exons_WS190.bed " .
    "| bedtools sort -i - > merged_exons_WS190.bed");
}

# Given the set of exons, groups all the genes overlapping together.
# Args:
#   exon_file - name of the file of exons to read (in BED6 format)
#   gff_output_file - name of the GFF file to output
# Side effects: writes a GFF file
sub group_exons {
  my($exon_file, $gff_output_file) = @_;

  my $uf = Graph::UnionFind->new;

  # first, pick a "canonical" gene for each exon
  open IN, "<$exon_file" || die;
  while (<IN>) {
    chomp;
    my @a = split /\t/;

    # this should be the names of all exons overlapping here
    my @genes = split /;/, $a[3];
    foreach my $g (@genes) {
      $uf->union($genes[0], $g);
    }
  }
  close IN;

  # write out file as GFF format, with canonical name
  open IN, "<$exon_file" || die;
  open OUT, ">$gff_output_file" || die;
  while (<IN>) {
    chomp;
    my @a = split /\t/;
    my $name = $a[3];
    my @genes = split /;/, $a[3];
    my $canonical_name = $uf->find($genes[0]);
    print OUT (join "\t", ($a[0], "WormBase", "exon",
      $a[1], $a[2], 0, $a[4], ".", $canonical_name));
    print OUT "\n";
  }
  close IN;
  close OUT;
}

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

# Converts from GFF format to BED format.
# Args:
#   gff_in - GFF file to convert
#   bed_out - BED file to write to
sub gff_to_bed {
  my($gff_in, $bed_out) = @_;

  my $wormbase_name = get_wormbase_names();

  # first, read in all records, indexed by gene
  open IN, "bedtools sort -i $gff_in |" || die;
  my %g = ();
  while (<IN>) {
    chomp;
    my @a = split /\t/;
    push @{$g{$a[8]}}, \@a;
  }
  close IN;

  # then, write out each gene
  open OUT, "| bedtools sort -i - > $bed_out" || die;
  foreach my $id (keys %g) {

    # sort by start location (may not be necessary)
    my @r = @{$g{$id}};
    @r = sort { $a->[3] <=> $b->[3] } @r;

    # aggregate all this information
    my $a = $r[0]->[3];
    my $b = $r[@r-1]->[4];
    my @sizes = ();
    my @starts = ();
    foreach (@r) {
      push @sizes, $_->[4] + 1 - $_->[3];
      push @starts, $_->[3] + 0 - $a;
    }

    # if possible, convert ID to be a Wormbase gene name
    if (defined $wormbase_name->{$id}) {
      $id = $wormbase_name->{$id};
    }

    print OUT join "\t", ($r[0]->[0], $a-1, $b, $id, 0,
      $r[0]->[6], $a-1, $b, "0,0,0", (@sizes+0),
      (join ",", @sizes), (join ",", @starts));
    print OUT "\n";
  }
  close OUT;
}

get_exons();
group_exons("merged_exons_WS190.bed", "merged_genes_WS190.gff");
gff_to_bed("merged_genes_WS190.gff", "merged_genes_WS190.bed");

# XXX the coordinates for these are off-by-one; not worrying
# about this for now, but still cleaning up these files
unlink("exons_WS190.bed");
unlink("merged_exons_WS190.bed");

