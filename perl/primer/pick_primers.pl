#!/usr/bin/perl -w
# Picks primers for a set of genes.
# Standard input: a set of gene identifiers.
# Standard output: a table with columns
#   gene
#   name (gene + _F or _R)
#   sequence

use strict;

# use Bio::DB::BigFile::Constants;
use Bio::DB::Fasta;

# use Bio::Tools::Run::Bowtie;
use Bio::Seq;
use Bio::Tools::GFF;
use Bio::Tools::Run::Primer3;


# various sources of data

# gene locations
# FIXME change this to "regions shared by all transcripts"?
# my $gene_gff = "/home/jburdick/gcb/git/data/seq/genes_tophat_WS220.gtf";
my $gene_bed_file = "/home/jburdick/gcb/git/data/seq/genes_WS220.bed";

# set this to wherever the .fa files of the genome are
my $genome_fasta = Bio::DB::Fasta->new(
  "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa");

# Reads in annotation in BED format.
sub read_annotation_bed {
  my($gene_bed) = @_;
  
  my %feature_by_gene = ();
  open IN, "<$gene_bed" || die;
  while (<IN>) {
    chomp;
    my @a = split /\t/;
    $feature_by_gene{$a[3]} = \@a;
  }
  close IN;

  return \%feature_by_gene;
}

my %feature_by_gene = %{ read_annotation_bed($gene_bed_file) };

# Gets sequence for a set of regions.
# Args:
#   bed_ref - the regions, as a ref to a list of BED-format fields
#   spacer_size - number of N's to put between regions
# Returns: the sequence, as a string, with N's in between regions
sub get_region_seq {
  my($bed_ref, $spacer_size) = @_;
  my @bed = @$bed_ref;
  my $chr = $bed[0];
  $chr =~ s/chr//;
  my $a = $bed[1];
  my $b = $bed[2];
  my $strand = $bed[5];
  my @block_sizes = split /,/, $bed[10];
  my @block_starts = split /,/, $bed[11];
  my $n = @block_sizes + 0;

  my $s = undef;

  # loop through the blocks
  foreach my $i1 (0..($n-1)) {
    # for minus strand, loop through in reverse order
    my $i = $strand eq "+" ? $i1 : (($n-1) - $i1);
    my $a1 = $a + $block_starts[$i];
    my $b1 = $a1 + $block_sizes[$i];
    my $s1 = $strand eq "+" ?
      $genome_fasta->seq($chr, $a1 => $b1) :
      $genome_fasta->seq($chr, $b1 => $a1);
    if (defined $s) {
      $s = $s . (join "", "N" x $spacer_size) . "\n" . $s1;
    } else {
      $s = $s1;
    }
  }

  return $s;
}

# Formats selected parts of a primer result as a table.
sub format_results {
  my($results, $name) = @_;

  # field names to use
  my %h = %{ $results->primer_results(1) };
  my @fields1 = qw/PRIMER_LEFT_SEQUENCE PRIMER_LEFT_TM PRIMER_RIGHT_SEQUENCE PRIMER_RIGHT_TM/;

  # XXX could use a library for set difference
  my %h1 = map {$_ => 1} @fields1;
  my @fields2 = sort (grep {not $h1{$_}} (keys %h));
  my @fields = (@fields1, @fields2);

  my @a = [ ("gene", "i", @fields) ];

  foreach my $i (0..($results->number_of_results-1)) {
    my %h = %{ $results->primer_results($i) };
    push @a, [ ($name, $i, map { $h{$_}; } @fields) ];
  }

  return @a;
}

# Designs primers for one gene.
# Returns: list of refs. to lists.
#   The first row is headers.
sub design_primers {
  my($gene) = @_;

  my $bed_line = $feature_by_gene{$gene};
  if (not defined $bed_line) {
    return undef;
  }
  my $s = get_region_seq($bed_line, 30);
  $s =~ s/\n//g;
#  my $seq = new Bio::Seq(-id => "pcr_target", -seq => $s, -alphabet => 'dna');

  my $primer3 = Bio::Tools::Run::Primer3->new;
#    -outfile => "primer3.tmp");

  # hack to add the sequence "by hand"
  $primer3->no_param_checks(1);
  $primer3->add_targets('SEQUENCE_TEMPLATE' => $s);

#  print join " ", %{ $primer3->arguments };

  my $results = $primer3->run;

  return format_results($results, $gene);
}


# my @a = keys(%feature_by_gene);
# print $a[123];
my @r = design_primers("Y66C5A.1");

foreach (@r) {
  print join "\t", @{ $_ };
  print "\n";
}

