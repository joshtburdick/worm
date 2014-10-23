#!/usr/bin/perl -w
# Picks primers for a set of genes.
# Args:
#   name of a .tsv file (whose name must end in ".tsv"), containing
#     a list of genes (the first line is assumed a header, and skipped)
# Output: creates a folder, containing:
#   primers.tsv
#   primers.fa - the primers, in FASTA format
#   bowtie.bam - the

use strict;

# use Bio::DB::BigFile::Constants;
use Bio::DB::Fasta;

use Bio::Seq;
use Bio::Tools::GFF;
use Bio::Tools::Run::Primer3;

# Configuration: various sources of data

# gene locations
# FIXME change this to "regions shared by all transcripts"?
# my $gene_gff = "/home/jburdick/gcb/git/data/seq/genes_tophat_WS220.gtf";
my $gene_bed_file = "/home/jburdick/gcb/git/data/seq/genes_WS220.bed";

# set this to wherever the .fa files of the genome are
my $genome_fasta = Bio::DB::Fasta->new(
  "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa");

# Bowtie2 index of genome (for checking uniqueness by mapping them)
my $bwt = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/Bowtie2Index/genome";

# end config section

# create output directory
my $in_file = $ARGV[0];
die if not ($in_file =~ /^(.*)\.tsv$/);
my $output_dir = $1;
system("mkdir -p $output_dir/");

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

# Processes the primer results.
# Returns:
#   ref to list of field names
#   ref to list of refs to lists of fields
#   the primers, in FASTA format
sub process_results {
  my($results, $name) = @_;

  # field names to use
  my %h = %{ $results->primer_results(1) };
  my @fields1 = qw/PRIMER_LEFT_SEQUENCE PRIMER_LEFT_TM PRIMER_RIGHT_SEQUENCE PRIMER_RIGHT_TM/;

  # XXX could use a library for set difference
  my %h1 = map {$_ => 1} @fields1;
  my @fields2 = sort (grep {not $h1{$_}} (keys %h));
  my @fields = (@fields1, @fields2);

  my @a = ();
  my $primer_fasta = "";

  foreach my $i (0..($results->number_of_results-1)) {
    my %h = %{ $results->primer_results($i) };
    push @a, [ ($name, $i, map { $h{$_}; } @fields) ];
    $primer_fasta = $primer_fasta .
      ">" . $name . "_" . $i . "_" . "L\n" . $h{"PRIMER_LEFT_SEQUENCE"} . "\n" .
      ">" . $name . "_" . $i . "_" . "R\n" . $h{"PRIMER_RIGHT_SEQUENCE"} . "\n";
  }

  return (["gene", "i", @fields], \@a, $primer_fasta);
}

# Designs primers for one gene.
# Returns the primer results object.
sub design_primers_1 {
  my($gene) = @_;
  print "$gene\n";

  my $bed_line = $feature_by_gene{$gene};
  if (not defined $bed_line) {
    return undef;
  }
  my $s = get_region_seq($bed_line, 30);
  $s =~ s/\n//g;

  my $primer3 = Bio::Tools::Run::Primer3->new;

  # hack to add the sequence "by hand"
  $primer3->no_param_checks(1);
  $primer3->add_targets('SEQUENCE_TEMPLATE' => $s);

  my $results = $primer3->run;
}

# Runs Bowtie on the .fa file, and counts how many times
# each primer matched.
# Side effects: writes out a .bam file
# Returns: ref to a hash, indexed by primer name, of matches
#  of each primer (currently the CIGAR strings.)
sub run_bowtie2 {

  system("bowtie2 $bwt -f $output_dir/primers.fa " .
    " --very-sensitive " .
#    " -D 20 -R 3 -N 1 -L 5 -i S,1,0.50 " .
    "| samtools view -bS - | samtools sort - $output_dir/primers");
  system("samtools index $output_dir/primers.bam");

  my %h = ();
  open IN, "samtools view $output_dir/primers.bam |" || die;
  while (<IN>) {
    chomp;
    my @a = split /\t/;
    my $s = $h{$a[0]};
    $s = defined $s ? "$s " . $a[5] : $a[5];
    $h{$a[0]} = $s;
  }
  close IN;

  return \%h;
}

# Writes out the results.
sub write_summary {
  my($outfile, $header_r, $primer_results, $bowtie_results) = @_;
  my @r = @$primer_results;

  open OUT, ">$outfile" || die;
  print OUT (join "\t", (@$header_r, "Bowtie_LEFT", "Bowtie_RIGHT")) . "\n";
  foreach (@$primer_results) {
    my @a = @$_;
    my $bowtie_left = $bowtie_results->{ $a[0] . "_" . $a[1] . "_L"};
    my $bowtie_right = $bowtie_results->{ $a[0] . "_" . $a[1] . "_R"};
    print OUT (join "\t", (@a, $bowtie_left, $bowtie_right)) . "\n";
  }
  close OUT;
}

if (undef) {
# my @a = keys(%feature_by_gene);
# print $a[123];
my @r = design_primers("Y66C5A.1");

foreach (@r) {
  print join "\t", @{ $_ };
  print "\n";
}
}

open IN, "<$in_file" || die;
$_ = <IN>;
open FASTA, ">$output_dir/primers.fa" || die;
my @results = ();
my @header = ();
while (<IN>) {
  chomp;
  my $gene = $_;

  my $p = design_primers_1($gene);
  next if not defined $p;
  my($header_r, $r, $fa) = process_results($p, $gene);
  @header = @$header_r;
  print FASTA $fa;
  push @results, @$r;
}

close IN;
close FASTA;

my $bowtie_results = run_bowtie2();

write_summary("$output_dir/primers.tsv",
  \@header, \@results, $bowtie_results);

