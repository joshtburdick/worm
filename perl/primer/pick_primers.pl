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

# same, for Tophat2
my $tophat_bwt = $bwt;
my $tophat_gtf = "/home/jburdick/gcb/git/data/seq/genes_tophat_WS220.gtf";
my $tophat_transcriptome_dir = "/var/tmp/data/tophat2/WS220_bt2/";

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
#     (presumably these are adjacent exons)
# Returns: two lists:
# - the sequence, as a list of strings, one per region
# - the junctions to include
sub get_region_seq {
  my($bed_ref) = @_;
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
  my @r = ();
  my @junctions = ();
  my $offset = 0;     # for tracking junctions
  foreach my $i1 (0..($n-1)) {
    # for minus strand, loop through in reverse order
    my $i = $strand eq "+" ? $i1 : (($n-1) - $i1);
    my $a1 = $a + $block_starts[$i];
    my $b1 = $a1 + $block_sizes[$i];
    my $s1 = $strand eq "+" ?
      $genome_fasta->seq($chr, $a1+1 => $b1) :
      $genome_fasta->seq($chr, $b1 => $a1+1);
    push @r, $s1;
    if ($i1 > 0) {
      push @junctions, $offset;
    }
    $offset += length($s1);
  }
  return (\@r, \@junctions);
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

  # FIXME: this should only include adjacent exons common to
  # all known transcripts
  my $bed_line = $feature_by_gene{$gene};
  if (not defined $bed_line) {
    return undef;
  }
  my($exons_ref, $junctions_ref) = get_region_seq($bed_line, 30);
  my $s = join "", @$exons_ref;

  my $primer3 = Bio::Tools::Run::Primer3->new(-outfile => "primer3.tmp");

  # hack to add the sequence "by hand"
  $primer3->no_param_checks(1);
  $primer3->add_targets('SEQUENCE_TEMPLATE' => $s,
    'SEQUENCE_OVERLAP_JUNCTION_LIST' => (join ' ', @$junctions_ref),
    'PRIMER_MIN_TM' => 65,
    'PRIMER_OPT_TM' => 68,
    'PRIMER_MAX_TM' => 71,
    'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION' => 8,
    'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION' => 8,
    'PRIMER_PRODUCT_SIZE_RANGE' => '50-75 75-100 100-125 125-200');

  my $results = $primer3->run;
}

# Runs Bowtie2 on the .fa file, and counts how many times
# each primer matched.
# Side effects: writes out a .bam file
# Returns: ref to a hash, indexed by primer name, of matches
#  of each primer (currently the CIGAR strings.)
sub align_primers {

#  system("bowtie2 $bwt -f $output_dir/primers.fa " .
#    " --very-sensitive " .
#    " -D 20 -R 3 -N 1 -L 5 -i S,1,0 --mp 2,2 " .
#    "| samtools view -bS - | samtools sort - $output_dir/primers");

  # trying Tophat2
  system("tophat2 " .
    " --GTF $tophat_gtf " .
    " --transcriptome-index $tophat_transcriptome_dir/ " .
    " --transcriptome-only --b2-very-sensitive " .
    " --segment-length 10 --min-anchor-length 3 " .
    " --min-coverage-intron 1 " .
    " --output-dir $output_dir/tophat2 " .
    " $tophat_bwt $output_dir/primers.fa ");
  system("samtools index $output_dir/tophat2/accepted_hits.bam");

  # things to keep track of for each alignment
  my %loc = ();
  my %cigar = ();
  open IN, "samtools view $output_dir/tophat2/accepted_hits.bam |" || die;
  while (<IN>) {
    chomp;
    my @a = split /\t/;

    # store the locations
    push @{ $loc{$a[0]} }, $a[3];

    # and the cigar strings
    my $s = $cigar{$a[0]};
    $s = defined $s ? "$s " . $a[5] : $a[5];
    $cigar{$a[0]} = $s;
  }
  close IN;

  my %r = ('loc' => \%loc, 'cigar' => \%cigar);
  return \%r;
}

# Writes out the results.
sub write_summary {
  my($outfile, $header_r, $primer_results, $align_results) = @_;
  my @r = @$primer_results;

  open OUT, ">$outfile" || die;
  print OUT (join "\t", (@$header_r,
    "align_LEFT", "align_RIGHT")) . "\n";
  foreach (@$primer_results) {
    my @a = @$_;

    # get genomic size
    # FIXME: this isn't including the CIGAR string size
    my @align_left =
      @{ $align_results->{'loc'}->{ $a[0] . "_" . $a[1] . "_L"} };
    my @align_right =
      @{ $align_results->{'loc'}->{ $a[0] . "_" . $a[1] . "_R"} };
    my $approx_genomic_size = "";
    if (@align_left == 1 && @align_right == 1) {
      $approx_genomic_size = abs($align_left[0] - $align_right[0]);
    }

    # get number of mismatches for each
    my $cigar_left =
      $align_results->{'cigar'}->{ $a[0] . "_" . $a[1] . "_L"};
    my $cigar_right =
      $align_results->{'cigar'}->{ $a[0] . "_" . $a[1] . "_R"};
    if (not defined $cigar_left) { $cigar_left = ""; }
    if (not defined $cigar_right) { $cigar_right = ""; }

    print OUT (join "\t", (@a,
      $cigar_left, $cigar_right)) . "\n";
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
  next if $p->number_of_results() < 1;

  my($header_r, $r, $fa) = process_results($p, $gene);
  @header = @$header_r;
  print FASTA $fa;
  push @results, @$r;
}

close IN;
close FASTA;

my $align_results = align_primers();

write_summary("$output_dir/primers.tsv",
  \@header, \@results, $align_results);

