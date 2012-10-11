#!/usr/bin/perl -w
# Runs Tophat2 on colorspace reads which are stored in an .xsq file
# (including some gratuitous format conversions.)

use strict;

use File::Temp;

my $csfasta_dir = $ARGV[0];
die "$csfasta_dir doesn't look like xsqconvert's output\n" if not -e "$csfasta_dir/Libraries/";
my $output_dir = $ARGV[1];
die if not defined $output_dir;

# The iGenome build (with colorspace index.)
my $bwt_index = "~/data/seq/Caenorhabditis_elegans/UCSC/ce6/Sequence/ColorspaceBowtieIndex/genome";
my $gene_gtf = "~/data/seq/Caenorhabditis_elegans/UCSC/ce6/Annotation/Genes/genes.gtf";

my $bin_dir = "/murrlab/software/bin";

# this index of the transcriptome should get reused across runs
# my $transcriptome_index = "/var/tmp/Caenorhabditis_elegans_ce6_genes";
my $transcriptome_index = "~/data/seq/Caenorhabditis_elegans/UCSC/ce6/TranscriptomeColorspaceBowtieIndex/genes";

# extract the reads into a temporary directory
# ??? this may not be part of this script...
# my $csfasta_dir = File::Temp->newdir( CLEANUP => undef );
# system("$bin_dir/xsqconvert -c -o $csfasta_dir $xsq_file");
# exit(0);

foreach my $f (<$csfasta_dir/Libraries/*/F3/reads/*.csfasta>) {

  print "$f\n";

  die if not $f =~ /.*Libraries\/([^\/]+)\/F3\/reads\/([^\/]+)_F3\.csfasta/;
  my $library_name = $1;
  my $base_name = $2;

  next if (($library_name eq "Unassigned") || ($library_name eq "Unclassified"));

  print "$library_name $base_name\n";

  my $output_name = $base_name;
  $output_name =~ s/ /_/g;   # avoid spaces in names, as they cause issues

  # XXX rude hack
  next if $output_name eq "Murray050912_Pool1Pool2_02_cehm27";

#  system("mkdir -p $csfasta_dir/modified_csfasta/$output_name");

  rename_csfasta("$csfasta_dir/Libraries/$library_name/", $base_name, $csfasta_dir);
  run_tophat2($csfasta_dir, "$output_dir/$output_name");
}

# Extracts csfasta files, with read name.
sub rename_csfasta {
  my($in_dir, $base_name, $out_dir) = @_;

  system("rm $out_dir/F3.csfasta");
  system("rm $out_dir/F3.qual");
  system("rm $out_dir/F5-RNA.csfasta");
  system("rm $out_dir/F5-RNA.qual");

  system("sed -e 's/_F3\$//' < \"$in_dir/F3/reads/$base_name" .
    "_F3.csfasta\" > $out_dir/F3.csfasta");
  system("sed -e 's/_F3\$//' < \"$in_dir/F3/reads/$base_name" .
    "_F3.QV.qual\" > $out_dir/F3.qual");

  system("sed -e 's/_F5-RNA\$//' < \"$in_dir/F5-RNA/reads/$base_name" .
    "_F5-RNA.csfasta\" > $out_dir/F5-RNA.csfasta");
  system("sed -e 's/_F3\$//' < \"$in_dir/F5-RNA/reads/$base_name" .
    "_F5-RNA.csfasta\" > $out_dir/F5-RNA.qual");
}

# Runs Tophat2 on one library.
# Args:
#   reads_dir - directory containing .csfasta/.qual files (with converted read names)
#   output_dir - directory to write TopHat output into
# Side effects: runs Tophat2.
sub run_tophat2 {
  my($reads_dir, $output_dir) = @_;

#  my $reads_base = "$reads_dir/$library/" . (join '_', ($experiment, $lane, $library));

  system("mkdir -p $output_dir");

# note: not including "--GTF $gene_gtf", as we're using a prebuilt transcriptome index
system(<<END);
$bin_dir/tophat2 --transcriptome-index $transcriptome_index --GTF $gene_gtf --no-novel-juncs --library-type fr-secondstrand --mate-inner-dist 200 --mate-std-dev 200 --color --integer-quals --bowtie1 --num-threads 6 --output-dir $output_dir $bwt_index $reads_dir/F3.csfasta $reads_dir/F5-RNA.csfasta $reads_dir/F3.qual $reads_dir/F5-RNA.qual
END
;

}

