#!/usr/bin/perl -w
# Measures transcription which isn't otherwise annotated.
# Requires bedtools.

# where to write output
my $output_base = "otherCoverage";

# files defining "annotated" coverage
my $known_exons = "../../../../git/data/seq/exons_WS220.bed";
my $known_gene_bounds = "../../../../git/data/seq/merged_genes_WS220.bed"; 

# Gets the names of all the samples in a directory, assuming
# they're named in a standard way.
sub get_sample_names {
  my($files_ref) = @_;
  my %h = ();
  foreach my $f (@$files_ref) {
    my $sample = $f;

    # just consider base pathname
    $sample =~ /\/([^\/]+)\.bam/;
    $sample = $1;
    $sample =~ s/0\d_//;
    $h{$sample} = 1;
  }

  return keys(%h);
}

# Finds regions covered by reads, which aren't included in a set of
# other regions.
# Args:
#   bam_files_ref - reference to a list of .bam files to combine
#   regions_bed - name of a file
#   output_bed - where to write output
# Side effects: writes out a .bed file of unannotated regions.
sub write_unannotated_coverage_1 {
  my($bam_files_ref, $regions_bed, $output_bed) = @_;
  my @bam_files = @$bam_files_ref;

  # command to use to read from files, or a single file
  my $view_cmd = (@bam_files > 1) ?
    "samtools merge - " . (join ' ', @bam_files) :
    "cat " . $bam_files[0];

  # XXX bam->bed conversion may be slightly slow
  system($view_cmd . " | bedtools bamtobed -split -i - " .
    " | bedtools subtract -s -A -a - -b $regions_bed " .
    " | bedtools merge -s -n > $output_bed");
}

# Writes out unannotated coverage for all files in a particular directory.
# Args:
#   bam_path - directory containing .bam files
#   output_path - directory in which to write output
# Side effects: writes out unannotated coverage
sub write_unannotated_coverage {
  my($bam_path, $output_path) = @_;

  system("mkdir -p $output_path/outside_exons");
  system("mkdir -p $output_path/outside_genes");

  my @files = <$bam_path/*.bam>;
  my @samples = get_sample_names(\@files);
print @samples;
  foreach my $sample (@samples) {
    print "[running on sample $sample]\n";

    my @files = grep /$sample\.bam$/, @files;
    print "files are " . (join " ", @files) . "\n";

    write_unannotated_coverage_1(\@files, $known_exons,
      "$output_path/outside_exons/$sample.bed");
    write_unannotated_coverage_1(\@files, $known_gene_bounds,
      "$output_path/outside_genes/$sample.bed");
  }
}

my $bam_base = "/murrlab/seq/tophat2/WS220_20140111/";

write_unannotated_coverage("$bam_base/Murray050912/", $output_base);
write_unannotated_coverage("$bam_base/Murray092812/", $output_base);
write_unannotated_coverage("$bam_base/20110922/", $output_base);

