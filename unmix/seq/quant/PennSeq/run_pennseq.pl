#!/usr/bin/perl -w
# Runs PennSeq on various samples.

use POSIX;

my $genome_fasta = "/var/tmp/data/tophat2/WS220/genome.fa";
my $header = "/var/tmp/data/tophat2/WS220/genome.dict";

my $isoform_compatible_matrix_file = "/home/jburdick/gcb/git/sort_paper/splicing/pennseq/refGene_WS220_preprocessed.txt";

# Gets the names of all the samples in a directory.
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

# Runs PennSeq on one sample's .bam files.
# Args:
#   file_names_ref - ref. to a list of .bam files
#     (as absolute pathnames)
#   options - options to pass to Cufflinks
#     (FIXME currently ignored)
#   output_base - where to put the output file
# Side effects: runs PennSeq
sub run_on_sample {
  my($file_names_ref, $options, $output_base) = @_;

  my $input_bam_files = join ",", @{$file_names_ref};

  # run PennSeq; note that we don't have to merge the files
  system("./pennseq_bam.pl -b $input_bam_files -i $isoform_compatible_matrix_file -o $output_base.tsv");
}

sub run_on_dir {
  my($bam_path, $options, $output_dir) = @_;

  system("mkdir -p $output_dir");

  my @files = <$bam_path/*.bam>;
  my @samples = get_sample_names(\@files);
print @samples;
  foreach my $sample (@samples) {
    print "[running on sample $sample]\n";

    my @files = grep /$sample\.bam$/, @files;
    print "files are " . (join " ", @files);
    print "\n";

    run_on_sample(\@files,
      $options,
      "$output_dir/$sample");
  }

}

run_on_dir("/murrlab/seq/tophat2/WS220_20140111/Murray050912/",
  "",
  "~/tmp/pennseq/Murray_050912");

run_on_dir("/murrlab/seq/tophat2/WS220_20140111/Murray092812/",
  "",
  "~/tmp/pennseq/Murray_092812");

if (1) {
run_on_dir("/murrlab/seq/tophat2/WS220/20110922/",
  "",
  "~/tmp/pennseq/Murray_20110922");
}

