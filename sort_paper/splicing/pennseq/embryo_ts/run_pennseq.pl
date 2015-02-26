#!/usr/bin/perl -w
# Runs Pennseq on the embryonic timeseries, for one gene.

use POSIX;

# my $pennseq = "/home/jburdick/gcb/software/pennseq/pennseq.pl";
my $pennseq = "./pennseq_1chrom.pl";

# Gets the names of all the samples in a directory.
sub get_sample_names {
  my($files_ref) = @_;
  my %h = ();
  foreach my $f (@$files_ref) {

    # extract the sample time
    die $f if not $f =~ /\/(EE_50-\d+_modENCODE).*\.bam/;
    $sample = $1;
    $h{$sample} = 1;
  }

  return keys(%h);
}

# Runs Pennseq on one sample's .bam files.
# Args:
#   file_names_ref - ref. to a list of .bam files
#     (as absolute pathnames)
#   output_file - where to write output
# Side effects: runs Pennseq on one sample.
sub run_on_sample {
  my($file_names_ref, $output_file) = @_;
  if (@{$file_names_ref} <= 1) {
    die "expected to be merging multiple files\n";
  }

#  my $tmp_sam = "/media/jburdick/disk2/tmp/merge.sam";

  my $input_bam_files = join " ", @{$file_names_ref};

  # merge .bam files, and run Pennseq
  system("samtools merge - $input_bam_files | samtools view - " .
    "| $pennseq -i lin-3_WS220_preprocessed.txt " .
    " -o $output_file -c IV");

#  unlink($tmp_sam);
}

sub run_on_dir {
  my($bam_path, $output_dir) = @_;
  system("mkdir -p $output_dir");

  my @files = <$bam_path/*.bam>;
  my @samples = get_sample_names(\@files);
print @samples;
  foreach my $sample (@samples) {
    print "[running on sample $sample]\n";

    my @files = grep /$sample/, @files;
    print "files are " . (join " ", @files);
    print "\n";

    run_on_sample(\@files, "$output_dir/$sample.txt");
  }

}

system("mkdir -p output");
run_on_dir("/murrlab/seq/tophat2/embryo_ts/", "output/");

