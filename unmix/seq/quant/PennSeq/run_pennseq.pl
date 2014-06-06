#!/usr/bin/perl -w
# Runs MISO on various samples.

use POSIX;

my $genome_fasta = "/var/tmp/data/tophat2/WS220/genome.fa";
my $header = "/var/tmp/data/tophat2/WS220/genome.dict";

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

# Runs MISO on one sample's .bam files.
# Args:
#   file_names_ref - ref. to a list of .bam files
#     (as absolute pathnames)
#   options - options to pass to Cufflinks
#     (FIXME currently ignored)
#   output_path - where to put Cufflinks' output directory
# Side effects: runs MISO
sub run_miso_on_sample {
  my($file_names_ref, $options, $output_path) = @_;

  my $tmp_bam = "/var/tmp/miso_tmp.bam";

  my $input_bam_files = join " ", @{$file_names_ref};

  # merge .bam files into a temporary (reheadered) .bam file
  # XXX "-u" option to "samtools merge" doesn't seem to be working
  if (@{$file_names_ref} > 1) {
    system("samtools merge $tmp_bam $input_bam_files");
  }
  else {
    system("cat $input_bam_files > $tmp_bam");
  }
  system("samtools index $tmp_bam");

  # run MISO
#  system("$cuff_bin --output-dir $output_path -p 7 --GTF $gtf_file " .
#    " --library-type fr-secondstrand --quiet " .
#    " $tmp_bam");
  system("./run_MISO.pl $tmp_bam $output_path");

  unlink($tmp_bam);
  unlink($tmp_bam . ".bai");
}

sub run_miso_on_dir {
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

    run_miso_on_sample(\@files,
      $options,
      "$output_dir/$sample");
  }

}

run_miso_on_dir("/murrlab/seq/tophat2/WS220_20140111/Murray050912/",
  "",
  "~/tmp/miso/Murray_050912");

run_miso_on_dir("/murrlab/seq/tophat2/WS220_20140111/Murray092812/",
  "",
  "~/tmp/miso/Murray_092812");

if (1) {
run_miso_on_dir("/murrlab/seq/tophat2/WS220/20110922/",
  "",
  "~/tmp/miso/Murray_20110922");
}

