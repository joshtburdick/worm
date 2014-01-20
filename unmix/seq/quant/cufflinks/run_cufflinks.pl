#!/usr/bin/perl -w
# Runs Cufflinks on various samples.

use POSIX;

my $output_dir = "/home/jburdick/tmp/cufflinks/";

my $cuff_bin = "/murrlab/software/bin/cufflinks";

# XXX undo the previous strand-flipping
my $strand_flip_filter = "../../../../perl/seq/strand_flip_filter.pl first";

my $gtf_file = "/var/tmp/data/tophat2/WS220/transcriptome.gff";
# my $gtf_file = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Annotation/Genes/genes.gtf";

# my $header = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/WholeGenomeFasta/genome.dict";
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

# Runs Cufflinks on one sample's .bam files.
# Args:
#   file_names_ref - ref. to a list of .bam files
#     (as absolute pathnames)
#   cufflinks_options - options to pass to Cufflinks
#   output_path - where to put Cufflinks' output directory
#   flip - true if we need to flip the strand of the first read
# Side effects: runs Cufflinks
sub run_cufflinks_on_sample {
  my($file_names_ref, $cufflinks_options, $output_path, $flip) = @_;

  my $tmp_bam = tmpnam() . ".bam";

  my $input_bam_files = join " ", @{$file_names_ref};

  # XXX possibly flip the strand of the first read
  my $flip_cmd = $flip ? $strand_flip_filter : " cat ";

  # merge .bam files into a temporary (reheadered) .bam file
  # XXX "-u" option to "samtools merge" doesn't seem to be working
  if (@{$file_names_ref} > 1) {
    system("samtools merge - $input_bam_files " .
      "| samtools view - | $flip_cmd " .
      "| samtools view -bS -T $genome_fasta - " .
      "| samtools reheader $header - > $tmp_bam");
  }
  else {
    system("samtools view $input_bam_files | $flip_cmd " .
      "| samtools view -bS -T $genome_fasta - " .
      "| samtools reheader $header - > $tmp_bam");
  }

  # run Cufflinks
  system("$cuff_bin --output-dir $output_path -p 7 --GTF $gtf_file " .
    " --library-type fr-secondstrand --quiet " .
    " $tmp_bam");

  # remove the .gff file
  unlink("$output_dir/transcripts.gtf");

  unlink($tmp_bam);
}

sub run_cufflinks_on_dir {
  my($bam_path, $cufflinks_options, $output_dir, $flip) = @_;
  system("mkdir -p $output_dir");

  my @files = <$bam_path/*.bam>;
  my @samples = get_sample_names(\@files);
print @samples;
  foreach my $sample (@samples) {
    print "[running on sample $sample]\n";

    my @files = grep /$sample\.bam$/, @files;
    print "files are " . (join " ", @files);
    print "\n";

    run_cufflinks_on_sample(\@files,
      $cufflinks_options,
      "$output_dir/$sample", $flip);
  }

}

run_cufflinks_on_dir("/murrlab/seq/tophat2/WS220_20140111/Murray050912/",
  " --library-type fr-secondstrand ",
  "~/tmp/cufflinks/Murray_050912", 1);

run_cufflinks_on_dir("/murrlab/seq/tophat2/WS220_20140111/Murray092812/",
  " --library-type fr-secondstrand ",
  "~/tmp/cufflinks/Murray_092812", 1);

if (1) {
run_cufflinks_on_dir("/murrlab/seq/tophat2/WS220/20110922/",
  " --library-type fr-secondstrand ",
  "~/tmp/cufflinks/Murray_20110922", undef);
}

