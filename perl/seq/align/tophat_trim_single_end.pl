#!/usr/bin/perl -w
# Version of tophat_trim.pl, for use with single-end data.

use strict;

my $tophat2_dir = "/murrlab/software/tophat-2.0.10.Linux_x86_64/";

# die if not (@ARGV == 2);

# genome and transcriptome indices
my $bwt_base_dir = "/var/tmp/data/tophat2/WS220/";
my $bwt_index = "$bwt_base_dir/genome";
my $transcriptome_index = "$bwt_base_dir/transcriptome";

# clean up files afterwards?
my $remove_tmp = 1;

# options to pass to Tophat2 in all cases
my $tophat2_common_options = " --color --quals --integer-quals --bowtie1 " .
  "--transcriptome-index $transcriptome_index " .
  "--min-intron-length 35 --max-intron-length 5000 " .
  "--no-coverage-search " .
  "--num-threads 6 ";

# Gets read names from a .bam file.
# Args:
#   bam_file - file to get read names from
# Returns: ref. to a hash, with read names as key (and value 1)
sub get_read_names {
  my($bam_file) = @_;

  # if file name is undefined, return an empty hash
  if (!defined($bam_file)) {
    my %h = ();
    return \%h;
  }

  my %reads = ();
  open IN, "samtools view $bam_file | cut -f 1 |" || die;
  while (<IN>) {
    chomp;
    $reads{$_} = 1;
  }
  close IN;

  return \%reads;
}

# Utility to filter out some reads from a pair of
#   .csfasta / .qual files.
# Args:
#   reads_ref - reference to a hash of reads
#   n - number of bases which should be present in clipped files
#   in_f, in_q - names of input csfasta and quality files
#   out_f, out_q - names of output csfasta and quality files
sub filter_reads_clipping_1 {
  my($reads_ref, $n, $in_f, $in_q, $out_f, $out_q) = @_;
  my %reads = %{ $reads_ref };

  # number of characters to include in .csfasta
  my $n1 = $n + 1;

  open INF, "<$in_f" || die;
  open INQ, "<$in_q" || die;
  open OUTF, "| cut -c1-$n1 > $out_f" || die;
  open OUTQ, "| cut -f1-$n -d' ' > $out_q" || die;

  # the .csfasta file
#  $_ = <INF>;
  while (<INF>) {
    next if /^#/;
    die if not />(.*)\n$/;
    if (!$reads{$1}) {
      print OUTF $_;
      $_ = <INF>;
      print OUTF $_;
    }
    else {
      $_ = <INF>;
    }
  }

  # the .qual file
  while (<INQ>) {
    next if /^#/;
    die if not />(.*)$/;
    if (!$reads{$1}) {
      print OUTQ $_;
      $_ = <INQ>;
      print OUTQ $_;
    }
    else {
      $_ = <INQ>;
    }
  }

  close INF;
  close INQ;
  close OUTF;
  close OUTQ;
}

# Does that filtering for a set of read files.
sub filter_reads_clipping {
  my($bam_file, $f3_size, $in_dir, $out_dir) = @_;

  my $reads_ref = get_read_names($bam_file);
  system("mkdir -p $out_dir");

  filter_reads_clipping_1($reads_ref, $f3_size,
    "$in_dir/F3.csfasta", "$in_dir/F3.qual",
    "$out_dir/F3.csfasta", "$out_dir/F3.qual");
}

# Runs Tophat2 with some set of options.
sub run_tophat {
  my($input_base, $output_dir, $options) = @_;

  my $cmd = "$tophat2_dir/tophat2 $tophat2_common_options " .
    "$options --output-dir $output_dir " .
    "$bwt_index " .
    "$input_base/F3.csfasta $input_base/F3.qual ";
  print $cmd;
  system($cmd);
}

# Does the mapping for a pair of files.
# Args:
#   input_base - directory containing F3 and F5-RNA reads
#   output_base - base directory in which to write output
sub run_tophat_clipping {
  my($input_base, $output_base) = @_;

  system("mkdir -p $output_base");

  # XXX hack
  system("mkdir -p $output_base/reads");
  system("ln -s $input_base.csfasta $output_base/reads/F3.csfasta");
  system("ln -s $input_base.qual $output_base/reads/F3.qual");

  # slower, hopefully thorough version
  run_tophat("$output_base/reads", "$output_base/transcriptome1",
    "--no-novel-juncs " .
    "--min-anchor-length 4 --splice-mismatch 1 --segment-length 17 " . 
    "--read-edit-dist 6 --read-mismatches 4 --no-sort-bam ");

  filter_reads_clipping(
    "$output_base/transcriptome1/accepted_hits.bam",
    28, "$output_base/reads", "$output_base/reads1");

  run_tophat("$output_base/reads1",
    "$output_base/remainder",
    "--no-novel-juncs --min-anchor-length 4 --splice-mismatch 1 " .
    "--no-sort-bam ");

  # sort, and merge .bam output
  system("samtools sort -m 2500000000 " .
    "$output_base/transcriptome1/accepted_hits.bam " .
    "$output_base/transcriptome1/sorted");
  system("samtools sort -m 2500000000 " .
    "$output_base/remainder/accepted_hits.bam " .
    "$output_base/remainder/sorted");
  system("samtools merge $output_base/merged.bam " .
    "$output_base/transcriptome1/sorted.bam " .
    "$output_base/remainder/sorted.bam ");

  # clean up these files
  if ($remove_tmp) {
    unlink("$output_base/transcriptome1/accepted_hits.bam");
    unlink("$output_base/transcriptome1/sorted.bam");
    unlink("$output_base/remainder/accepted_hits.bam");
    unlink("$output_base/remainder/sorted.bam");
  }
}

# Runs this on all files in the input directory.
die if not (@ARGV == 1);
my $input_dir = "/home/jburdick/data/seq/reads/20110922";
my $output_dir = $ARGV[0];

system("mkdir -p $output_dir");

foreach my $f (<$input_dir/*.csfasta>) {
  my $file_base = $f;
  $file_base =~ s/\.csfasta//;
  die if not ($file_base =~ /\/([^\/]+)$/);
  my $base_name = $1;
  print "running on $file_base, writing to $output_dir/$base_name\n";
  run_tophat_clipping($file_base, "$output_dir/$base_name");
}

