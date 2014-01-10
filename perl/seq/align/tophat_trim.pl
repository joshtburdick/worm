#!/usr/bin/perl -w
# Runs Tophat2 with reads trimmed iteratively.

use strict;

my $tophat2_dir = "/murrlab/software/tophat-2.0.10.Linux_x86_64/";

# die if not (@ARGV == 2);

# base name of directory containing F3 and F5-RNA files
my $input_basename = $ARGV[0];

# where to write output to
my $output_dir = $ARGV[1];

# genome and transcriptome indices
my $bwt_base_dir = "/var/tmp/data/tophat2/WS220/";
# my $bwt_base_dir = "/var/tmp/data/tophat2/WS220/";
my $bwt_index = "$bwt_base_dir/genome";
my $transcriptome_index = "$bwt_base_dir/transcriptome";

# clean up files afterwards?
my $remove_tmp = undef;

# options to pass to Tophat2 in all cases
my $tophat2_common_options = " --color --quals --integer-quals --bowtie1 " .
  "--transcriptome-index $transcriptome_index " .
  "--GTF $transcriptome_index.gff " .
  "--min-intron-length 35 --max-intron-length 5000 " .
  "--no-coverage-search " .
  "--num-threads 7 ";


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
  my($bam_file, $f3_size, $f5_size, $in_dir, $out_dir) = @_;

  my $reads_ref = get_read_names($bam_file);
  system("mkdir -p $out_dir");

  filter_reads_clipping_1($reads_ref, $f3_size,
    "$in_dir/F3.csfasta", "$in_dir/F3.qual",
    "$out_dir/F3.csfasta", "$out_dir/F3.qual");
  filter_reads_clipping_1($reads_ref, $f5_size,
    "$in_dir/F5-RNA.csfasta", "$in_dir/F5-RNA.qual",
    "$out_dir/F5-RNA.csfasta", "$out_dir/F5-RNA.qual");
}

# Runs Tophat2 with some set of options.
sub run_tophat {
  my($input_base, $output_dir, $options) = @_;

  my $cmd = "$tophat2_dir/tophat2 $tophat2_common_options " .
    "$options --output-dir $output_dir " .
    "$bwt_index " .
    "$input_base/F3.csfasta $input_base/F5-RNA.csfasta " .
    "$input_base/F3.qual $input_base/F5-RNA.qual ";
  print $cmd;
  system($cmd);
}

# Runs Tophat with various numbers of bases clipped.
sub tophat_tune_clipping {
  my $input_base = "/var/tmp/data/reads/03_F21D5.9_2e6/";
  my $reads_tmp = "/var/tmp/tophat_clip_tune/reads";
  system("mkdir -p $reads_tmp");

  foreach my $clip (0,14,12,10,8) {
    my $f3_clip = $clip;
    my $f5_clip = $clip;

    my $output_dir = "/var/tmp/tophat_clip_tune/clip_" .
      $f3_clip . "_" . $f5_clip;
  
    filter_reads_clipping(undef, 50-$f3_clip, 35-$f5_clip,
      $input_base, $reads_tmp);

    run_tophat($reads_tmp, $output_dir,
      "--transcriptome-only --min-anchor-length 3 ");
  }
}

# Does the mapping for a pair of files.
# Args:
#   input_base - directory containing F3 and F5-RNA reads
#   output_base - base directory in which to write output
sub run_tophat_clipping {
  my($input_base, $output_base) = @_;

  system("mkdir -p $output_base");

  # slower, hopefully thorough version
#  run_tophat($input_base, "$output_base/transcriptome1",
#    "--no-mixed --segment-length 26 --min-anchor-length 3 " . 
#    "--read-edit-dist 6 --read-mismatches 4 ");

# faster version of same
#  run_tophat($input_base, "$output_base/transcriptome1",
#    "--transcriptome-only --no-mixed --min-anchor-length 3 --read-edit-dist 6 --read-mismatches 4 ");

  filter_reads_clipping(
    undef,   # "$output_base/transcriptome1/accepted_hits.bam",
    40, 29, $input_base, "$output_base/reads1");

  run_tophat("$output_base/reads1",
    "$output_base/remainder",
    "--no-novel-juncs --segment-length 23 --min-anchor-length 3",
    "--read-edit-dist 6 --read-mismatches 4 ");
}


# tophat_tune_clipping();

# for testing
run_tophat_clipping("/var/tmp/data/reads/03_F21D5.9_1e5/",
  "/var/tmp/tophat_trim_output/6");


