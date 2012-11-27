#!/usr/bin/perl -w
# Imports a directory's worth of .wig.gz files into the bigWig format.

use strict;

# input and output directories
# my $wig_gz_path = "/home/jburdick/data/modencode/ftp/data.modencode.org/C.elegans";
# my $bigwig_path = "/home/jburdick/data/modencode/bigWig/C.elegans";
my $wig_gz_path = "/murrlab/jburdick/src/tf/meme/";
my $bigwig_path = "/murrlab/seq/igv/motif/";

# directory to convert, relative to the above
my $data_path = "meme/";

# needed for bedGraphToBigWig, obtainable using, e.g., "samtools faidx"
my $chromosome_sizes_file = "/murrlab/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.sizes";

import_dir("$wig_gz_path/$data_path", "$bigwig_path/$data_path");

# Writes out a file, just keeping the maximum score for
# each location.
# Args:
#   wig_gz_file - the .wig.gz file to read
sub take_max {
  my($wig_gz_file, $out_file) = @_;

  # hash storing maximum for each location
  my %m = ();

  # first, get maximum score for each location
  open IN, "<$wig_gz_file" || die;
  while (<IN>) {
    chomp;
    my($chr, $a, $b, $score) = split /\t/;
    my $k = "$chr:$a";
    if (!defined $m{$k} || $score >= $m{$k}) {
      $m{$k} = $score;
    }
  }
  close IN;

  # then, write out maximum score for each location
  # (making sure to write out each only once)
  open IN, "<$wig_gz_file" || die;
  open OUT, ">$out_file" || die;
  while (<IN>) {
    chomp;
    my($chr, $a, $b, $score) = split /\t/;
    my $k = "$chr:$a";
    if (defined $m{$k}) {
      print OUT (join "\t", ($chr, $a, $b, $m{k}));
      print OUT;
      $m{$k} = undef;    # to avoid writing something out twice
  }
  close IN;
  close OUT;
}

# Imports all the files in one directory.
sub import_dir {
  my($wig_gz_dir, $bigwig_dir) = @_;

  system("mkdir -p $bigwig_dir");

  foreach my $f (`ls $wig_gz_dir`) {
    chomp $f;
    next if not $f =~ /\.wig\.gz$/;

    print STDERR "converting $f\n";

    # name of file to write
    my $bw_file = $f;
    $bw_file =~ s/\.wig\.gz$/\.bw/;

    # we need a temporary file, since the UCSC tools don't read STDIN, and we
    # need to remove the header line
    my $tmp_file = $f;
    $tmp_file =~ s/\.wig\.gz/\.tmp.bedGraph/;

    system("gunzip -c $wig_gz_dir/$f | fgrep -v track > $bigwig_dir/$tmp_file");

    system("bedGraphToBigWig $bigwig_dir/$tmp_file $chromosome_sizes_file $bigwig_dir/$bw_file");

    unlink("$bigwig_dir/$tmp_file");
  }

}

