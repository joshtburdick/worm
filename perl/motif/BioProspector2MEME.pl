#!/usr/bin/perl -w
# Extracts motifs from BioProspector output, and converts to MEME format.
# Arguments:
#   input_dir - directory containing BioProspector result files to read
#   output_file - file in which to write results
#   prefix - prefix to add onto the output motif names

use strict;

use IO::File;

convert_bp_dir($ARGV[0], $ARGV[1], $ARGV[2]);

# XXX testing
#convert_bp_dir("/home/jburdick/gcb/git/unmix/seq/cluster/hierarchical/hier.ts.50clusters/BioProspector/",
#  "./bp_meme.txt",
#  "foo");

# Converts motifs in a directory of BioProspector results.
sub convert_bp_dir {
  my($bp_dir, $output_file, $prefix) = @_;

  my $of = new IO::File($output_file, "w");

$of->print(<<'END');
MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from /home/jburdick/gcb/git/tf/motif/Ce_WS220.order2markov.txt):
A 0.32700 C 0.17300 G 0.17300 T 0.32700 

END
  ;

  foreach my $f (<$bp_dir/*.txt>) {
    die if not ($f =~ /\/([^\/]+)\.txt$/);
    my $file_short_name = $1;
# print "$file_short_name ";
    convert_motifs_bp_to_meme($f, $of,
      $prefix . "_" . $file_short_name);
  }

  $of->close;
}

# Converts motifs in one BioProspector file.
# Args:
#   bp_file - BioProspector output file to read
#   of - filehandle to write results to
#   prefix - prefix for results
sub convert_motifs_bp_to_meme {
  my($bp_file, $of, $prefix) = @_;
  my $motif_num = undef;
  my $motif_width = undef;
  my $num_sites = undef;

  open IN, "<$bp_file" or die;

  # check of file format
  $_ = <IN>; $_ = <IN>; $_ = <IN>;
  die if not /BioProspector Search Result/;

  while (<IN>) {

    if (/Motif #(\d+):/) {
      $motif_num = $1;
      $_ = <IN>;
      $_ = <IN>;
      die if not /Width \((\d+),.*Sites (\d+)$/;
      $motif_width = $1;
      my $num_sites = $2;
# print "$motif_num $motif_width\n";
      my $name = $prefix . "_" . $motif_num;

      print $of "MOTIF $name\n\n";
# ??? what should E be?
      print $of "letter-probability matrix: alength= 4 w= $motif_width nsites= $num_sites E= 1\n";

      $_ = <IN>; $_ = <IN>;
      foreach my $i (1..$motif_width) {
        $_ = <IN>;
        my @r = split /\s+/;
        print $of join " ", ($r[1] / 100, $r[2] / 100, $r[3] / 100, $r[4] / 100);
        print $of "\n";
      }
      print $of "\n";
    }
  }

  close IN;
}

