#!/usr/bin/perl -w
# Utility to count motifs in a particular region.


my $mode = $ARGV[0];
die if not ($mode eq "--list" || $mode eq "--count");

# name of the motif (e.g. Otx1_DBD_2)
my $motif_name = $ARGV[1];

# the motif score, e.g. 40; this is typically >= 30
my $motif_score = $ARGV[2];

# the motif conservation level (e.g. 0.9)
my $cons_level = $ARGV[3];

# the upstream distance (in kb, e.g. 1; this has a 500 bp minimum)
my $upstream_dist = $ARGV[4];

# lists motifs in BAM format
if ($mode eq "--list") {
  system("samtools view /murrlab/seq/igv/motif/known/" . $motif_name . ".bam " .
  "-q " . $motif_score .
  "-L ~/gcb/git/perl/motif/tmp/phastCons7way_" . $cons_level . ".bed " .
  " ");
}

if ($mode eq "--count") {
  my $cmd = ("samtools view /murrlab/seq/igv/motif/known/" . $motif_name . ".bam " .
  " -q " . $motif_score .
  " -L ~/gcb/git/perl/motif/tmp/phastCons7way_" . $cons_level . ".bed " .
  " -b | bedtools coverage -abam - " .
  " -b ~/gcb/git/tf/motif/motifCount/regions/upstream_" . $upstream_dist . "kb_WS220.bed | cut -f4,7 ");
#  print($cmd);
  system($cmd);
}


