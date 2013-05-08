#!/usr/bin/perl -w
# Runs Tomtom on a pair of meme-format motif files.

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

my $meme_bin = "~/meme/bin";
my $bg_file = "/home/jburdick/gcb/git/tf/motif/Ce_WS220.order1markov.txt";

# ~/meme/bin/tomtom -bfile /home/jburdick/gcb/git/tf/motif/Ce_WS220.order1markov.txt -dist pearson  meme.hier_ts_200.txt meme.hier_ts_200.txt

system("$meme_bin/tomtom -bfile $bg_file -dist pearson -verbosity 1 -text $file1 $file2");

