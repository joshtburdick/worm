#!/usr/bin/perl -w
# Combines many MEME-format files from many directories
# (e.g. from running MEME) into one large minimal
# MEME-format file, with each motif named according to
# the directory it's in.
# (Arguably, I should have renamed these in the first place.)
# Writes to standard output.

use strict;

my $base_dir = $ARGV[0];
my $prefix = $ARGV[1];

die if not defined $base_dir;
die if not defined $prefix;

my $meme_bin = "/home/jburdick/meme/bin/";


foreach my $f (sort <$base_dir/*>) {
  die if not ($f =~ /\/([^\/]+)$/);
  my $base_name = $1;

  open IN, "$meme_bin/meme2meme $f/meme.txt |" || die;
  while (<IN>) {
    chomp;
    if (/^MOTIF (.*)$/) {
      my $name = $prefix . "_" . $base_name . "_" . $1;
      print "MOTIF $name\n";
    }
    else {
      print "$_\n";
    }

  }

  close IN;
}


