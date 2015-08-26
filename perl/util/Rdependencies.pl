#!/usr/bin/perl -w
# Finds which files depend on which other files in
# some R scripts (by default, all R scripts in the
# current directory.) Writes to standard output.
#
# XXX Note that this is based on a pretty superficial
# parsing of the R source files.

use strict;


# Gets dependencies for one file.
# Args:
#   f - the file being parsed
sub get_dependencies {
  my($f) = @_;
  my @deps = ();

  # scan for dependencies
  open IN, "<$f" || die;
  while (<IN>) {
    chomp;
    if (/source\("[^"]+.r"\)/) {
      push @deps, $1;
    }
    if (/load\("[^"]+.Rdata"\)/) {
      push @deps, $1;
    }
  }
  close IN;

  # print 
  if (@deps > 0) {
    print "$f : " . (join " ", @deps);
  }
}

# run this for all .r files in the current directory



