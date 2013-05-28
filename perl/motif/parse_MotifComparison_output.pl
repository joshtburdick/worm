#!/usr/bin/perl -w
# Parses output of MotifComparison.

use strict;

$_ = <>;
die if not /^#INCLUSive GFF File/;

$_ = <>;
$_ = <>;
$_ = <>;
$_ = <>;
$_ = <>;
die if not /^#query:query_consensus/;
$_ = <>;

# loop through the actual matches
while (<>) {
  chomp;
  my @a = split /\t/;

  # possibly stop
  if (not defined $a[2]) {
    last;
  }

  # ignoring, for now, whether these are significant
  $a[0] =~ s/^##//;

  die if not ($a[0] =~ /([^:]+):([^:]+)/);
  my $query_name = $1;
  my $query_consensus = $2;

  die if not ($a[1] =~ /([^:]+):([^:]+)/);
  my $match_name = $1;
  my $match_consensus = $2;

  die if not ($a[2] =~ /([^;]+);\(([^,]+),([^\)]+)\)$/);
  my $score = $1;
  my $offset1 = $2;
  my $offset2 = $3;

  print join "\t", ($query_name, $query_consensus,
    $match_name, $match_consensus, $score, $offset1, $offset2);
  print "\n";
}

$_ = <>;
die if not /^#SUMMARY/;

