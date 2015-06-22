#!/usr/bin/perl -w
# Converts a .cdt-format TSV file to JSON format.

use strict;
use JSON;

# get array names
$_ = <>;
my @a = split /\t/;
my @arrayName = splice @a, 4;

# skip rest of header
$_ = <>;
$_ = <>;

# other things to keep track of: gene names, descriptions...
my @geneName = ();
my @description = ();

# gene name to number index
# (??? slightly wasteful of space to include this here; possibly
# the script should re-create it instead?)
my %geneNameToIndex = ();

# expression data, as list of refs to lists
my @x = ();

# read in data
my $i = 0;
while (<>) {
  chomp;
  my @a = split /\t/;

  # convert strings that look like numbers to numbers
  # XXX which columns are numbers shouldn't be hard-coded
  foreach my $j (4..(@a-1)) {
    $a[$j] = $a[$j] + 0;
  }

  push @geneName, $a[1];
  push @description, $a[2];
  $geneNameToIndex{$a[1]} = $i;

  @a = splice(@a, 4);
  push @x, \@a;

  $i++;
}

my %h = ('arrayName' => \@arrayName,
  'geneName' => \@geneName,
  'descr' => \@description,
  'x' => \@x,
  'geneNameToIndex' => \%geneNameToIndex);

# print out @r, encoded as JSON
print "a = " . encode_json(\%h) . ";\n";

