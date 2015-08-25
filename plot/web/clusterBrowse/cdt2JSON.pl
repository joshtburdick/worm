#!/usr/bin/perl -w
# Converts a .cdt-format TSV file to JSON format.

use strict;
use JSON;

open IN, "<../../../cluster//hierarchical/hier.300.clusters/cluster.cdt" || die;

# get array names
$_ = <IN>;
my @a = split /\t/;
my @arrayName = splice @a, 9;

# skip rest of header
$_ = <IN>;
$_ = <IN>;

# other things to keep track of: gene names, descriptions, cluster...
my @geneName = ();
my @description = ();
my @cluster = ();

# gene name to number index
# (??? slightly wasteful of space to include this here; possibly
# the script should re-create it instead?)
my %geneNameToIndex = ();

# expression data, as list of refs to lists
my @x = ();

# read in data
my $i = 0;
while (<IN>) {
  chomp;
  my @a = split /\t/;

  # convert strings that look like numbers to numbers
  # XXX which columns are numbers shouldn't be hard-coded
  foreach my $j (9..(@a-1)) {
    if ($a[$j] =~ /\d+/) {
      $a[$j] = $a[$j] + 0;
    }
    else {
      # these should be greyed out
      # XXX the JSON module won't print out NaN;
      # we'll hack around that later
      $a[$j] = "__NaN__";
    }
  }

  push @geneName, $a[2];
  push @description, $a[3];
  push @cluster, $a[5];
  $geneNameToIndex{$a[2]} = $i;

  @a = splice(@a, 9);
  push @x, \@a;

  $i++;
}

my %h = ('arrayName' => \@arrayName,
  'geneName' => \@geneName,
  'descr' => \@description,
  'cluster' => \@cluster,
  'x' => \@x,
  'geneNameToIndex' => \%geneNameToIndex);

# print out @r, encoded as JSON
my $s = encode_json(\%h);
$s =~ s/\"__NaN__\"/NaN/g;
print "a = " . $s . ";\n";

