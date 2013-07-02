#!/usr/bin/perl -w
# Counts total coverage in some regions (using Bio::DB::BigFile.)
# The file of regions can include duplicated regions, and this
# will print the total in each named region.
# (Presumably this will need normalizing by region size; however
# that's not done here.)

use strict;

# FIXME: make these command-line args
my $regions_file = "/home/jburdick/gcb/git/tf/motif/motifCount/upstreamRegionsWS220_5kb_nogenes_cons_0.5.bed";
# my $regions_file = "/home/jburdick/gcb/git/tf/motif/motifCount/upstreamRegionsWS220_5kb_nogenes.bed";
my $bigwig_dir = "/murrlab/seq/igv/motif/meme/";

# where Jim Kent's tools are located
my $kent_bin = "/home/jburdick/bin/x86_64/";

# loop through the files
my $header_printed = undef;
foreach my $f (<$bigwig_dir/*.bw>) {
  die if not $f =~ /\/(.*)\.bw$/;
  my $name = $1;

  # count the regions (for now, just for one file)
  my $sum_ref = count_totals($regions_file, $f);

  # possibly print the header
  if (!$header_printed) {
    print "\t";
    print join "\t", keys(%$sum_ref);
    $header_printed = 1;
  }

  write_counts($sum_ref, $name);
}

# Counts the totals for one file, using the UCSC
# bigWigAverageOverBed tool.
# Args:
#   regions_file - BED format file of regions
#   bw_file - name of a BigWig-format file
# Returns: a reference to a hash of total counts for each
sub count_totals {
  my($regions_file, $bw_file) = @_;

  # temporary output file
  # XXX note that this isn't multithreaded
  my $tmpfile = "bigWigAverageOverBed_tmp.tsv";

  # compute the averages
  system("$kent_bin/bigWigAverageOverBed $bw_file $regions_file $tmpfile");

  my %sum = ();

  open IN, "<$tmpfile" || die;
  while (<IN>) {
    chomp;
    my($name, $size, $covered, $sum, $mean0, $mean) = split /\t/;
    $name =~ s/:\d+$//;
    
    my $x = $covered;
    if (defined($sum{$name})) {
      $sum{$name} += $x;
    }
    else {
      $sum{$name} = $x;
    }
  }

  close IN;
  unlink($tmpfile);

  return \%sum;
}

# Counts the totals for one file, using Bio::DB::BigFile .
# Args:
#   regions_file - BED format file of regions
#   bw_file - name of a BigWig-format file
# Returns: a reference to a hash of total counts for each
#
# Deprecated; this seems nontrivially slower than using
# bigWigAverageOverBed + post-processing the (admittedly huge)
# output.
sub count_totals_BioPerl {
  my($regions_file, $bw_file) = @_;

  use POSIX;

  use Bio::DB::BigFile;
  use Bio::DB::BigFile::Constants;

  # read the regions into an array
  open R, "<$regions_file" || die;
  my @regions = <R>;
  close R;

  # open the BigWig file
  my $wig = Bio::DB::BigFile->bigWigFileOpen($bw_file);

  # the total counts, indexed by name
  my %sum = ();

  # loop through the regions
  my $count = 0;
  foreach (@regions) {
    chomp;
    my($chr, $a, $b, $name, $score, $strand) = split /\t/;
    $chr =~ s/chr//;

# print "counting in $chr:$a-$b\n";

    # strip off any trailing ":number"
    $name =~ s/:\d+$//;

    # get the coverage count (if there's nothing there, add 0)
    my $r = $wig->bigWigSingleSummary("$chr",$a=>$b,bbiSumCoverage,0);
    my $c = (defined $r ? ($r*($b-$a)) : 0);

    if (defined $sum{$name}) {
      $sum{$name} += $c;
    }
    else {
      $sum{$name} = $c;
    }

    print "$name\t$c\n";
#    last if $count++ > 10000;
  }

  return \%sum;
}

# Writes total counts to standard output.
sub write_counts {
  my($sum_ref, $name) = @_;
  my %h = %$sum_ref;

  my @k = sort (keys %h);
  print $name;

  # print out entries, sorted by the key
  foreach (@k) {
    print "\t";
    print $h{$_};
  }
  print "\n";
}

