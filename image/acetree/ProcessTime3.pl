#!/usr/bin/perl -w
# Gets time at which each image was acquired.
#
# This only works in cases in which there is a file called
# 'imageList.tsv'.

use strict;

use Date::Parse;
use XML::LibXML;

my $series = $ARGV[0];

# my $series = "/gpfs/fs0/u/azach/images//20140228_cwn-2_L2/";

die if not $series =~ /\/?([^\/]+)\/?$/;
my $series_name = $1;
print "[getting timestamps for series $series_name]\n";

# Gets image names for one series.
sub get_image_names {
  my($series_dir) = @_;

  # get list of image names
  my @image_names = ();
  open IN, "<$series_dir/dats/imageList.tsv" || die;
  $_ = <IN>;
  die if not $_ eq "time\timage\n";
  while (<IN>) {
    chomp;
    my @a = split /\t/;
    push @image_names, $a[1];
  }
  close IN;

  return @image_names;
}

# Gets mapping from series name to date.
sub get_times {
  my($series_dir) = @_;

  # read and parse the XML
  my $parser = XML::LibXML->new();
  local( $/, *FH ) ;
  open FH, "zcat \"$series_dir/dats/omxml.xml.gz\" |" || die;
  my $xml = <FH>;
  close FH;
  my $r = $parser->load_xml({string => $xml});

  # parse out just the Image and AcquisitionDate
  my %when_by_name = ();
  my @image_node = $r->getElementsByTagName("Image");
  foreach (@image_node) {
    my $name = $_->getAttribute("Name");
    my @ad = $_->getElementsByTagName("AcquisitionDate");
    my $when = $ad[0]->to_literal;
    $when_by_name{$name} = $when;
  }

  return \%when_by_name;
}

# Converts from time to relative time.
# Args:
#   names - ref. to a list of names of images
#   when_by_name - ref. to a hash, mapping name to time
# Returns: list of times in seconds, relative to the
#   first plane.
sub relative_time {
  my($names, $when_by_name) = @_;

  my $t0 = str2time($when_by_name->{$names->[0]});

  my @r = map { str2time($when_by_name->{$_}) - $t0 } @$names;

#  return @r;
}

# Writes out the times as a table
sub write_table {
  my($outfile_name, $series_name, $t_ref) = @_;
  my @t = @$t_ref;
  open OUT, ">$outfile_name" || die;

  my $t_prev = 0;
  foreach my $i (0..(@t-1)) {
    print OUT (join "\t",
      $series_name,
      sprintf("%03d", $i+1),
      $t[$i],
      $t[$i] - $t_prev); 
    print OUT "\n";
    $t_prev = $t[$i];
  }

  close OUT;
}

my @image_names = get_image_names($series);

my $when_by_name = get_times($series);

my @t = relative_time(\@image_names, $when_by_name);


write_table($series . "/dats/TIME" . $series_name . ".csv",
  $series_name, \@t);

