#!/usr/bin/perl -w
# Converts from .xsq to gzipped .csfasta / .qual files.
# XXX the XSQ file names, and output paths, are hard-coded.

use strict;

use File::Temp;

# my $output_dir = $ARGV[0];
# die if not (@ARGV == 1);

# XXX
my $output_dir = "/media/3f447a33-7182-4638-a6dc-16d2fbbe9654/reads/";

my $bin_dir = "/murrlab/software/bin";

# my $csfasta_dir = File::Temp->newdir();
my $csfasta_dir = "/media/disk2/jburdick/csfasta/";
print "writing temporary files to $csfasta_dir\n";

# Tweaks read names in a csfasta file.
sub rename_csfasta {
  my($in_dir, $base_name, $out_dir) = @_;

  system("mkdir -p $out_dir");

#  system("rm $out_dir/F3.csfasta.gz");
#  system("rm $out_dir/F3.qual.gz");
#  system("rm $out_dir/F5-RNA.csfasta.gz");
#  system("rm $out_dir/F5-RNA.qual.gz");

  system("sed -e 's/_F3\$//' < \"$in_dir/F3/reads/$base_name" .
    "_F3.csfasta\" | gzip -c --best > $out_dir/F3.csfasta.gz");
  system("sed -e 's/_F3\$//' < \"$in_dir/F3/reads/$base_name" .
    "_F3.QV.qual\" | gzip -c --best > $out_dir/F3.qual.gz");

  system("sed -e 's/_F5-RNA\$//' < \"$in_dir/F5-RNA/reads/$base_name" .
    "_F5-RNA.csfasta\" | gzip -c --best > $out_dir/F5-RNA.csfasta.gz");
  system("sed -e 's/_F5-RNA\$//' < \"$in_dir/F5-RNA/reads/$base_name" .
    "_F5-RNA.QV.qual\" | gzip -c --best > $out_dir/F5-RNA.qual.gz");
}

# Extracts one .xsq file, does read renaming, and
# compresses the read files.
sub extract_xsq {
  my($xsq_file, $output_dir) = @_;

  print "[extracting $xsq_file]\n";
  system("$bin_dir/xsqconvert -c -o $csfasta_dir $xsq_file");

  foreach my $f (<$csfasta_dir/Libraries/*/F3/reads/*.csfasta>) {

    print "[converting $f]\n";

    die if not $f =~ /.*Libraries\/([^\/]+)\/F3\/reads\/([^\/]+)_F3\.csfasta/;
    my $library_name = $1;
    my $base_name = $2;

    next if (($library_name eq "Unassigned") ||
      ($library_name eq "Unclassified"));

    print "$library_name $base_name\n";

    my $output_name = $base_name;
    $output_name =~ s/ /_/g;   # avoid spaces in names, as they cause issues

    rename_csfasta("$csfasta_dir/Libraries/$library_name/", $base_name,
      "$output_dir/$output_name");
  }

  # clean up, to avoid running out of disk space
  system("\\rm -rf $csfasta_dir/Libraries");
}


sub extract_all {

  my $xsq_prefix = "/murrlab/seq/reads/Murray050912/xsq/Murray050912_Pool1Pool2_0";
#  foreach my $i (1..6) {
#    my $f = $xsq_prefix . $i . ".xsq";
#      print "running on $f\n";
#    extract_xsq($f, $output_dir);
#  }

  my $p = "/murrlab/seq/reads/Murray_52831_092812/result/";
#  extract_xsq("$p/lane1/Murray_52831_092812_01.xsq", $output_dir);
#  extract_xsq("$p/lane2/Murray_52831_092812_02.xsq", $output_dir);
#  extract_xsq("/home/jburdick/xsq/Murray_52831_092812_03.edit.xsq",
#    $output_dir);
#  extract_xsq("/home/jburdick/xsq/Murray_52831_092812_04.edit.xsq",
#    $output_dir);
  extract_xsq("$p/lane5/Murray_52831_092812_05.xsq", $output_dir);
  extract_xsq("$p/lane6/Murray_52831_092812_06.xsq", $output_dir);
}

extract_all();

