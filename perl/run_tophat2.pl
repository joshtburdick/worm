#!/usr/bin/perl -w
# Runs Tophat2 on colorspace reads which are stored in an .xsq file
# (including some gratuitous format conversions.)

use strict;

use File::Temp;

my $xsq_file = $ARGV[0];

# directory to use for temporary files
my $csfasta_dir = $ARGV[1];

my $output_dir = $ARGV[2];

die if not (@ARGV + 0 == 3);

print "[running on xsq file $xsq_file, writing to $output_dir]\n";

# genome and transcriptome indices
my $bwt_index = "/var/tmp/data/tophat2/WS220/genome";
my $transcriptome_index = "/var/tmp/data/tophat2/WS220/transcriptome";

my $bin_dir = "/murrlab/software/bin";

# extract the reads into a temporary directory
# my $csfasta_dir = File::Temp->newdir( CLEANUP => undef );
system("mkdir -p $csfasta_dir");
system("$bin_dir/xsqconvert -c -o $csfasta_dir $xsq_file");

die "$csfasta_dir doesn't look like xsqconvert's output\n" if not -e "$csfasta_dir/Libraries/";

foreach my $f (<$csfasta_dir/Libraries/*/F3/reads/*.csfasta>) {

  print "$f\n";

  die if not $f =~ /.*Libraries\/([^\/]+)\/F3\/reads\/([^\/]+)_F3\.csfasta/;
  my $library_name = $1;
  my $base_name = $2;

  next if (($library_name eq "Unassigned") ||
    ($library_name eq "Unclassified"));

  print "$library_name $base_name\n";

  my $output_name = $base_name;
  $output_name =~ s/ /_/g;   # avoid spaces in names, as they cause issues

  rename_csfasta("$csfasta_dir/Libraries/$library_name/", $base_name, $csfasta_dir);
#   run_tophat2($csfasta_dir, "$output_dir/$output_name");
  system("/murrlab/jburdick/gcb/git/perl/seq/align/tophat_trim.pl " .
    " $csfasta_dir $output_dir/$output_name ");

}

# removing csfasta directory; otherwise we run out of disk space
# (and potentially re-map a bunch of reads)
system("rm -rf $csfasta_dir");

# Extracts csfasta files, with read name.
sub rename_csfasta {
  my($in_dir, $base_name, $out_dir) = @_;

  system("rm $out_dir/F3.csfasta");
  system("rm $out_dir/F3.qual");
  system("rm $out_dir/F5-RNA.csfasta");
  system("rm $out_dir/F5-RNA.qual");

  system("sed -e 's/_F3\$//' < \"$in_dir/F3/reads/$base_name" .
    "_F3.csfasta\" > $out_dir/F3.csfasta");
  system("sed -e 's/_F3\$//' < \"$in_dir/F3/reads/$base_name" .
    "_F3.QV.qual\" > $out_dir/F3.qual");

  system("sed -e 's/_F5-RNA\$//' < \"$in_dir/F5-RNA/reads/$base_name" .
    "_F5-RNA.csfasta\" > $out_dir/F5-RNA.csfasta");
  system("sed -e 's/_F5-RNA\$//' < \"$in_dir/F5-RNA/reads/$base_name" .
    "_F5-RNA.QV.qual\" > $out_dir/F5-RNA.qual");
}

# Runs Tophat2 on one library.
# Args:
#   reads_dir - directory containing .csfasta/.qual files (with converted read names)
#   output_dir - directory to write TopHat output into
# Side effects: runs Tophat2.
sub run_tophat2 {
  my($reads_dir, $output_dir) = @_;

#  my $reads_base = "$reads_dir/$library/" . (join '_', ($experiment, $lane, $library));

  system("mkdir -p $output_dir");

system(<<END);
$bin_dir/tophat2 --transcriptome-index $transcriptome_index --no-novel-juncs --library-type fr-secondstrand --mate-inner-dist 200 --mate-std-dev 200 --color --integer-quals --bowtie1 --no-sort-bam --num-threads 7 --output-dir $output_dir $bwt_index $reads_dir/F3.csfasta $reads_dir/F5-RNA.csfasta $reads_dir/F3.qual $reads_dir/F5-RNA.qual
END
;

}

