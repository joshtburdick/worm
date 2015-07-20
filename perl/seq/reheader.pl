#!/usr/bin/perl -w
# Creates a new SAM header, and reheaders a directory of .bam files,
# for submitting them to the SRA.

use strict;

my $new_header_file = "   ";

# Creates a new header file.
sub create_header {



}

# Reheaders one directory's worth of files.
sub reheader_dir {
  my($new_header_file, $input_dir, $output_dir) = @_;

  # loop through the files
  foreach (<$input_dir/*.bam>) {
    chomp;
    /\/([^\/]+)\.bam$/;
    my $a = $1;

    print "$a\n";

    print("samtools reheader $new_header $input_dir/$a.bam > $output_dir/$a.bam");
    print("samtools index $output_dir/$a.bam");
  }
}














