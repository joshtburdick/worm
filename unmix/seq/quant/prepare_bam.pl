#!/usr/bin/perl -w
# Various processing of read files, for submitting to the SRA.
# - changes header
# - renames slightly
# - computes MD5 sums

use strict;

my $new_header_file = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/WholeGenomeFasta/sra_header.sam";

# Reheaders one file.
sub reheader_file {
  my($input_file, $output_file) = @_;

  print("samtools reheader $new_header_file $input_file > $output_file\n");

# presumably we don't need to sort or index it
#    "| samtools sort -m 2000000000 - $output_file_1");
#    print("samtools index $output_dir/$a.bam");
}

open IN, "<file.to.experiment.tsv" || die;
$_ = <IN>;
while (<IN>) {
  chomp;
  my($i, $a, $b) = split /\t/;

  reheader_file($a, $b);
}

system("md5sum *.bam > md5_sums.txt)");

