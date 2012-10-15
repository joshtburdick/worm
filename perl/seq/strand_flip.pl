#!/usr/bin/perl -w
# Given SoLID paired-end reads, flips strand of second read.

use strict;

# for samtools conversion -> BAM
my $genome_fa = "~/data/seq/Caenorhabditis_elegans/UCSC/ce6/Sequence/WholeGenomeFasta/genome.fa";

sub flip_file {
  my($in_file, $out_file) = @_;

  open IN, "samtools view $in_file |" || die;
  open OUT, "|samtools view -bT $genome_fa - > $out_file";

  while (<IN>) {
    chomp;
    my @a = split /\t/;

    my $flags = ($a[1] / 0x10) & 0x0f;

    # if this is the "last segment in the template"...
    if ($a[1] & 0x80) {
      # flip its strand
      $a[1] ^= 0x30;
    }

    print OUT (join "\t", @a);
    print OUT "\n";
  }

  close IN;
  close OUT;
}

# XXX
#flip_file("/murrlab/seq/tophat2/Murray050912/ce6_genes/Murray050912_Pool1Pool2_01_F21D5.9/accepted_hits.bam",
#  "/murrlab/seq/tophat2/Murray050912/strand_flip/01_F21D5.9.bam");

flip_file("/murrlab/seq/tophat2/Murray_52831_092812/lane1/Murray_52831_092812_01_ceh26n/accepted_hits.bam", "~/ceh26n.bam");

if (undef) {
foreach my $lane (2..6) {
  my $dir_paths = "/murrlab/seq/tophat2/Murray050912/ce6_genes/Murray050912_Pool1Pool2_0"
    . $lane . "_*";
  foreach my $dir (glob $dir_paths) {
    chomp $dir;
    print "$dir\n";

    die if not ($dir =~ /_(0[1-6]_.*)$/);
    my $id = $1;

    my $out_file = "/murrlab/seq/tophat2/Murray050912/strand_flip/$id.bam";
    flip_file("$dir/accepted_hits.bam", $out_file);
    system("samtools index $out_file");
  }
}
}

