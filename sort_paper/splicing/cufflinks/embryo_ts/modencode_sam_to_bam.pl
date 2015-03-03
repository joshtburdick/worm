#!/usr/bin/perl -w
# Converts some .sam.gz files to BAM.

chdir("/media/jburdick/disk2/data/ftp/");

foreach my $sam (<sam_gz_1/*.sam.gz>) {
  print "$sam\n";

  my $bam_base = $sam;
  $bam_base =~ s/^sam_gz_1/bam/;
  $bam_base =~ s/\.sam\.gz$//;
  my $bam = "$bam_base.bam";

  system "gunzip -c $sam | samtools view -bS - " .
    "| samtools sort -m 200000000 - $bam_base";
  system "samtools index $bam";
}

