#!/usr/bin/perl -w
# Converts some .sam.gz files to BAM.

chdir("/media/jburdick/disk2/data/ftp/");

# Converts a SAM file to a BAM file, with these alterations:
# - any D in the CIGAR string is converted to an N (as
#   these seem to correspond to introns, and that's what
#   Cufflinks, at least, expects)
# - any read with a Y1 tag gets a corresponding XS tag
#   added, indicating the strand (based on alignment)
# Args:
#   sam_gz_in  - the input file
#   bam_out - the output file
sub modencode_sam_gz_to_bam {
  my($sam_gz_in, $bam_out) = @_;

  die if not $bam_out =~ /(^.+)\.bam$/;
  my $bam_base = $1;

  open IN, "gunzip -c $sam_gz_in |" || die;
  open OUT, "| samtools view -bS - | samtools sort -m 1000000000 - $bam_base" || die;

  while (<IN>) {

    # pass through header unchanged
    if (/^@/) {
      print OUT;
      next;
    }

    chomp;
    my @a = split /\t/;

    # change "D" to "N" in the CIGAR string
    $a[5] =~ s/D/N/g;

    print OUT (join "\t", @a);

    # possibly add an XS tag
    if (defined $a[12] && 
        $a[12] =~ /Y1:Z:.*_\d+_\d+_([+-])_wb220/) {
      print OUT "\tXS:A:$1\n";
    }
    else {
      print OUT "\n";
    }
  }

  close IN;
  close OUT;

  system "samtools index $bam_out";
}

foreach my $sam (<sam_gz_1/*.sam.gz>) {
  print "$sam\n";

  my $bam = $sam;
  $bam =~ s/^sam_gz_1/bam/;
  $bam =~ s/\.sam\.gz$/\.bam/;

  modencode_sam_gz_to_bam($sam, $bam);

  system "samtools index $bam";
}

