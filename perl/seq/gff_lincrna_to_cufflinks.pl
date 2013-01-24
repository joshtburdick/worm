#!/usr/bin/perl -w
# Converts GFF files from Bartel lab to what Cufflinks
# seems to expect.

# ??? include liftOver in this?

my $dir = "/home/jburdick/gcb/data/expression/lincRNA/";

convert("$dir/ancRNAs_W3PSeq3_ce10.gtf",
  "$dir/ancRNAs_cuff.gff");
convert("$dir/lincRNAs_W3PSeq3_ce10.gtf",
  "$dir/lincRNAs_cuff.gff");

# XXX put this elsewhere?
#system("cat /home/jburdick/data/modencode/mRNA/integratedTranscriptsMinimalWS220.gff3 " .
#  " $dir/ancRNAs_cuff.gff $dir/lincRNAs_cuff.gff > " .
#  " /home/jburdick/data/expression/WS220_genes.gff");

sub convert {
  my($infile, $outfile) = @_;
  system("gffread $infile -E -m chromosome_mapping.txt -o $outfile");
}

# um, should've read the gffread documentation
sub convert_deprecated {
  my($infile, $outfile) = @_;

  open IN, "<$infile" || die;
  open OUT, ">$outfile" || die;

  while (<IN>) {
    chomp;
    my($chr, $source, $type, $start, $end, $score,
      $strand, $phase, $attr) = split /\t/;
    next if not ($type =~ /^(transcript|exon|CDS)$/);

    if ($type eq "transcript") {
      die if not ($attr =~
        /gene_id \"([^\"]+)\"; transcript_id \"([^\"])\";/);
      $attr = "ID=$1;gene_name=$2";
    }
    
    if ($type eq "exon" || $type eq "CDS") {
      die if not ($attr =~
        /gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\";/);
      $attr = "Parent=$1";
    }

    print OUT (join "\t", ($chr, $source, $type, $start, $end,
      $score, $strand, $phase, $attr));
    print OUT "\n";
  }

  close IN;
  close OUT;
}




