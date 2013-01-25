#!/usr/bin/perl -w
# Converts a GFF file from modENCODE to the GFF3 format
# which Cufflinks expects (along with other Tuxedo-suite tools.)

# my $infile = "/home/jburdick/data/modencode/mRNA/Aggregate_1003.integrated_transcripts..gff3.gz";
# my $outfile = "/home/jburdick/data/modencode/mRNA/integratedTranscriptsMinimalWS220.gff3";

my $infile =
  "/home/jburdick/data/modencode/mRNA/integratedTranscriptsWS170.gff3.gz";
my $outfile = "/home/jburdick/data/modencode/mRNA/integratedTranscriptsMinimalWS190.gff3";

lift_over($infile);
get_minimal_gff("tmp.ce6.gff", $outfile);
unlink("tmp.ce4.gff");
unlink("tmp.ce6.gff");

system(
  "cat $outfile " .
  "~/gcb/data/expression/lincRNA/ancRNAs_cuff_ce6.gff " .
  "~/gcb/data/expression/lincRNA/lincRNAs_cuff_ce6.gff " .
  "| grep -v \"^\#\" " .
  " > ~/data/expression/WS190_genes.gff");

# Parses a line of attributes.
sub parse_attr {
  my ($s) = @_;

#  my @a = split /([^=;]+)=([^=;]+);/, "$s;";
  my @a = split /[=;]/, $s;   # XXX this is a bit fast-and-loose

  return @a;
}

# Simplifies a name (or set of names.)
sub simplify_IDs {
  my ($s) = @_;
  my @ids = split /,/, $s;

  return map { s/Transcript:Aggregate_1003_//; $_; } @ids;
}

# Runs "liftOver" to convert from one build to another.
# Includes some chromosome renaming.
sub lift_over {
  my($infile) = @_;

  open IN, "gunzip -c $infile |" || die;
  open OUT, ">tmp.ce4.gff" || die;
  while (<IN>) {
    next if /^\#/;
    next if /^$/;
    next if /bundle_of_reads_supporting/;

    print OUT "chr$_";
  }
  close IN;
  close OUT;

  system("liftOver -gff tmp.ce4.gff " .
    " ~/gcb/data/UCSC/ce4/ce4ToCe6.over.chain tmp.ce6.gff tmp.unMapped");
}

# Gets just the minimal features from a GFF file.
# This only includes features called "transcript",
# "exon", and "CDS".
sub get_minimal_gff {
  my($infile, $outfile) = @_;

#  open IN, "gunzip -c $infile |" || die;
  open IN, "<$infile" || die;
  open OUT, ">$outfile" || die;

  while(<IN>) {
    chomp;
    next if /^\#/;
    next if /^$/;

    my($chr, $source, $type, $start, $end, $score,
      $strand, $phase, $attr) = split /\t/;
    next if not ($type =~ /^(transcript|exon|CDS)$/);
#    $chr =~ s/^chr//;     # for now, including this
    my %a = parse_attr($attr);

    if ($type eq "transcript") {
      my @ids = simplify_IDs($a{"ID"});
  #print ($a{"ID"});
  #print "IDS = ";
  #print (join "...", @ids);
  #print "\n";
      my $id = $ids[0];

      my $gene_name = $a{"overlapping_wormbase_transcript"};
      if (not defined $gene_name) {
        $gene_name = $id;
      }
      $attr = "gene_id=$gene_name ; transcript_id=$id;";

      print OUT (join "\t", ($chr, $source, $type, $start, $end,
        $score, $strand, $phase, $attr));
      print OUT "\n";
    }

    if ($type eq "exon" || $type eq "CDS") {
      my @parents = simplify_IDs($a{"Parent"});

      # in this case, write out one line per parent
      foreach (@parents) {
        print OUT (join "\t", ($chr, $source, $type, $start, $end,
          $score, $strand, $phase, "Parent=$_"));
        print OUT "\n";
      }
    }
  }
}

close IN;
close OUT;

