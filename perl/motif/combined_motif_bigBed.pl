#!/usr/bin/perl -w
# Creates a bigBed file of all the motifs.
# ??? is a 9-column .bed file legit? if not, I can
# tack on the remaining columns; skipping for now.
# This writes out:
# - a merged .bam file (which is fast, and smallish), and
# - a .bigBed file (which is larger, but colorized)

# FIXME: this should rename the "MA____"-style motif names
# from JASPAR into something more readable.

use strict;

use Convert::Color;
use Convert::Color::HSV;
use Convert::Color::RGB8;

my $motif_dir = "/murrlab/seq/igv/motif/known/";
my $output_base = "/murrlab/seq/igv/motif/known";

my $genome_sizes = "/var/tmp/data/tophat2/WS220/genome_with_MtDNA.sizes";
my $tmp_base = "/media/disk2/jburdick/motif_tmp";

# Gets the list of unique motif names.
sub get_motif_names {
  my $collapsed_motif_list = "../../tf/motif/motifFilter.tsv";
  my %h = ();

  open IN, "<$collapsed_motif_list" || die;
  $_ = <IN>;    # skip first line

  # loop through the list of motifs, keeping canonical name
  while (<IN>) {
    chomp;
    my @a = split /\t/;

    # XXX for now, just skipping this motif; it's not clear
    # why it doesn't have a .bam file
    next if $a[2] eq "MA0199.1";

    # skip de novo motifs
    if (!($a[2] =~ /hier\./)) {
      $h{ $a[2] } = 1;
    }
  }
  close IN;

  # return the keys, sorted ignoring case
  return sort { lc($a) cmp lc($b) } (keys %h);
}

# get motif names
my @motif_names = get_motif_names();

print "[merging motifs...]\n";
my @motif_files = map { "$_.bam"; } @motif_names;

system("cd $motif_dir; samtools merge $output_base.bam " . 
  (join " ", @motif_files));
system("samtools index $output_base.bam");

exit(0);

# definition of how to color these
my %motif_color = ();
foreach my $i (0..(@motif_names-1)) {
  foreach my $score (30..60) {

    my $c = Convert::Color::HSV->new(
      $i * (360/@motif_names), 1 - ($score / 30 - 1), 0.9);
    my $c1 = $c->convert_to("rgb8");
    my $rgb = join ",", $c1->rgb8();

    $motif_color{ $motif_names[$i] . "_" . $score } = $rgb;
  }
}
print join " ", %motif_color;

print("\n[colorizing]\n");
open IN, "bedtools bamtobed -i $output_base.bam |" || die;
open OUT, ">$tmp_base.1.bed" || die;
while (<IN>) {
  chomp;
  my @a = split /\t/;
  my $score = $a[4];
  if ($score > 60) {
    $score = 60;
  }
  my $c = $motif_color{ $a[3] . "_" . $score };
  die $a[3] . "_" . $score if not defined $c;

  print OUT join "\t", (@a, $a[1], $a[2], $c);
  print OUT "\n";
}
close IN;
close OUT;

print "[converting to bigBed]\n";
system("bedToBigBed $tmp_base.1.bed $genome_sizes $output_base.bigBed");
unlink("$tmp_base.1.bed");

