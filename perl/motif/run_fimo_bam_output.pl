#!/usr/bin/perl -w
# Runs FIMO for some set of motifs, writing output in .bam format.

use strict;

use POSIX qw/log10 floor/;

use Bio::Seq;

my $meme_bin = "~jburdick/meme/bin/";
my $fasta = "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa";

# sizes of chromosomes, for SAM->BAM conversion
my $fasta_sizes = $fasta . ".sizes";

# my $meme_matrix_path = "/home/jburdick/gcb/data/tf/meme/motif_databases/";
my $meme_matrix_path = "";  # XXX temporary hack
# my $meme_matrix_path = "/home/jburdick/gcb/git/tf/motif/meme_file/";
# my $meme_matrix_path = "/home/jburdick/gcb/git/tf/motif/";
# my $meme_matrix_path = "/home/jburdick/gcb/git/tf/motif/shuffle/meme_files/";

my @meme_files = (
  qw/jolma2013 JASPAR_CORE_2009 dmmpmm2009 dpinteract/,
  qw/fly_factor_survey flyreg.v2 hallikas2006 homeodomain/,
  qw/idmmpmm2009 JASPAR_FAM_2008 JASPAR_PHYLOFACTS_2008/,
  qw/JASPAR_POLII_2008 jolma2010 macisaac_theme.v1 macisaac_yeast.v1/,
  qw/prodoric regtransbase scpd_matrix uniprobe_mouse uniprobe_worm/,
  qw/wei2010_human_mws wei2010_mouse_mws wei2010_mouse_pbm zhao2011/);

# just including a subset of these
# @meme_files = (
#   qw/jolma2013_shuffled JASPAR_CORE_2009_insects_shuffled 
# JASPAR_CORE_2009_nematodes_shuffled JASPAR_CORE_2009_vertebrates_shuffled/);
@meme_files = qw/jolma2013 JASPAR_CORE_2009_insects JASPAR_CORE_2009_nematodes JASPAR_CORE_2009_vertebrates/;

# ??? just trying to understand why MA0199.1 wasn't there
# @meme_files = qw/JASPAR_CORE_2009_insects/;

if (undef) {
foreach my $f (@meme_files) {
  print "[running on MEME database file $f.meme]\n";
  fimo_to_bedGraph("$f.meme", "/media/disk2/jburdick/meme/tmp2/");
}
}

# fimo_to_bedGraph("hughes_motif.meme",
#   "/home/jburdick/tmp/fimo_hughes/");

# fimo_to_bedGraph("/home/jburdick/gcb/data/tf/bhambhani2014/bhambhani2014.meme",
#   "/home/jburdick/tmp/fimo_bhambhani");
# fimo_to_bedGraph("/murrlab/jburdick/src/tf/meme/TCF_LEF.meme",
#   "/home/jburdick/tmp/fimo_TCF_LEF");

# fimo_to_bedGraph(
#   "/home/jburdick/gcb/git/data/tf/hughes_motif_20141202.meme",
#   "/home/jburdick/tmp/fimo_hughes_20141202.meme");

# fimo_to_bedGraph(
#   "/home/jburdick/gcb/git/tf/motif/shuffle/meme_files/jolma2013_shuffled.meme",
#   "/home/jburdick/tmp/fimo/jolma2013_shuffled/");
fimo_to_bedGraph(
  "/home/jburdick/gcb/git/tf/motif/shuffle/meme_files/hughes_motif_20141202_shuffled.meme",
  "/home/jburdick/tmp/fimo/hughes_motif_20141202_shuffled/");

# Utility to convert from .sam to .bam .
# Also erases the .sam file, and indexes the .bam file.
sub sam_to_bam {
  my($sam_file) = @_;
  my $base_name = $sam_file;
  $base_name =~ s/.sam$//;
  my $bam_file = "$base_name.bam";

  system("samtools view -b -t $fasta_sizes $sam_file > $bam_file");
  unlink($sam_file);

  system("samtools sort $bam_file $base_name");
  system("samtools index $bam_file");
}

# Runs FIMO, converting its output to bedGraph files (which nonetheless
# end in .wig, which seems to be what igvtools expects.)
sub fimo_to_bedGraph {
  my($meme_file, $output_path) = @_;

print("$meme_matrix_path/$meme_file\n");
  # if the database isn't there, skip it
  return if not -e "$meme_matrix_path/$meme_file";

  print STDERR "[running FIMO on $meme_file]\n";
  system("mkdir -p $output_path");

# cat("$meme_bin/fimo --text $meme_matrix_path/$meme_file $fasta\n");
  open IN, "$meme_bin/fimo --thresh 1e-3 --text $meme_matrix_path/$meme_file $fasta |" || die;

  $_ = <IN>;
  die if not /#pattern.*name.*sequence.*name.*start.*stop.*strand.*score.*p-value.*q-value.*matched sequence/;

  # current motif we're writing out
  my $current_motif = undef;

  while (<IN>) {
    chomp;
    next if /^#/;        # skip header
    my($motif_name, $chr, $start, $end, $strand, $log_odds,
      $p_value, $q_value, $site) = split /\t/;

    if (!(defined $current_motif) || !($motif_name eq $current_motif)) {
      if (defined $current_motif) { 
        close OUT;

        # compress current file
        sam_to_bam("$output_path/$current_motif.sam");
      }
      $current_motif = $motif_name;
      print STDERR "[motif $current_motif]\n";
      open OUT, ">$output_path/$current_motif.sam";
#      open OUT, "| samtools view -b -t $fasta_sizes > $output_path/$current_motif.bam" || die;
    }

    # height of the wiggle here reflects how good the match is
    # XXX I don't know quite what this p-value means
    my $mapq = floor( -10.0 * log10($p_value) );

#    print OUT join "\t", ($chr, $start, $start+1, "$mapq\n");

    # if this is on "-" strand, write that, and revcomp the sequence
    my $flag = ($strand eq "+" ? 0 : 16);
    if ($strand eq "-") {
      my $s = Bio::Seq->new(-seq=>$site);
      my $s1 = $s->revcom();
      $site = $s1->seq();
    }

    my $length = 1 + $end - $start;

    # ??? why was I doing this?
    # subtracting 2 from start, so this lines up with genomic sequence
    # (was $start-2)
    print OUT join "\t", ($current_motif, $flag, $chr, $start,
      $mapq, $length . "M", "*", 0, 0, $site, "*\n");
  }

  close IN;
  close OUT;

  # compress the last file
  sam_to_bam("$output_path/$current_motif.sam");
}

