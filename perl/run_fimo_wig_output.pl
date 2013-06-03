#!/usr/bin/perl -w
# Runs FIMO for some set of motifs, writing output in .wig format.

use strict;

use POSIX qw/log10 floor/;

my $meme_bin = "~jburdick/meme/bin/";
# my $fasta = "/murrlab/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa";
my $fasta = "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa";

# my $meme_matrix_path = "/home/jburdick/gcb/data/tf/meme/motif_databases/";
# my $meme_matrix_path = ".";  # XXX temporary hack
my $meme_matrix_path = "/home/jburdick/gcb/git/tf/motif/meme_file/";

my @meme_files = (
  qw/jolma2013 JASPAR_CORE_2009 dmmpmm2009 dpinteract/,
  qw/fly_factor_survey flyreg.v2 hallikas2006 homeodomain/,
  qw/idmmpmm2009 JASPAR_FAM_2008 JASPAR_PHYLOFACTS_2008/,
  qw/JASPAR_POLII_2008 jolma2010 macisaac_theme.v1 macisaac_yeast.v1/,
  qw/prodoric regtransbase scpd_matrix uniprobe_mouse uniprobe_worm/,
  qw/wei2010_human_mws wei2010_mouse_mws wei2010_mouse_pbm zhao2011/);

# fimo_to_bedGraph("TCF_LEF.meme", "TCF_LEF");
# fimo_to_bedGraph("pha4_hardcoded.meme", "pha4_hardcoded");

#foreach my $f (@meme_files) {
#  print "[running on MEME database file $f.meme]\n";
#  fimo_to_bedGraph("$f.meme", "/home/jburdick/tmp/meme");
#}

fimo_to_bedGraph("meme_hier_ts_200clusters.meme", "/home/jburdick/tmp/meme");


# Runs FIMO, converting its output to bedGraph files (which nonetheless
# end in .wig, which seems to be what igvtools expects.)
sub fimo_to_bedGraph {
  my($meme_file, $output_path) = @_;

  # if the database isn't there, skip it
  return if not -e "$meme_matrix_path/$meme_file";

  print STDERR "[running FIMO on $meme_file]\n";
  system("mkdir -p $output_path");

  open IN, "$meme_bin/fimo --text $meme_matrix_path/$meme_file $fasta |" || die;

  $_ = <IN>;
  die if not /#pattern.*name.*sequence.*name.*start.*stop.*strand.*score.*p-value.*q-value.*matched sequence/;

  # current motif we're writing out
  my $current_motif = undef;

  while (<IN>) {
    chomp;
    my($motif_name, $chr, $start, $end, $strand, $log_odds,
      $p_value, $q_value, $site) = split /\t/;

    if (!(defined $current_motif) || !($motif_name eq $current_motif)) {
      if (defined $current_motif) { 
        close OUT;
      }
      $current_motif = $motif_name;
      print STDERR "[motif $current_motif]\n";
      open OUT, ">$output_path/$current_motif.wig";
#      print OUT "track type=bedGraph\n";
    }

    # height of the wiggle here reflects how good the match is
    # XXX I don't know quite what this p-value means
    my $mapq = floor( -10.0 * log10($p_value) );

    # note that we only write out the leftmost base of this,
    # as it simplifies computing some statistics
    print OUT join "\t", ($chr, $start, $start+1, "$mapq\n");
  }

  close IN;
  close OUT;

  # FIXME: gzip all output files?
}

