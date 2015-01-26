#!/usr/bin/perl -w
# Counts the motifs near genes.
# (Also counts averages of other things.)

use strict;

if (1) {
  foreach my $i (1..5) {
    count_motifs_chip("/media/disk2/jburdick/chip_bw/",
      "regions/upstream_" . $i . "kb_WS220.bed",
      "TF_chip_" . $i . "kb");
    }
}

if (undef) {
write_motif_counts("/murrlab/seq/igv/motif/meme/",
  "upstream_liftOver_WS220.bed",
  "./knownMotif_5kbUp.tsv", 3);
write_motif_counts("/home/jburdick/tmp/meme_denovo_bw/",
  "upstream_liftOver_WS220.bed",
  "./deNovoMotif_5kbUp.tsv", 3);
}

if (undef) {
  write_motif_counts("/home/jburdick/tmp/meme_shuffled_bw/",
    "upstream_liftOver_WS220.bed",
    "./shuffledKnownMotif_5kbUp.tsv", 3);
  write_motif_counts("/home/jburdick/tmp/meme_shuffled_bw/",
    "upstream_liftOver_WS220_0.5cons.bed",
    "./shuffledKnownMotif_5kbUp_0.5cons.tsv", 3);
}


# write_motif_counts("/home/jburdick/tmp/meme_denovo_bw/",
#   "upstreamRegionsWS220_5kb_nogenes.bed",
#   "./foo.tsv", 3);

# write_motif_counts("/murrlab/seq/igv/motif/meme/",
#  "upstreamRegionsWS220_5kb_nogenes.bed",
#   "./motifs_5kb_nogenes.tsv", 3);

# write_motif_counts("/murrlab/seq/igv/motif/meme/",
#  "upstreamRegionsWS220_5kb_nogenes_cons_0.5.bed",
#   "/home/jburdick/tmp/motifs_5kb_nogenes_cons_0.5.tsv", 3);

# count_motifs("/murrlab/seq/igv/histone.chip/Early-Embryos/", "histone_EE_1kbUp.tsv", 5);
# count_motifs("/murrlab/seq/igv/histone.chip/Early-Embryos/", "histone_EE_1kbUp.tsv", 5);

# count_motifs_chip("/murrlab/seq/igv/chip.TF/", "TF");
# count_motifs_chip("/murrlab/seq/igv/histone.chip/", "histone.chip");

# count_motifs("/murrlab/seq/igv/histone.chip/Early-Embryos/", "histone_EE_1kbUp.tsv", 5);
# count_motifs("/murrlab/seq/igv/histone.chip/Late-Embryos/", "histone_LE_1kbUp.tsv", 5);
#count_motifs("/murrlab/seq/igv/chip.TF/Early-Embryos/", "chip.TF_EE_1kbUp.tsv", 5);
#count_motifs("/murrlab/seq/igv/chip.TF/Late-Embryos/", "chip.TF_LE_1kbUp.tsv", 5);
#count_motifs("/murrlab/seq/igv/chip.TF/Fed-L1-stage-larvae/", "chip.TF_FedL1.tsv", 5);

# count_motifs("/murrlab/seq/igv/chip.TF/Fed-L1-stage-larvae/", "chip.TF_FedL1.tsv", 5);

sub count_motifs_chip {
  my($base_dir, $regions_file, $output_dir) = @_;
  my @stages = `ls $base_dir/`;

  system("mkdir -p $output_dir");
  foreach my $stage (@stages) {
    chomp $stage;
    print "stage $stage\n";
#    count_motifs("$base_dir/$stage/",
#      "upstream_liftOver_WS220.bed",
#      "regions/WS220_5000_bp_upstream.bed",
#      "TF_chip/" . $name . "_" . $stage . "_5kbUp.tsv", 5);
    write_motif_counts("$base_dir/$stage/",
      $regions_file,
      $output_dir . "/" . $stage . ".tsv", 5);
  }
}

sub count_motifs {
  my($bw_dir, $regions_bed, $output_file, $column) = @_;

#  my $bed_file = "upstreamRegionsWS220_10kb.bed";
  my $bed_file = $regions_bed;
  my $tmp_dir = "./";

  # read in names from BED file
  my @r = ();
  open IN, "cut -f4 $bed_file |" || die;
  while (<IN>) {
    chomp;
    push @r, [$_];
  }
  close IN;

  # read in data
  my @names = ();

  open LOG, ">log.txt" || die;

  foreach my $bw_file (<$bw_dir/*.bw>) {
    die if not ($bw_file =~ /\/?([^\/]+)\.bw$/);
    my $base_name = $1;
    print LOG "$bw_file $base_name\n";
    print "$bw_file    $base_name\n";
    push @names, $base_name;

    system("touch tmp.tsv");
    system "bigWigAverageOverBed $bw_file $bed_file tmp.tsv";

    # read in just one column
    open IN, "cut -f$column tmp.tsv |" || die;
    foreach my $i (0..(@r-1)) {
      $_ = <IN>;
#      die if not defined $_;
      # XXX now skipping cases in which something fails
      if (not defined $_) {
        $_ = "";
      }
      chomp;
      push @{$r[$i]}, $_;
    }
    close IN;
    system "unlink tmp.tsv";
  }

  # write output
  open OUT, ">$output_file" || die;
  print OUT "gene\t";
  print OUT (join "\t", @names);
  print OUT "\n";
  foreach (@r) {
    print OUT (join "\t", @$_);
    print OUT "\n";
  }
  close OUT;
}

# Counts motifs, and writes out results as it goes
# (with one line per motif, "horizontally", for
# lower memory usage.)
sub write_motif_counts {
  my($bw_dir, $regions_bed, $output_file, $column) = @_;

  open OUT, ">$output_file" || die;

  # write out names from BED file
  my @genes = ();
  open IN, "cut -f4 $regions_bed |" || die;
  while (<IN>) {
    chomp;
    push @genes, $_;
    print OUT "\t$_";
  }
  print OUT "\n";
  close IN;

  # write out each motif
  foreach my $bw_file (<$bw_dir/*.bw>) {
    die if not ($bw_file =~ /\/?([^\/]+)\.bw$/);
    my $base_name = $1;

    # count motifs
    system("touch tmp.tsv");
    system "bigWigAverageOverBed $bw_file $regions_bed tmp.tsv";

    # read in just one column
    # XXX note that this assumes the average ou
    print OUT $base_name;
    open IN, "cut -f$column tmp.tsv |" || die;
    while (<IN>) {
      die if not defined $_;
      chomp;
      print OUT "\t$_";
    }
    print OUT "\n";
    close IN;
    system "unlink tmp.tsv";
  }

  close OUT;
}

