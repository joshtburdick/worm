#!/usr/bin/perl -w
# Counts the motifs near genes.

use strict;

count_motifs("/murrlab/seq/igv/motif/meme/",
  "regions/WS220_5000_bp_upstream.bed",
  "./motifs_5kbUpstream.tsv", 3);
# count_motifs("/murrlab/seq/igv/histone.chip/Early-Embryos/", "histone_EE_1kbUp.tsv", 5);
# count_motifs("/murrlab/seq/igv/histone.chip/Late-Embryos/", "histone_LE_1kbUp.tsv", 5);
#count_motifs("/murrlab/seq/igv/chip.TF/Early-Embryos/", "chip.TF_EE_1kbUp.tsv", 5);
#count_motifs("/murrlab/seq/igv/chip.TF/Late-Embryos/", "chip.TF_LE_1kbUp.tsv", 5);
#count_motifs("/murrlab/seq/igv/chip.TF/Fed-L1-stage-larvae/", "chip.TF_FedL1.tsv", 5);

# count_motifs("/murrlab/seq/igv/chip.TF/Fed-L1-stage-larvae/", "chip.TF_FedL1.tsv", 5);

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
  foreach my $bw_file (<$bw_dir/*>) {
    die if not ($bw_file =~ /\/?([^\/]+)\.bw$/);
    my $base_name = $1;
    print "$base_name ";
    push @names, $base_name;

    system "bigWigAverageOverBed $bw_file $bed_file $tmp_dir/$base_name.tsv";

    # read in just one column
    open IN, "cut -f$column $tmp_dir/$base_name.tsv |" || die;
    foreach my $i (0..(@r-1)) {
      $_ = <IN>;
      chomp;
      push @{$r[$i]}, $_;
    }
    close IN;
    system "unlink $tmp_dir/$base_name.tsv";
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

