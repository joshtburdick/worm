#!/usr/bin/perl -w
# Runs FIRE on all clusters.
# Note that this assumes gene names are "clone IDs" like F21D5.9,
# rather than gene names like pha-4.

# XXX shouldn't be hardwired
my $output_dir = "/home/jburdick/gcb/git/sort_paper/tf/motif/FIRE/output/";

my $FIREDIR = "/murrlab/software/FIRE-1.1a";

# Runs FIRE on a cluster.
sub run_fire {
  my ($one_cluster_file) = @_;

  # run FIRE
  system("export FIREDIR=$FIREDIR; cd /murrlab/software/FIRE-1.1a/; nice perl fire.pl " .
    "--expfile=$one_cluster_file --species=worm --exptype=discrete");
}

# Writes out one clustering.
sub write_clustering {
  my($cluster_file, $cluster_to_write, $output_file) = @_;

  open IN, "<$cluster_file" || die;
  open OUT, ">$output_file" || die;

  # copy header line
  $_ = <IN>;
  print OUT $_;

  # write out whether each gene is in the cluster
  while (<IN>) {
    chomp;
    my ($gene, $c) = split /\t/;
    my $a = ($c == $cluster_to_write ? "1" : "0");
    print OUT "$gene\t$a\n";
  }

  close IN;
  close OUT;
}

# Runs FIRE on all clusters in a clustering.
sub run_fire_all {
  my($cluster_file, $output_dir) = @_;

  system("mkdir -p $output_dir");

  # get list of clusters to run this on
  my %clusters = ();
  open IN, "<$cluster_file" || die;
  $_ = <IN>;
  while (<IN>) {
    chomp;
    my($gene, $c) = split /\t/;
    $clusters{$c} = 1;
  }
  close IN;

  foreach my $c (sort (keys %clusters)) {
    print "[running on cluster $c]\n";
    write_clustering($cluster_file, $c, "$output_dir/$c");
    run_fire("$output_dir/$c");
    unlink("$output_dir/$c");
  }

}

run_fire_all("spencerEmbryonic.tsv", "$output_dir/spencerEmbryonic");
run_fire_all("facs.tsv", "$output_dir/facs");


