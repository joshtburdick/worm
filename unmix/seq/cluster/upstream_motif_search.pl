#!/usr/bin/perl -w
# Looks for motifs upstream of genes in clusters.
# (Actually, currently just writes out .fa files of
# sequence upstream of genes.)

use Bio::DB::Fasta;

use List::MoreUtils qw/ uniq /;

# where the worm genome .fa files are
my $fasta = Bio::DB::Fasta->new(
  "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa");

# the upstream region for each gene
my $upstream_region_bed_file =
#  "/home/jburdick/gcb/git/tf/motif/motifCount/upstreamRegionsWS220_1kb.bed";
  "/home/jburdick/gcb/git/tf/motif/motifCount/upstreamRegionsWS220_5kb_nogenes.bed";

foreach my $a (qw/hier.50clusters hier.100clusters hier.200clusters hier.ts.50clusters hier.ts.100clusters hier.ts.200clusters/) {
  print "[writing .fa for clustering $a]\n";
  write_cluster_gene_fa("hierarchical/$a", "hierarchical/$a/5kb.intergenic/fa");
}
foreach my $a (qw/wnet_pow11_mch0.1 wnet_pow11_mch0.2
    wnet.ts_pow11_mch0.1 wnet.ts_pow11_mch0.2/) {
  print "[writing .fa for clustering $a]\n";
  write_cluster_gene_fa("WGCNA/$a", "WGCNA/$a/5kb.intergenic/fa");
}

# Reads in the clustering, as a hash.
sub get_clustering {
  my($cluster_file) = @_;
  my %h = ();

  open IN, "<$cluster_file" || die;
  $_ = <IN>;
  die if not /gene.*cluster/; 
  while (<IN>) {
    chomp;
    s/\"//g;
    my($skipped, $gene, $cluster) = split /\t/;
# print "$gene $cluster\n";
    $h{$gene} = $cluster;
  }
  close IN;

  return %h;
}

# Creates a file of upstream regions.
# FIXME add a switch indicating whether to get all / conserved sequence
sub write_cluster_upstream_regions {
  my($clustering_ref, $cluster, $fa_file) = @_;
  my %clustering = %$clustering_ref;

  open IN, "<$upstream_region_bed_file" || die;
#  open OUT, "| /home/jburdick/gcb/git/tf/motif/get_upstream_conserved.pl > $fa_file" || die;
  open OUT, "| /home/jburdick/gcb/git/tf/motif/get_dna_in_region.pl > $fa_file" || die;

  while (<IN>) {
    chomp;
    my($chr, $a, $b, $name, $score, $strand) = split /\t/;

    # only include genes in this cluster
    my $c1 = $clustering{$name};
    next if not (defined $c1);
    next if not ($c1 eq $cluster);

    $chr =~ s/chr//;

    # force region to be >= 500 bp
    # (although conservation may mean less sequence is printed)
    my $length = $b - $a;
    if ($length < 500) {
      if ($strand eq "+") {
        $a = $b - 500;
        if ($a < 1) {
          $a = 1;
        }
      }
      if ($strand eq "-") {
        $b = $a + 500;
      }
    }

    # print this region to the filter (which will write the DNA)
    print OUT join "\t", ($chr, $a, $b, $name, $score, $strand);
    print OUT "\n";
  }

  close IN;
  close OUT;
}

# Writes out file of upstream sequence for one cluster.
# Args:
#   dir - directory which contains a tab-separated file "clusters.tsv",
#     which has three columns:
#       row name - (not used)
#       gene - name of a gene
#       cluster - number indicating which cluster that gene is in
# Side effects: creates a directory in that directory, called "meme",
#   containing result of running meme on genes in each cluster.
sub write_cluster_gene_fa {
  my($dir, $output_dir) = @_;

  system("mkdir -p $output_dir");

  my %clustering = get_clustering("$dir/clusters.tsv");
  my @clusters = sort(uniq(values %clustering));

  $| = 1;   # for printing progress
  foreach my $c (@clusters) {
    print "\b" x 80;
    print "[writing .fa upstream of genes in cluster $c]";

    my $upstream_fa = "$output_dir/$c.fa";

    # write file of sequences upstream of genes in cluster
    write_cluster_upstream_regions(\%clustering, $c, $upstream_fa);
  }

  $| = 0;
  print "\n";
}

