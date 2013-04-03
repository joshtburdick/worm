#!/usr/bin/perl -w
# Looks for motifs upstream of genes in clusters.

use Bio::DB::Fasta;

use List::MoreUtils qw/ uniq /;

# where the worm genome .fa files are
my $fasta = Bio::DB::Fasta->new(
  "/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa");

# the upstream region for each gene
my $upstream_region_bed_file =
  "/home/jburdick/gcb/git/tf/motif/motifCount/upstreamRegionsWS220_5kb_nogenes.bed";

# background (for now, the entire genome)
my $bg_markov_model = "~/gcb/git/tf/motif/Ce_WS220.order1markov.txt";

# path to MEME binaries
my $meme_path = "/home/jburdick/meme/";

find_motifs_in_cluster("hierarchical/hier.50clusters");

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
sub write_cluster_upstream_regions {
  my($clustering_ref, $cluster, $fa_file) = @_;
  my %clustering = %$clustering_ref;

  open IN, "<$upstream_region_bed_file" || die;
  open OUT, ">$fa_file" || die;

  while (<IN>) {
    chomp;
    my($chr, $a, $b, $name, $score, $strand) = split /\t/;

    # only include genes in this cluster
    my $c1 = $clustering{$name};
    next if not (defined $c1);
    next if not ($c1 eq $cluster);

    $chr =~ s/chr//;

    print OUT ">$name\n";
    if ($strand eq "+") {
      print OUT $fasta->seq($chr, $a => $b) . "\n";
    }
    elsif ($strand eq "-") {
      print OUT $fasta->seq($chr, $b => $a) . "\n";
    }
    else {
      die "unknown strand $strand\n";
    }
  }

  close IN;
  close OUT;
}

# Looks for motifs in one clustering.
# Args:
#   dir - directory which contains a tab-separated file "clusters.tsv",
#     which has three columns:
#       row name - (not used)
#       gene - name of a gene
#       cluster - number indicating which cluster that gene is in
# Side effects: creates a directory in that directory, called "meme",
#   containing result of running meme on genes in each cluster.
sub find_motifs_in_cluster {
  my($dir) = @_;
  system("mkdir -p $dir/meme");

  my %clustering = get_clustering("$dir/clusters.tsv");

  my @clusters = sort(uniq(values %clustering));

  foreach my $c (@clusters) {
    print "[running meme on cluster $c]\n";

    my $upstream_fa = "$dir/meme/$c.fa";

    my $output_dir = "$dir/meme/$c";

    # write file of sequences upstream of genes in cluster
    write_cluster_upstream_regions(\%clustering, $c, $upstream_fa);

    # run meme
    system("\\rm -rf $output_dir");
    system("$meme_path/bin/meme $upstream_fa -o $output_dir " .
      " -bfile $bg_markov_model -maxsize 100000000 -maxw 15 " .
      " -dna -revcomp -nmotifs 3");

    # unlink($upstream_fa);
  }
}

