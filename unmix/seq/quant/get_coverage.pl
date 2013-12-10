#!/usr/bin/perl -w
# Computes coverage for all BAM files.

# gene bounds
my $gene_bounds = "../../../data/seq/merged_genes_WS220.bed";

# just the exons of those
my $gene_bounds_split = "../../../data/seq/merged_genes_split_WS220.bed";

if (1) {
compute_coverage($gene_bounds,
  "/murrlab/seq/tophat2/WS220/20110922/",
  "rawCoverage/WS220/20110922/", "-s");
compute_coverage($gene_bounds,
  "/murrlab/seq/tophat2/WS220/20110922/",
  "rawCoverage/WS220/20110922_as/", "-S");
}


if (1) {
compute_coverage($gene_bounds,
  "/murrlab/seq/tophat2/WS220/Murray050912/",
  "rawCoverage/WS220/Murray050912/", "-s");
compute_coverage($gene_bounds,
  "/murrlab/seq/tophat2/WS220/Murray050912/",
  "rawCoverage/WS220/Murray050912_as/", "-S");
}

if (1) {
compute_coverage($gene_bounds,
  "/murrlab/seq/tophat2/WS220/Murray_52831_092812/",
  "rawCoverage/WS220/Murray_52831_092812/", "-s");
compute_coverage($gene_bounds,
  "/murrlab/seq/tophat2/WS220/Murray_52831_092812/",
  "rawCoverage/WS220/Murray_52831_092812_as/", "-S");
}

# timeseries; this shouldn't be stranded
if (undef) {
  compute_coverage("geneBounds_WS220.tsv",
    "/media/disk2/jburdick/embryo_timeseries",
    "rawCoverage_WS220/embryo_ts", "");
}

sub compute_coverage {
  my($bed_file, $bam_dir, $out_dir, $str) = @_;

  system("mkdir -p $out_dir");

  my $read_counts_file = $out_dir;
  $read_counts_file =~ s/\/+$//;
  $read_counts_file = $read_counts_file . ".readCounts.tsv";

  open COUNTS, ">$read_counts_file" || die;
  print COUNTS "sample\ttotal.reads\n";

  foreach (<$bam_dir/*.bam>) {
    my $file = $_;
    die if not $file =~ /\/?([^\/]+)\.bam/;
    my $sample_name = $1;
    my $output_file = "$out_dir/$sample_name.tsv.gz";

print "$file\n";
    # XXX temporary hack
#    next if not (($file =~ "\/03_") || ($file =~ "\/04_"));

    print "[writing coverage for $file to $output_file]\n";
    #  system ("coverageBed -abam $file -b $bed_file > $output_file");
    # for now, only including uniquely-mapping reads
#    system ("samtools view -q 1 -b $file | bedtools coverage $str -split -abam - -b $bed_file | gzip -c > $output_file");

    # only include reads which intersect at least one exon
    # XXX this is finicky
    system ("samtools view -q 1 -b $file | bedtools intersect $str -u -abam - -b $gene_bounds_split | bedtools coverage $str -split -abam - -b $gene_bounds |cut -f1-6,13-16 | gzip -c > $output_file");

    my $read_count = `samtools view -q 1 -c $file`;
    print COUNTS "$sample_name\t$read_count";
  }

  close COUNTS;
}

