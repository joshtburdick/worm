#!/usr/bin/perl -w
# Computes coverage for all BAM files.

# confusingly, "geneBounds_flip.tsv" is the "sense" strand quantification

if (undef) {
compute_coverage("geneBounds_flip.tsv",
  "/murrlab/seq/tophat2/20110922/",
  "20110922");
compute_coverage("geneBounds.tsv",
  "/murrlab/seq/tophat2/20110922/",
  "20110922_as");
}

if (undef) {
compute_coverage("geneBounds_flip.tsv",
  "/murrlab/seq/tophat2/Murray050912/strand_flip",
  "Murray_050912");
compute_coverage("geneBounds.tsv",
  "/murrlab/seq/tophat2/Murray050912/strand_flip",
  "Murray_050912_as");
}
if (1) {
compute_coverage("geneBounds_flip.tsv",
  "/murrlab/seq/tophat2/Murray_52831_092812/strand_flip",
  "Murray_52831_092812");
compute_coverage("geneBounds.tsv",
  "/murrlab/seq/tophat2/Murray_52831_092812/strand_flip",
  "Murray_52831_092812_as");
}

sub compute_coverage {
  my($bed_file, $bam_dir, $out_dir) = @_;

  system("mkdir -p $out_dir");

  open COUNTS, ">read_counts.tsv" || die;
  print COUNTS "sample\ttotal.reads\n";

  foreach (<$bam_dir/*.bam>) {
    my $file = $_;
    die if not $file =~ /\/?([^\/]+)\.bam/;
    my $sample_name = $1;
    my $output_file = "$out_dir/$sample_name.tsv.gz";

    # XXX temporary hack
print "$file\n";
#    next if not (($file =~ "\/03_") || ($file =~ "\/04_"));

    print "[writing coverage for $file to $output_file]\n";
    #  system ("coverageBed -abam $file -b $bed_file > $output_file");
    # for now, only including uniquely-mapping reads
    #  system ("samtools view -q 255 -b $file | coverageBed -s -split -abam stdin -b $bed_file | gzip -c > $output_file");
    system ("samtools view -q 1 -b $file | coverageBed -s -split -abam stdin -b $bed_file | gzip -c > $output_file");
    my $read_count = `samtools view -q 1 -c $file`;
    print COUNTS "$sample_name\t$read_count";
  }

  close COUNTS;
}

