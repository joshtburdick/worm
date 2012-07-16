#!/usr/bin/perl -w
# Computes coverage for all BAM files.

# confusingly, this is the "sense" strand quantification
compute_coverage("geneBounds_flip.tsv", "/murrlab/seq/tophat2/Murray050912/strand_flip", "Murray_050912");

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
