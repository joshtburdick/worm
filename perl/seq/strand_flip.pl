#!/usr/bin/perl -w
# Given SOLiD paired-end reads, flips strand of one read.

use strict;

# this should be a glob which expands to all the Tophat
# output directories to be included.
my $dir_paths = $ARGV[0];

# where to write output
my $out_dir = $ARGV[1];

# Which strand should be flipped
# - For most samples, this should be "first".
# - For unamplified samples (e.g. heat-shock samples),
#   this should be "last".
# my $flip_direction = $ARGV[2];
# XXX hacking this based on sample name

# for samtools conversion -> BAM
# my $genome_fa = "~/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/WholeGenomeFasta/genome.fa";
my $genome_fa = "/var/tmp/data/tophat2/WS220/genome.fa";

sub flip_file {
  my($in_file, $out_file) = @_;

  # XXX hack
  my $flip_direction = "first";
  if ($in_file =~ /(ges1|lit1|pop1)_/) {
    $flip_direction = "last";
  }

  open IN, "samtools view $in_file |" || die;
  open OUT, "|samtools view -bT $genome_fa - > $out_file";

  while (<IN>) {
    chomp;
    my @a = split /\t/;

    my $flags = ($a[1] / 0x10) & 0x0f;

    # if this is the "last segment in the template"...
    if (($flip_direction eq "last") && ($a[1] & 0x80)) {
      # flip its strand
      $a[1] ^= 0x30;
    }

    # ... or if it's the first
    if (($flip_direction eq "first") && (!($a[1] & 0x80))) {
      # flip its strand
      $a[1] ^= 0x30;
    }

    print OUT (join "\t", @a);
    print OUT "\n";
  }

  close IN;
  close OUT;
}

# XXX
#flip_file("/murrlab/seq/tophat2/Murray050912/ce6_genes/Murray050912_Pool1Pool2_01_F21D5.9/accepted_hits.bam",
#  "/murrlab/seq/tophat2/Murray050912/strand_flip/01_F21D5.9.bam");

#  my $dir_paths = "/murrlab/seq/tophat2/Murray050912/ce6_genes/Murray050912_Pool1Pool2_0"
#    . $lane . "_*";
#  my $dir_paths = "/murrlab/seq/tophat2/Murray_52831_092812/lane$lane/Murray_52831_092812_*";

foreach my $dir (glob $dir_paths) {
  chomp $dir;
#  print "$dir\n";

  die if not ($dir =~ /_(0[1-6](?:\.edit)?_.*)$/);
#  die if not ($dir =~ /\/([^\/]+)$/);
  my $id = $1;

#    my $out_file = "/murrlab/seq/tophat2/Murray050912/strand_flip/$id.bam";
#  my $out_file = "/murrlab/seq/tophat2/Murray_52831_092812/strand_flip/$id.bam";
  my $out_file = "$out_dir/$id.bam";

  # XXX another hack: removing ".edit" from the names of
  # output files
  $out_file =~ s/\.edit//;

  # skip cases which are already done
  next if -e $out_file;

# XXX was "accepted_hits.bam"
print("$dir/merged.bam $id $out_file\n");
    flip_file("$dir/merged.bam", "tmp.bam");
    my $out_prefix = $out_file;
    $out_prefix =~ s/\.bam$//;
    system("samtools sort -m 1000000000 tmp.bam $out_prefix");
    system("samtools index $out_file");
    unlink("tmp.bam");
}

