#!/usr/bin/perl -w
# Aligns the modENCODE embryonic timeseries.

my $ftp_host = "data.modencode.org";
my $ftp_dir = "C.elegans/mRNA/RNA-seq/raw-seqfile_fastq/";

my $tmp_dir = "/var/tmp/embr_ts_align/";
my $output_dir = "/media/disk2/jburdick/embryo_timeseries/";

# Tophat-specific stuff
my $bwt_index = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/Bowtie2Index/genome";
my $gene_gtf = "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Annotation/Genes/genes.gtf";
# my $gene_gtf = "/home/jburdick/data/expression/WS220_genes.gff";

# my $tophat2_bin = "/home/jburdick/bin/tophat2";
my $tophat2_bin =
  "/murrlab/software/tophat-2.0.8.Linux_x86_64/tophat2";
my $bin_dir = "/murrlab/software/bin";

# this index of the transcriptome should get reused across runs
my $transcriptome_index = "/var/tmp/Ce_transcriptome_bowtie2";

use Net::FTP;

system("mkdir -p $tmp_dir");
system("mkdir -p $output_dir");

# connect to server, and get file listing
my $ftp = Net::FTP->new($ftp_host) or die "couldn't connect to $ftp_host\n";
$ftp->login;
$ftp->binary;
$ftp->cwd($ftp_dir);
my @files = $ftp->ls;
$ftp->close;

foreach my $f (@files) {
  my $name = short_ee_filename($f);
  if (defined $name) {

    # skip things that are (presumably) already aligned
    next if (-e "$output_dir/$name.bam");

    print "[aligning $name]\n";

    # connect to server, and get this file
    $ftp = Net::FTP->new($ftp_host) or die "couldn't connect to $ftp_host\n";
    $ftp->login;
    $ftp->binary;
    $ftp->cwd($ftp_dir);
    $ftp->get($f, "$tmp_dir/combined.fastq.gz") or
      die("ftp failed: " . $ftp->message);
    $ftp->close;

    split_reads("$tmp_dir/combined.fastq.gz",
      "$tmp_dir/reads");

    tophat_align("$tmp_dir/reads", "$tmp_dir/tophat_out");

    system("mv $tmp_dir/tophat_out/accepted_hits.bam $output_dir/$name.bam");
# exit(0);    # XXX debugging
    unlink("$tmp_dir/combined.fastq.gz");
    unlink("$tmp_dir/reads.1.fastq");
    unlink("$tmp_dir/reads.2.fastq");
    system("\\rm -rf $tmp_dir/tophat_out");
  }
}

$ftp->quit;

# Checks whether a file name is an early embryonic timeseries,
# and if it is, computes a simplified filename.
# Args: a filename
# Returns: a simpler filename, or undef if it doesn't look like
#   an early embryonic timeseries.
sub short_ee_filename {
  my($f) = @_;

  if ($f =~ /total-RNA:Developmental-Stage=N2-EE-(50-\d+)\#Strain=N2\#.*(modENCODE_\d+):(\d+).fastq.gz/) {
    my $short_name = "EE_" . $1 . "_" . $2 . "_" . $3;
    return $short_name;
  }
  else {
    return undef;
  }
}

# Splits reads in a FASTQ file into two reads.
# presumably left and right.
# Args:
#   f - compressed FASTQ (.fastq.gz) file of
#     150-bp reads (which presumably is actually
#     two 75-bp reads)
#   out.base - base name for the output files
#     which will have ".1.fastq" and ".2.fastq"
#     tacked onto this.
# Side effects: writes two fastq files.
sub split_reads {
  my($f, $output_base) = @_;
  my $n = 76;

  open IN, "gunzip -c $f |" || die;
  open OUT1, ">$output_base.1.fastq" || die;
  open OUT2, ">$output_base.2.fastq" || die;

  while (!eof(IN)) {
    my $name1 = <IN>;
    chomp $name1;
    my $read = <IN>;
    chomp $read;
    my $name2 = <IN>;
    chomp $name2;
    my $qual = <IN>;
    chomp $qual;
    $name1 = substr($name1, 1, 100);
    $name2 = substr($name2, 1, 100);
    die if not ($name1 eq $name2);

    print OUT1 "@" . $name1 . "\n";
    print OUT1 substr($read, 0, $n) . "\n";
    print OUT1 "+" . $name1 . "\n";
    print OUT1 substr($qual, 0, $n) . "\n";

    print OUT2 "@" . $name1 . "\n";
    print OUT2 substr($read, $n, $n) . "\n";
    print OUT2 "+" . $name1 . "\n";
    print OUT2 substr($qual, $n, $n) . "\n";
  }

  close IN;
  close OUT1;
  close OUT2;
}

# Aligns one file using TopHat.
# Args:
#   fastq_base - the base name of the file to align
#   output_dir - directory in which to write output
# Side effects: runs TopHat, puts results in file.
sub tophat_align {
  my($fastq_base, $output_dir) = @_;

  system("$tophat2_bin -G $gene_gtf " .
    "--transcriptome-index $transcriptome_index " .
#    "--transcriptome-only --transcriptome-index $transcriptome_index " .
#    "--segment-length 38 " .
    "--num-threads 5 " . # --no-novel-juncs "
    "--output-dir $output_dir $bwt_index " .
    $fastq_base . ".1.fastq " . $fastq_base . ".2.fastq");
}

