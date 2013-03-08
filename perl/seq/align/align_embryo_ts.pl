#!/usr/bin/perl -w
# Aligns the modENCODE embryonic timeseries.

my $ftp_host = "data.modencode.org";
my $ftp_dir = "C.elegans/mRNA/RNA-seq/raw-seqfile_fastq/";

my $tmp_dir = "/var/tmp/fastq/";

my $output_dir = "/var/tmp/embr_timeseries/";

use Net::FTP;

system("mkdir -p $tmp_dir");
system("mkdir -p $output_dir");

# connect to server
my $ftp = Net::FTP->new($ftp_host) or die "couldn't connect to $ftp_host\n";
$ftp->login;
$ftp->binary;
$ftp->cwd($ftp_dir);


# get files
my @files = $ftp->ls;

my $i = 0;
foreach my $f (@files) {
  my $name = short_ee_filename($f);
  if (defined $name) {
    print "$f $name\n";
    $i++;
  }

  if ($i > 3) {
    exit(0);
  }
}

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


$ftp->quit;

