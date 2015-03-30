#!/usr/bin/perl -w
# Runs MISO on the embryo timeseries.
# XXX for now, running in single-end mode (since we're only
# looking at the one splicing event)

use strict;

# XXX for now, using previous version of MISO

# where the MISO binaries are
my $miso_dir = "/home/jburdick/lib/python/misopy-0.5.3/misopy/";

# events to look for (FIXME change this to just the events we're
# interested in)
my $events_dir = "/home/jburdick/gcb/git/unmix/seq/quant/MISO/commonshortest/SE.WS220.index";

my $bam_dir = "/media/jburdick/disk2/data/ftp/bam/";

# loop through the BAM files
foreach my $f (<$bam_dir/*.bam>) {
  die if not $f =~ /\/([^\/]+).bam$/;
  my $name = $1;

  print "[running on $f]\n";

  system("cd /home/jburdick/gcb/git/unmix/seq/quant/MISO/; ./run_MISO.pl $f /home/jburdick/tmp/miso/$name");
}

# system("python $miso_dir/miso.py --run $events_dir $bam_file --output-dir output2/ --read-len 100");

