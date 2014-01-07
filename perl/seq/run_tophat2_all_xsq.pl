#!/usr/bin/perl -w
# Runs this on all XSQ files.

my $tmp_dir = "/media/disk2/jburdick/csfasta/";

my $output_dir = "/media/disk2/jburdick/tophat2/WS220/";

system("mkdir -p $output_dir");

if (1) {
  my $xsq_prefix = "/murrlab/seq/reads/Murray050912/xsq/Murray050912_Pool1Pool2_0";
  
  foreach my $i (1..6) {
    my $f = $xsq_prefix . $i . ".xsq";
      print "running on $f\n";
      system("../run_tophat2.pl $f " .
        " $tmp_dir $output_dir");
  }
}

if (1) {
  my $p = "/murrlab/seq/reads/Murray_52831_092812/result/";

  system("../run_tophat2.pl $p/lane1/Murray_52831_092812_01.xsq " .
    "$tmp_dir $output_dir");
  system("../run_tophat2.pl $p/lane2/Murray_52831_092812_02.xsq " .
    "$tmp_dir $output_dir");
  system("../run_tophat2.pl /home/jburdick/xsq/Murray_52831_092812_03.edit.xsq " .
    "$tmp_dir $output_dir");
  system("../run_tophat2.pl /home/jburdick/xsq/Murray_52831_092812_04.edit.xsq " .
    "$tmp_dir $output_dir");
  system("../run_tophat2.pl $p/lane5/Murray_52831_092812_05.xsq " .
    "$tmp_dir $output_dir");
  system("../run_tophat2.pl $p/lane6/Murray_52831_092812_06.xsq " .
    "$tmp_dir $output_dir");
}

