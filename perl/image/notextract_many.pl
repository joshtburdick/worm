#!/usr/bin/perl -w
# Runs notextract on many movies.

my $tools_dir = "/gpfs/fs0/l/murr/tools3/";

while (<>) {
  chomp;
  print "running notextract on $_\n";  
  system("$tools_dir/Notextract.sh $_");
}

