#!/usr/bin/perl -w
# Copies files for the cluster browser to one location.

use strict;

my $cwd = `pwd`;

# where to write the output
my $output_dir = "/home/jburdick/build/sortWeb/";

system("mkdir -p $output_dir");

# copy over the cluster browser
foreach my $f ("js/SvgHeatmap.js", "js/SvgLabel.js", "js/vectorMath.js",
  "clusterBrowse/clusterBrowse.html", "clusterBrowse/clusterBrowse.js",
  "clusterBrowse/data.js") {
  print("$f\n");
  system("(cd ..; cp --parents $f $output_dir; cd clusterBrowse )");
}

# special case: the "help" file
system("cp clusterBrowseHelp.html $output_dir/index.html");

# copy over per-cluster annotation
system("cd ~/gcb/git/sort_paper/plot/web/; cp --recursive clusters $output_dir");



