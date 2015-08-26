#!/usr/bin/perl -w
# Runs an R script, and saves the resulting
# definitions in an .Rdata file.

my $sourceFile = $ARGV[0];

# check that name looks like an R source file
die if !($sourceFile =~ /\.r$/);

my $rdataFile = $sourceFile;
$rdataFile =~ s/\.r$/\.RData/;

my $cmd = "R --quiet --no-save -e " .
  "'source(\"$sourceFile\"); save.image(file=\"$rdataFile\")'" .
  " > /dev/null";
system $cmd;

print "[built $rdataFile]\n";

