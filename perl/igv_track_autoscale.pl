#!/usr/bin/perl -w
# Alters an IGV session file so that:
# - "autoscale" is false
# - scale is set to some constant (currently 200)

while (<>) {
  s/autoscale=\"true\"/autoscale=\"false\"/;
  if (/^(.*)maximum=\"[^ ]+\" minimum=\"[^ ]+\"(.*)$/) {
    $_ = "$1 maximum=\"200\" minimum=\"0\" $2";
  }

  print;
}

