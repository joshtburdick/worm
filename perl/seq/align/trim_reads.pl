#!/usr/bin/perl -w
# Trims some reads.

use strict;

my $a = "../03_F21D5.9_2e6";

# Trims csfasta.
sub trim_csfasta {
  my($in, $out, $n) = @_;

  open IN, "<$in" || die;
  open OUT, ">$out" || die;

  while (<IN>) {
    print OUT;
    $_ = <IN>;
    chomp;
    $_ = substr($_, 0, length($_) - $n);
    print OUT "$_\n";
  }

  close IN;
  close OUT;
}

# Trims qualities.
sub trim_qual {
  my($in, $out, $n) = @_;

  open IN, "<$in" || die;
  open OUT, ">$out" || die;

  while (<IN>) {
    print OUT;
    $_ = <IN>;
    chomp;
    my @a = split /\s+/;
    splice @a, @a - $n;
    print OUT join " ", @a;
    print OUT "\n";
  }

  close IN;
  close OUT;
}

trim_csfasta("$a/F3.csfasta", "F3.csfasta", 8);
trim_qual("$a/F3.qual", "F3.qual", 8);

trim_csfasta("$a/F5-RNA.csfasta", "F5-RNA.csfasta", 6);
trim_qual("$a/F5-RNA.qual", "F5-RNA.qual", 6);


