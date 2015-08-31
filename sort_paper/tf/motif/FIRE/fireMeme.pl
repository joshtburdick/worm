#!/usr/bin/perl -w
# Converts FIRE's output to MEME format.

use strict;

# Converts one motif to MEME format.
# Args:
#   motif.name - name of the motif
#   m - a motif as a regexp (as returned by FIRE)
# Returns: that motif, in MEME format. (Note that the
# "nsites" and "E" numbers here are not correct.)
sub FIRE_to_meme {
  my($motifName, $m) = @_;
  my %r = ("A"=>"", "C"=>"", "G"=>"", "T"=>"");

  foreach my $letters (split /([ACGT]|\[[ACGT]+\])/, $m) {

    # count up which letters are present
    my %counts = ();
    my $s = 0;
    if (!($letters eq "")) {
      foreach my $b (qw/A C G T/) {
        if ($letters =~ $b) {
          $counts{$b} = 1;
          $s++;
        }
      }

      # tack onto results
      foreach my $b (qw/A C G T/) {
        $r{$b} .= (defined $counts{$b} ? ($counts{$b} / $s) : 0) . " ";
      }
    }
  }

  return join "\n",
    ("MOTIF $motifName", "", $r{"A"}, $r{"C"}, $r{"G"}, $r{"T"}, "", "");
}



print FIRE_to_meme("testMotif", "[ACG]ATCG[AG][AC]ATG[ACGT]");


