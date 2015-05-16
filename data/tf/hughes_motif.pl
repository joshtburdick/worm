#!/usr/bin/perl -w
# Imports motifs from the Hughes lab into MEME format.
# For now, just using the filenames as motif names.

my $input_dir = "~/gcb/data/tf/hughes/Ce_1.01/pwms_all_motifs/";

my $output_file = "hughes_Ce_1.01.meme";

# Translates one file to the format MEME expects.
sub as_meme {
  my($f, $name) = @_;
  my @r = ();

  open IN, "<$f" || die;
  $_ = <IN>;

  die if not ($_ =~ /Pos.*A.*C.*G.*T/);
  while (<IN>) {
    chomp;
    my @a = split /\s+/;
    shift @a;
    # XXX a little normalization
    my $P = $a[0] + $a[1] + $a[2] + $a[3];

    push @r, (join " ", $a[0]/$P, $a[1]/$P, $a[2]/$P, $a[3]/$P);
  }
  close IN;

  my $w = @r;

  # handle this special case
  if ($w < 1) {
    return "";
  }

  # XXX note that nsites here is totally made up
  my $s = "MOTIF $name\n\nletter-probability matrix: alength= 4 w= $w nsites= 100 E= 1\n";

  return $s . (join "\n", @r) . "\n\n";
}


open OUT, ">$output_file" || die;
print OUT <<END;
MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from /home/jburdick/gcb/git/tf/motif/Ce_WS220.order2markov.txt):
A 0.32700 C 0.17300 G 0.17300 T 0.32700

END
  ;

foreach my $f (<"$input_dir"/*.txt>) {
  die if not ($f =~ /\/([^\/]+)\.txt$/);
  my $name = $1;

  print OUT as_meme($f, $name);
}

close OUT;

