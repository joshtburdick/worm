#!/usr/bin/perl -w
# Converts from MEME format to TRANSFAC format.
# Reads from the given file, writes to STDOUT.
# TRANSFAC format is as described at 
# http://www.benoslab.pitt.edu/stamp/help.html

convert_MEME_to_TRANSFAC($ARGV[0]);

# Converts one file.
sub convert_MEME_to_TRANSFAC {
  my($meme_input) = @_;

  open IN, "<$meme_input" or die;

  # skip past header
  $_ = <IN>;
  die if not /^MEME/;
  while (not ($_ =~ /Background letter frequencies/)) {
    $_ = <IN>;
  }
  $_ = <IN>;

  while ($_ = <IN>) {
    # note that some files contain spaces after this identifier
    if ($_ =~ /^MOTIF (\S+)/) {
      my $id = $1;
      $_ = <IN>;
      $_ = <IN>;

      # skip past the log-odds matrix, if it's there,
      # until letter-probability matrix
      if (/^log-odds matrix/) {
        while ($_ = <IN>) {
          last if /^letter-probability matrix/;
        }        
      }

# ??? some files need this; not sure all of them do

#      if not (/letter-probability matrix: alength= (\d+) w= (\d+) nsites= (\d+) E= /) {
#        $_ = <IN>;
#      }

      die if not /letter-probability matrix: alength= (\d+) w= (\d+) nsites= (\d+) E= /;
      my $w = $2;
      my $nsites = $3;
#      my $e = $4;
      print"DE $id\n";

      # copy lines in motif; presumably MotifSuite normalizes these
      foreach my $n (1..$w) {
        $_ = <IN>;
        my @a = split /\s+/, " $_";
        print (join "\t", $n, $a[1], $a[2], $a[3], $a[4]);
        print "\n";
      }

      print "XX\n";
    }
  }

  close IN;
}

