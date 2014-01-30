#!/usr/bin/perl -w
# Makes merged BAM files.

use strict;

# Gets experiment names.
sub get_experiment_names {
  my %h = ();
  open IN, "</home/jburdick/gcb/git/unmix/seq/quant/experimentNames.tsv";
  $_ = <IN>;
  while (<IN>) {
    chomp;
    s/\"//g;
    my($id, $name) = split /\t/;
    $name =~ s/\//_/g;
    $h{$id} = $name;
  }
  close IN;

  return \%h;
}

my %experiment_names = %{ get_experiment_names() };
# print join " ", %experiment_names;

# Makes a merged list of all the .bam files.
sub make_merged_list {
  my($input_dir, $output_dir) = @_;

  system("mkdir -p $output_dir");

  my %file_by_gene = ();
  foreach my $f (sort (glob "$input_dir/*.bam")) {
    print $f;
    print "\n";
    die if not $f =~ /$input_dir\/(?:0\d_)?(.+)\.bam$/;
    my $name = $1;

    # convert to a human-readable name
    my $name1 = $experiment_names{$name};
print "$name\n";
    die if not defined $name1;

    push @{ $file_by_gene{$name1} }, $f;
  }

  foreach my $name (keys %file_by_gene) {
    open OUT, ">$output_dir/$name.bam.list" || die;
    foreach my $file (sort @{ $file_by_gene{$name}}) {
      print OUT "$file\n";
    }
    close OUT;
  }
}

# make_merged_list("/murrlab/seq/tophat2/Murray050912/strand_flip/",
#   "/murrlab/seq/igv/expression/Murray050912/");
my $output_dir = "/murrlab/seq/igv/expression/embryo_FACS/";

make_merged_list("/murrlab/seq/tophat2/WS220_20140111/20110922",
  $output_dir);
make_merged_list("/murrlab/seq/tophat2/WS220_20140111/Murray050912",
  $output_dir);
make_merged_list("/murrlab/seq/tophat2/WS220_20140111/Murray050912",
  $output_dir);

