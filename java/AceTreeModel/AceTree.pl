#!/usr/bin/perl -w
my $dir = `pwd`;
chomp $dir;
print $dir . "\n";

open(OUT,">LocalReferenceModel.xml") || die "Couldn't write to file $dir/LocalReferenceModel.xml: Check permissions\n";
open(IN, "ReferenceModel.xml") || die "Couldn't read $dir/ReferenceModel.xml\n";

foreach(<IN>){
    chomp;
    if(/(.+)\.(\/.+)/){
	print OUT $1 . $dir . $2;
	print OUT "\n";
    }else{
	print OUT;
	print OUT "\n";
    }
}

# jtb: altering this so that standard output is printed
system("java -jar -Xmx500m AceTree_recomp.jar LocalReferenceModel.xml");

