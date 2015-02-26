#!/usr/bin/perl
# 20150226, jtb: hacking this to only run on one chromosome
#   at a time

use Getopt::Long;
use Pod::Usage;
my $inf;
my $sam;
my $result;
my $chromosome_to_scan;

# GetOptions('i=s'=>\$inf,'s=s'=>\$sam,'o=s'=>\$result,'c=s'=>\$chromosome_to_scan);
GetOptions('i=s'=>\$inf,'o=s'=>\$result,'c=s'=>\$chromosome_to_scan);
if((!($inf))||(!($result))) {
    pod2usage();
}
my %judge;
my %geneinf;
my %genestr;
open(GENEINF,"$inf");
while(<GENEINF>){
    chomp();
    my @a=split("\t");
    $judge{$a[0]}++;
    if($judge{$a[0]}==1){
	$geneinf{$a[0]}=$_;
    }
    if($judge{$a[0]}>1){
	$genestr{$a[0]}=$genestr{$a[0]}.$_.";";
    }
}
close(GENEINF);
undef %judge;
my %samfile;
#my $tmpct=0;
my $fpkmnum=0;
print "Data Scanning...\n";
# open(SAM,"$sam");
while(<>){
    chomp();
    my @a=split("\t");
    #$tmpct++;
    #print "$tmpct\n";
    $fpkmnum=$fpkmnum+0.5;
    if ((!defined($chromosome_to_scan)) ||
        ($a[2] eq $chromosome_to_scan)) {
      $samfile{$a[2]}{$a[3]}{$a[0]}=$a[5];
    }
}
# close(SAM);
my $fpkmbase=0;
if($fpkmnum>0){
    $fpkmbase=1000000000/$fpkmnum;
}
print "Scan Finished!\n";
open(OUT,">$result");
print OUT "Gene_name\tLocation\tIsoform_name\tNumber_of_Fragments\tFPKM\tRelative_Abundance\n";
foreach my $gene (keys %geneinf){
    my %rdpair;
    my %rdchr;
    my %rdstart;
    my %rdend;
    my %rdcigar;
    my @agene=split("\t",$geneinf{$gene});
    foreach my $samchr (keys %samfile){
	if($samchr eq $agene[1]){
	    foreach my $samstart (keys %{$samfile{$samchr}}){
		if($samstart>=$agene[3] && $samstart<=$agene[4]){
		    foreach my $samname (keys %{$samfile{$samchr}{$samstart}}){
			$rdpair{$samname}++;
			$rdchr{$samname}=$samchr;
			$rdstart{$samname}{$rdpair{$samname}}=$samstart;
			$rdcigar{$samname}{$rdpair{$samname}}=$samfile{$samchr}{$samstart}{$samname};
		    }
		}
	    }
	}
    }
    my @rdlist=();
    foreach my $rdname  (keys %rdpair ){
	#print "gaga\t$rdname\n";
	if($rdpair{$rdname}>0){
	    @rdlist=(@rdlist,$rdname);
	}
    }
    print "$gene\t$agene[1]:$agene[3]-$agene[4]";
    if($#rdlist==-1){
	print "\twith 0 fragments\n";
	next;
    }
    #undef %rdpair;
    my $line=0;
    my $isochr;
    my @isolist;
    my @isostart;
    my @isoend;
    my %isoindex;
    @isolist=split(",",$agene[5]);
    $isochr=$agene[1];
    my @astr=split(";",$genestr{$gene});
    foreach my $str (@astr){
	my @bstr=split("\t",$str);
	@isostart=(@isostart,$bstr[3]+1);
	@isoend=(@isoend,$bstr[4]);
	my @aiso=split(",",$bstr[5]);
	for(my $i=0;$i<=$#aiso;$i++){
	    $isoindex{$isolist[$i]}{$line}=$aiso[$i];
	}
	$line++;
    }
    my %isolocindex;
    for(my $i=0;$i<=$#isolist;$i++){
	for(my $j=0;$j<=$#isostart;$j++){
	    if($isoindex{$isolist[$i]}{$j}==1){
		for(my $k=$isostart[$j];$k<=$isoend[$j];$k++){
		    $isolocindex{$isolist[$i]}{$k}=1;
		}
	    }
	}
    }
#####################
    my $rdct=0;
    my @rdlistf=();
    my %z;
    my %tranz;
    my %unique;
    my %fgstart;
    my %fgend;
#    print "2\n";
    foreach my $rdname (@rdlist){
	#print "gaga\t$rdname\n";
	my @tmp1=split(/[A-Z]/,$rdcigar{$rdname}{1});
	my @tmp2=split(/[A-Z]/,$rdcigar{$rdname}{2});
	$rdend{$rdname}{1}=$rdstart{$rdname}{1};
	$rdend{$rdname}{2}=$rdstart{$rdname}{2};
	for(my $i=0;$i<=$#tmp1;$i++){
	    $rdend{$rdname}{1}=$rdend{$rdname}{1}+$tmp1[$i];
	}
	for(my $i=0;$i<=$#tmp2;$i++){
	    $rdend{$rdname}{2}=$rdend{$rdname}{2}+$tmp2[$i];
	}
	if($rdpair{$rdname}==2){
	    if($rdstart{$rdname}{1}>=$rdstart{$rdname}{2}){
		$fgstart{$rdname}=$rdstart{$rdname}{2};
	    } else {
		$fgstart{$rdname}=$rdstart{$rdname}{1};
	    }
	    if($rdend{$rdname}{1}>=$rdend{$rdname}{2}){
		$fgend{$rdname}=$rdend{$rdname}{1};
	    } else {
		$fgend{$rdname}=$rdend{$rdname}{2};
	    }
	}
	if($rdpair{$rdname}==1){
	    $fgstart{$rdname}=$rdstart{$rdname}{1};
	    $fgend{$rdname}=$rdend{$rdname}{1};
	}
	if(($rdchr{$rdname} eq $isochr) && (($rdstart{$rdname}{1}>=$isostart[0] && $rdstart{$rdname}{1}<=$isoend[$#isoend]) || ($rdstart{$rdname}{2}>=$isostart[0] && $rdstart{$rdname}{2}<=$isoend[$#isoend]))){
	    my @type1=();
	    my @type2=();
	    my @b1=split("",$rdcigar{$rdname}{1});
	    my @b2=split("",$rdcigar{$rdname}{2});
	    for(my $i=0;$i<=$#b1;$i++){
		if($b1[$i] eq "N" || $b1[$i] eq "M" || $b1[$i] eq "I" || $b1[$i] eq "D"){
		    @type1=(@type1,$b1[$i]);
		}
	    }
	    for(my $i=0;$i<=$#b2;$i++){
		if($b2[$i] eq "N" || $b2[$i] eq "M" || $b2[$i] eq "I" || $b2[$i] eq "D"){
		    @type2=(@type2,$b2[$i]);
		}
	    }
	    if($#tmp1==$#type1 && $#tmp2==$#type2){
		my $base1=$rdstart{$rdname}{1}-1;
		my $base2=$rdstart{$rdname}{2}-1;
		my $ct1=0;
		my $ct2=0;
		my @indi1=();
		my @indi2=();
		my @comp=();
		for(my $i=0;$i<=$#isostart;$i++){
		    @indi1=(@indi1,0);
		    @indi2=(@indi2,0);
		}
		for(my $i=0;$i<=$#isolist;$i++){
		    @comp=(@comp,1);
		}
		for(my $i=0;$i<=$#type1;$i++){
		    if($type1[$i] eq "M" || $type1[$i] eq "I"){
			for(my $j=1;$j<=$tmp1[$i];$j++){
			    my $loc=$base1+$j;
			    for(my $k=0;$k<=$#isostart;$k++){
				if($indi1[$k]==1){
				    next;
				}
				if($loc>=$isostart[$k] && $loc<=$isoend[$k]){
				    $indi1[$k]=1;
				}
			    }
			}
			$base1=$base1+$tmp1[$i];
		    }
		    if($type1[$i] eq "N"){
			$base1=$base1+$tmp1[$i];
		    }		
		}
		my $unique_r1=0;
		for(my $i=0;$i<=$#indi1;$i++){
		    if($indi1[$i]==1){
			$unique_r1++;
		    }
		}
		my $unique1_r1=0;
		for(my $i=0;$i<=$#indi1;$i++){
		    if($indi1[$i]==1){
			$unique1_r1++;
			for(my $j=0;$j<=$#isolist;$j++){
			    if($isoindex{$isolist[$j]}{$i}==0){
				$comp[$j]=0;
			    }
			}
		    }
		    if($indi1[$i]==0 && $unique_r1>1 && $unique1_r1>=1 && $unique1_r1<$unique_r1){
			for(my $j=0;$j<=$#isolist;$j++){
			    if($isoindex{$isolist[$j]}{$i}==1){
				$comp[$j]=0;
			    }
			}
		    }
		}
		for(my $i=0;$i<=$#type2;$i++){
		    if($type2[$i] eq "M" || $type2[$i] eq "I"){
			for(my $j=1;$j<=$tmp2[$i];$j++){
			    my $loc=$base2+$j;
			    for(my $k=0;$k<=$#isostart;$k++){
				if($indi2[$k]==1){
				    next;
				}
				if($loc>=$isostart[$k] && $loc<=$isoend[$k]){
				    $indi2[$k]=1;
				}			    
			    }
			}
			$base2=$base2+$tmp2[$i];
		    }
		    if($type2[$i] eq "N"){
			$base2=$base2+$tmp2[$i];
		    }
		}
		my $unique_r2=0;
		for(my $i=0;$i<=$#indi2;$i++){
		    if($indi2[$i]==1){
			$unique_r2++;
		    }
		}
		my $unique1_r2=0;
		for(my $i=0;$i<=$#indi2;$i++){
		    if($indi2[$i]==1){
			$unique1_r2++;
			for(my $j=0;$j<=$#isolist;$j++){
			    if($isoindex{$isolist[$j]}{$i}==0){
				$comp[$j]=0;
			    }
			}
		    }
		    if($indi2[$i]==0 && $unique_r2>1 && $unique1_r2>=1 && $unique1_r2<$unique_r2){
			for(my $j=0;$j<=$#isolist;$j++){
			    if($isoindex{$isolist[$j]}{$i}==1){
				$comp[$j]=0;
			    }
			}
		    }
		}
		if($unique_r1>0 || $unique_r2>0){
		    $rdct++;
		    @rdlistf=(@rdlistf,$rdname);
		    for(my $i=0;$i<=$#isolist;$i++){
			$z{$rdname}{$isolist[$i]}=$comp[$i];
			$tranz{$rdname}{$isolist[$i]}=$comp[$i];
			if($comp[$i]==1){
			    $unique{$rdname}++;
			}
		    }
		}
		#print "@comp\n";
	    }
	}
    }
    if($rdct==0){
	next;
    }
    print "\twith $rdct fragments\n";
    my %dis_fg;
    my %dis_iso;
    my %dis_iso_known;
    my %denom_known;
    my %denom;
    for(my $i=0;$i<=$#isolist;$i++){
	foreach my $rdname (@rdlistf){
	    my $tmpstart;
	    my $tmpend;
	    for(my $j=0;$j<=$#isostart;$j++){
		if($j==0){
		    if($fgstart{$rdname}<$isostart[$j]){
			$tmpstart=$isostart[$j];
		    }
		    if($fgstart{$rdname}>=$isostart[$j] && $fgstart{$rdname}<=$isoend[$j]){
			$tmpstart=$fgstart{$rdname};
		    }
		    if($fgend{$rdname}>=$isostart[$j] && $fgend{$rdname}<=$isoend[$j]){
			$tmpend=$fgend{$rdname};
		    }
		}
		if($j>0 && $j<$#isostart){
		    if($fgstart{$rdname}<$isostart[$j] && $fgstart{$rdname}>$isoend[$j-1]){
			$tmpstart=$isostart[$j];
		    }
		    if($fgstart{$rdname}>=$isostart[$j] && $fgstart{$rdname}<=$isoend[$j]){
			$tmpstart=$fgstart{$rdname};
		    }
		    if($fgend{$rdname}<$isostart[$j] && $fgend{$rdname}>$isoend[$j-1]){
			$tmpend=$isoend[$j-1];
		    }
		    if($fgend{$rdname}>=$isostart[$j] && $fgend{$rdname}<=$isoend[$j]){
			$tmpend=$fgend{$rdname};
		    }
		}
		if($j==$#isostart){
		    if($fgstart{$rdname}>=$isostart[$j] && $fgstart{$rdname}<=$isoend[$j]){
			$tmpstart=$fgstart{$rdname};
		    }
		    if($fgstart{$rdname}<$isostart[$j] && $fgstart{$rdname}>$isoend[$j-1]){
			$tmpstart=$isostart[$j];
		    }
		    if($fgend{$rdname}>=$isostart[$j] && $fgend{$rdname}<=$isoend[$j]){
			$tmpend=$fgend{$rdname};
		    }
		    if($fgend{$rdname}<$isostart[$j] && $fgend{$rdname}>$isoend[$j-1]){
			$tmpend=$isoend[$j-1];
		    }
		    if($fgend{$rdname}>$isoend[$j]){
			$tmpend=$isoend[$j];
		    }		    
		}
	    }
	    $fgstart{$rdname}=$tmpstart;
	    $fgend{$rdname}=$tmpend;
	    if($z{$rdname}{$isolist[$i]}==1){
		for(my $j=0;$j<=$#isostart;$j++){
		    if($isoindex{$isolist[$i]}{$j}==1){
			if($tmpstart>=$isostart[$j] && $tmpstart<=$isoend[$j]){
			    if($tmpend>=$isostart[$j] && $tmpend<=$isoend[$j]){
				for(my $k=$tmpstart;$k<=$tmpend;$k++){
				    $dis_fg{$rdname}{$k}=1;
				    $dis_iso{$isolist[$i]}{$k}++;
				    $denom{$isolist[$i]}++;
				}
			    }
			    if($tmpend>$isoend[$j]){
				for(my $k=$tmpstart;$k<=$isoend[$j];$k++){
				    $dis_fg{$rdname}{$k}=1;
				    $dis_iso{$isolist[$i]}{$k}++;
				    $denom{$isolist[$i]}++;
				}
			    }
			}
			if($tmpstart<$isostart[$j]){
			    if($tmpend>=$isostart[$j] && $tmpend<=$isoend[$j]){
				for(my $k=$isostart[$j];$k<=$tmpend;$k++){
				    $dis_fg{$rdname}{$k}=1;
				    $dis_iso{$isolist[$i]}{$k}++;
				    $denom{$isolist[$i]}++;
				}
			    }
			    if($tmpend>$isoend[$j]){
				for(my $k=$isostart[$j];$k<=$isoend[$j];$k++){
				    $dis_fg{$rdname}{$k}=1;
				    $dis_iso{$isolist[$i]}{$k}++;
				    $denom{$isolist[$i]}++;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    my @alpha;
    my @alpha_old;
    for(my $i=0;$i<=$#isolist;$i++){
	$alpha[$i]=1;
    }
    my %isolength;
    for(my $i=0;$i<=$#isolist;$i++){
	for(my $j=0;$j<=$#isostart;$j++){
	    if($isoindex{$isolist[$i]}{$j}==1){
		$isolength{$isolist[$i]}=$isolength{$isolist[$i]}+$isoend[$j]-$isostart[$j]+1;
	    }
	}
	#print "$gene\t$isolist[$i]\t$isolength{$isolist[$i]}\n";
    }
    print "$gene\tComputing...\n";
    my %hfunction;
    my %numer;
    my %numer_known;
    for(my $i=0;$i<=$#isolist;$i++){
	foreach my $rdname(keys %dis_fg){
	    #print "$rdname\n";
	    foreach my $loci(keys %{$dis_fg{$rdname}}){
		if($isolocindex{$isolist[$i]}{$loci}==1){
		    $numer{$rdname}{$isolist[$i]}=$numer{$rdname}{$isolist[$i]}+$dis_iso{$isolist[$i]}{$loci};
		}
	    }
	    if($denom{$isolist[$i]}>0){
		$hfunction{$rdname}{$isolist[$i]}=$numer{$rdname}{$isolist[$i]}/$denom{$isolist[$i]};
	    }
	}
    }
    my $diff=1;
    my $ct=0;
    while($diff>0.0001){
	my ($up,$down);
	if($ct==0){
	    for(my $i=0;$i<=$#isolist;$i++){
		foreach my $rdname (@rdlistf){
		    ($up,$down) = (0,0);

# jtb: stopping if denominator is zero (or negative)
        my $denom = $isolength{$isolist[$i]}-300+1;
        if ($denom <= 0) {
          die "isoform " . $isolist[$i] . ": denom = " . $denom . "\n";
        }

		    if($z{$rdname}{$isolist[$i]}==1){
			$up = $alpha[$i]/($isolength{$isolist[$i]}-300+1);
		    }
		    for(my $j=0;$j<=$#isolist;$j++){
			if($z{$rdname}{$isolist[$j]}==1){
			    $down = $down + $alpha[$j]/($isolength{$isolist[$j]}-300+1); 
			}
		    }
		    if($down!=0){
			$tranz{$rdname}{$isolist[$i]}=$up/$down;
		    }
		}
	    }
	}
	if($ct>0){
	    for(my $i=0;$i<=$#isolist;$i++){
		foreach my $rdname (@rdlistf){
		    ($up,$down) = (0,0);
		    if($z{$rdname}{$isolist[$i]}==1){
			$up = $alpha[$i] * $hfunction{$rdname}{$isolist[$i]};
		    }
		    for(my $j=0;$j<=$#isolist;$j++){
			if($z{$rdname}{$isolist[$j]}==1){
			    $down = $down + $alpha[$j] * $hfunction{$rdname}{$isolist[$j]};
			}
		    }
		    if($down!=0){
			$tranz{$rdname}{$isolist[$i]}=$up/$down;
		    }
		}
	    }
	}
	for(my $i=0;$i<=$#isolist;$i++){
	    my $up1=0;
	    foreach my $rdname (@rdlistf){
		$up1=$up1+$tranz{$rdname}{$isolist[$i]};
	    }
	    $alpha_old[$i] = $alpha[$i];
	    $alpha[$i] = $up1/$rdct;
	}
	$ct++;
	for(my $t=0;$t<=$#alpha;$t++){
	    if($t==0){
		$diff = abs($alpha[$t] - $alpha_old[$t]);
	    }
	    if($t>0){
		if(abs($alpha[$t]-$alpha_old[$t])>$diff){
		    $diff = abs($alpha[$t] - $alpha_old[$t]);
		}
	    }
	}
	#print "$ct\t@alpha\t$diff\t$gene\n";
    }
    my $sumalpha;
    for(my $i=0;$i<=$#alpha;$i++){
	$sumalpha = $sumalpha + $alpha[$i];
    }
    if($sumalpha==0){
	next;
    }
    for(my $i=0;$i<=$#alpha;$i++){
	$alpha[$i]=$alpha[$i]/$sumalpha;
    }
    my $sumtheta=0;
    my @theta; 
    for(my $i=0;$i<=$#alpha;$i++){
	$sumtheta=$sumtheta+$alpha[$i]/$isolength{$isolist[$i]};
    }
    for(my $i=0;$i<=$#alpha;$i++){
	$theta[$i]=$alpha[$i]/($isolength{$isolist[$i]}*$sumtheta);
	my $fpkm= $fpkmbase * $rdct * $alpha[$i] / $isolength{$isolist[$i]};
	print OUT "$gene\t$agene[1]:$agene[3]-$agene[4]\t$isolist[$i]\t$rdct\t$fpkm\t$theta[$i]\n";
    }
    print "$gene\tComputation Finished With $ct Iterations!\n";
}

=head1 SYNOPSIS

    -s ---SAM_FILE is the mapping result file generated by mapping tools such as TopHat. It is in SAM format (see http://samtools.sourceforge.net/).

    -i ---ISOFORM_Compatible_Matrix is from the output of preprocess.pl

    -o ---Output_Result_FILE
