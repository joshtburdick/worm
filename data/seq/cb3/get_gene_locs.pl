#!/usr/bin/perl -w
# Gets C. briggsae gene locations, and upstream regions.

# first, filter the BED file of genes
open IN, "<../../../../data/UCSC/cb3/database/xenoRefFlat.bed" || die;
open OUT, "| bedtools sort -i - > cb3_Ce_genes.bed" || die;
open OUT2, "| bedtools sort -i - > cb3_Ce_gene_bounds.bed" || die;
while (<IN>) {
    chomp;
    my @a = split /\t/;

    # only include Ce orthologs
    my $gene = $a[3];
    next if not $gene =~ /CELE_(.*)/;
    $a[3] = $1;
    
    print OUT (join "\t", @a) . "\n";
    print OUT2 (join "\t", (map { $a[$_]; } (0..5))) . "\n";
}
close IN;
close OUT;
close OUT2;

# compute upstream regions
system("bedtools slop -i cb3_Ce_gene_bounds.bed -g cb3.fa.sizes " .
       " -i cb3_Ce_gene_bounds.bed -l 3000 -r 0 -s " .
       "| bedtools subtract -s -a - -b cb3_Ce_gene_bounds.bed " .
       "> cb3_Ce_gene_upstream_3kb.bed");

