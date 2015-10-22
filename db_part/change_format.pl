#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";
open OUT, ">$ARGV[0].vcf";

my @split;
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC57BL6\tS129P2\n";
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	print OUT $split[0] . "\t" . $split[1] . "\t.\t" . $split[2] . "\t" . $split[3] . "\t.\t.\t.\t.\t0/0\t1/1\n"; 
}
close FH;
close OUT;
