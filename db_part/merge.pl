#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";
my @split;
open OUT, ">$ARGV[0].vcf";

print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC57BL6\tBALBc\tS129P2\n";

my $balb;
my $s129;
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if($split[7] eq "\\N") {
		print OUT $split[0] . "\t" . $split[1] . "\t.\t";
		print OUT $split[2] . "\t" . $split[3];
		$balb = $split[3];
		$s129 = "NA";
	} elsif($split[3] eq "\\N") {
		print OUT $split[5] . "\t" . $split[6] . "\t.\t";
		print OUT $split[7] . "\t" . $split[8];
		$balb = "NA";
		$s129 = $split[8];
	} else {
		print OUT $split[0] . "\t" . $split[1] . "\t.\t";
		if($split[2] eq $split[7]) {
			$balb = $split[3];
			$s129 = $split[8];
			print OUT $split[2] . "\t" . $split[3];
			if($balb ne $s129) {
				print OUT "," . $split[8];
			}
		} else {
			if(length($split[2]) > length($split[7])) {
				print OUT $split[2] . "\t" . $split[3] . "," . $split[8] . substr($split[2], length($split[7]));
				$balb = $split[3];
				$s129 = $split[8] . substr($split[2], length($split[7]));
			#	print $split[2] . "\t" . $split[3] . "\t" . $split[2] . "\t" . $split[8] . substr($split[2], length($split[7])) . "\n\n";
			} else {
				print OUT $split[7] . "\t" . $split[3] . substr($split[7], length($split[2])) . "," . $split[8];
				$balb = $split[3] . substr($split[7], length($split[2]));
				$s129 = $split[8];
			#	print $split[7] . "\t" . $split[3] .substr($split[7], length($split[2])) . "\t" . $split[7] . "\t" . $split[8] . "\n\n";
			}
		}
	}
	print OUT "\t.\t.\t.\t.\t0/0\t";
	if($balb eq "NA") {
		print OUT "0/0\t1/1\n";
	} elsif($s129 eq "NA") {
		print OUT "1/1\t0/0\n";
	} elsif($balb eq $s129) {
		print OUT "1/1\t1/1\n";
	} else {
		print OUT "1/1\t2/2\n";
	}
}
