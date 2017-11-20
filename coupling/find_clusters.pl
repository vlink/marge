#!/usr/bin/perl -w
#
use strict;
if(@ARGV < 3) {
	print STDERR "Input file and chromosome and treatment\n";
	exit;
}
print STDERR $ARGV[0] . "\n";
open FH, "<$ARGV[0]";
my $i = 0;
my $f = 0;
$_ = () for my(@split, @save, @header);
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if($f == 0) { 
		for(my $j = 1; $j < @split; $j++) {
			$header[$j-1] = $split[$j];
		}
		$f++;
		next;
	}
	for(my $j = 1; $j < @split; $j++) {
		$save[$i][$j - 1] = $split[$j];
	}
	$i++;
}
close FH;
my $extend = 1;
open OUT, ">$ARGV[0]_coupled_peaks.txt";
open BED, ">$ARGV[0]_coupled_peaks.bed";
my @name;
my @name2;
my %number;
for(my $i = 0; $i < @save - 1; $i++) {
	$extend = 1;
	while($i+$extend < @save - 1 && $save[$i][$i+$extend] > 0.8) {
		$extend++;
	}
	if($extend > 4) {
		$extend--;
		@name = split("-", $header[$i]);
		@name2 = split("-", $header[$i+$extend-1]);
		#print STDERR "done: " . $name[2] . "\t" . $name2[2] . "\n";
		#print STDERR "number of peaks: " . $extend . "\n";
		print OUT $ARGV[1] . ":" . ($name[2] - 134) . "-" . ($name2[2] + 134) . "\t" . ($name2[2] - $name[2]) . "\t" . $extend . "\n";
		#print STDERR $ARGV[1] . ":" . ($name[2] - 134) . "-" . ($name2[2] + 134) . "\t" . ($name2[2] - $name[2]) . "\t" . $extend . "\n";
		print BED $ARGV[1] . "\t" . ($name[2] - 134) . "\t" . ($name2[2] + 134) . "\t" . $extend . "\n";
		#print STDERR "\n\n";
		$number{$extend}++;
		$i = $i + $extend;
	}
}
close BED;
close OUT;
open OUT, ">$ARGV[0]_coupled_distribution.txt";
foreach my $key (sort {$a <=> $b} keys %number) {
	print OUT $key . "\t" . $number{$key} . "\n";
}
close OUT;
