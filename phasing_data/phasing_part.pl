#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my @split;
my $f = 0;
my $start;
my $stop;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if($f == 0) {
		$start = $split[3] - 50;
		$stop = $split[3] + 50;
		print $line . "\n";
		$f++;
		next;
	}
	print STDERR $start . "\t" . $stop . "\n";
	if($split[3] > $start && $split[3] < $stop) {
		print $line . "\n";
	}
	
}
