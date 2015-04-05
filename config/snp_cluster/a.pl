#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my $count = 1;
my @split;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	for(my $i = 2; $i < @split; $i++) {
		print "chr" . $split[1] . "\t" . $split[$i] . "\t" . ($split[$i] + 1) . "\tcluster" . $count . "\t1\t+\n";
	}
	$count++;
}
