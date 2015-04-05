#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my @split;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\s+', $line);
	print $line . "\t" . ($split[6] - $split[5]) . "\n";
}
