#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my %seen;
my @split;
my $print = 0;
foreach my $line (<FH>) {
	if(substr($line, 0, 1) eq ">") {
		$print = 0;
		@split = split('\t', $line);
		if(exists $seen{$split[1]}{$split[0]}) {
			$print = 1;
		} else {
			$seen{$split[1]}{$split[0]} = 1;
		}
	}
	if($print == 0) {
		print $line;
	}
}
