#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my @split;
my $count = 0;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t' , $line);
	if(@split == 3) {
		if($split[1] > 25) {
			$count++;
			next;
		}
		if($split[2] > 25) {
			$count++;
		}
	}
}
print "Additional splits\n";
print $count . "\n";
