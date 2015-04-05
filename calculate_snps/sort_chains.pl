#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my @split;
my %save;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\s+', $line);
	my @a = split("_", $split[2]);
	if(@a > 1) {
		next;
	}
	@a = split("_", $split[7]);
	if(@a > 1) { next; }
	$save{$split[2]}{$split[5]} = $line;
}	
close FH;

foreach my $key (sort {$a cmp $b} keys %save) {
	foreach my $k (sort {$a <=> $b} keys $save{$key}) {
		print $save{$key}{$k} . "\n";
	}
}
