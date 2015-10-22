#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

foreach my $line (<FH>) {
	chomp $line;
	my $command = "perl get_gene.pl -genome mm10 -method gene -strains reference, balbcj, nodshiltj, pwkphj, spreteij, wsbeij -homo -transcript " . $line . " > " . $line . ".txt";
	print $command . "\n";
	`$command`;
}
close FH;
