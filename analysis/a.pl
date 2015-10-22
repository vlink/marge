#!/usr/bin/perl -w

use strict;

my @split;

open FH, "<bg";
open REF, ">bg_reference.txt";
open NOD, ">bg_nod.txt";
open BALB, ">bg_balb.txt";
open SPRET, ">bg_spret.txt";
my $all = 26762;
my %bg;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	$bg{$split[1]}{$split[0]}{$split[2]} = $split[3];
}
close FH;

my $sum;

foreach my $strain (keys %bg) {
	foreach my $motif (keys %{$bg{$strain}}) {
		$sum = 0;
		print $strain . "\t" . $motif . "\n";
		foreach my $pos (keys %{$bg{$strain}{$motif}}) {
			$sum += $bg{$strain}{$motif}{$pos};
			if($strain eq "reference") {
				print REF $motif . "\t" . $pos . "\t" . $bg{$strain}{$motif}{$pos} . "\n";
			} elsif($strain eq "nodshiltj") {
				print NOD $motif . "\t" . $pos . "\t" . $bg{$strain}{$motif}{$pos} . "\n";

			} elsif($strain eq "balbcj") {
				print BALB $motif . "\t" . $pos . "\t" . $bg{$strain}{$motif}{$pos} . "\n";

			} elsif($strain eq "spreteij") {
				print SPRET $motif . "\t" . $pos . "\t" . $bg{$strain}{$motif}{$pos} . "\n";

			} else {
				print "FUCK\n";
			}
		}
		$bg{$strain}{$motif}{0} = $all - $sum;
			if($strain eq "reference") {
				print REF $motif . "\t" . 0 . "\t" . $bg{$strain}{$motif}{0} . "\n";
			} elsif($strain eq "nodshiltj") {
				print NOD $motif . "\t" . 0 . "\t" . $bg{$strain}{$motif}{0} . "\n";

			} elsif($strain eq "balbcj") {
				print BALB $motif . "\t" . 0 . "\t" . $bg{$strain}{$motif}{0} . "\n";

			} elsif($strain eq "spreteij") {
				print SPRET $motif . "\t" . 0 . "\t" . $bg{$strain}{$motif}{0} . "\n";

			} else {
				print "FUCK\n";
			}

	}
}

