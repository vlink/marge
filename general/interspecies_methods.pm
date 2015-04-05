#!/usr/bin/perl -w

package interspecies_methods;
use strict;

sub split_chain_file{
	$_ = 0 for my ($start1, $start2, $count1, $count2, $chr1, $chr2);
	$_ = () for my (@save, @split, @chain);
	my $chain;

	open FH, "<$_[0]";
	my $gap = $_[1];
	open OUT, ">$_[2]";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 5) eq "chain") {
			@split = split('\s+', $line);
			@chain = @split;
			$start1 = $split[5];
			$start2 = $split[10];
			$count1 = 0;
			$count2 = 0;	
			next;
		} 
		@split = split('\t', $line);
		if($line eq "") {
			next;
		}
		if(@split == 1) {
			print OUT "$chain[0] $chain[1] $chain[2] $chain[3] $chain[4] $start1 $chain[6] $chain[7] $chain[8] $chain[9] $start2 $chain[11] $chain[12]\n";
			foreach my $l (@save) {
				print OUT $l . "\n";
			}
			print OUT $split[0] . "\n\n";
			@save = ();
			next;
		}
		if($split[1] > $gap || $split[2] > $gap) {
			print OUT "$chain[0] $chain[1] $chain[2] $chain[3] $chain[4] $start1 " . ($start1 + $count1 + $split[0]) . " $chain[7] $chain[8] $chain[9] $start2 " . ($count2 + $start2 + $split[0]) . " $chain[12]\n";
			$start1 = $start1 + $count1 + $split[1] + $split[0];
			$start2 = $start2 + $count2 + $split[2] + $split[0];
			foreach my $l (@save) {
				print OUT $l . "\n";
			}
			print OUT $split[0] . "\n\n";
			@save = ();
			$count1 = 0;
			$count2 = 0;
		} else {
			$count1 += $split[0] + $split[1];
			$count2 += $split[0] + $split[2];
			push(@save, $line);
		}
	}
}
1;
