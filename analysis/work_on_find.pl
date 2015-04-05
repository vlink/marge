#!/usr/bin/perl -w

use strict;

open FH, "<$ARGV[0]";

my %save;
my %save_pos;
my @split;
my @second;
my %seq;

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	@second = split("_", $split[0]);
	$seq{$split[3]} = $split[2];
	if(exists $save{$second[1]}{$split[3]}{$second[0]}) {
		$save{$second[1]}{$split[3]}{$second[0]}++;
		$save_pos{$second[1]}{$split[3]}{$second[0]} .= "," . $split[1];
	} else {
		$save{$second[1]}{$split[3]}{$second[0]} = 1;
		$save_pos{$second[1]}{$split[3]}{$second[0]} = $split[1];
	}
}
close FH;

my $first = 0;
my $print;
my $equal = 0;
my $val;
my %cands;

foreach my $pos (keys %save) {
	$first = 0;
	print $pos . "\n";
	foreach my $motif (keys $save{$pos}) {
		if($motif =~ /fly/) { next; }
		if($motif =~ /Yeast/) { next; }
		$equal = 0;
		foreach my $k (sort {$a cmp $b} keys $save{$pos}{$motif}) {
			if($save{$pos}{$motif}{$k} > 10 && $first == 0) { last; }
			if($first == 0) { 
				$print = $motif; 
			}
			$print .= "\t" . $k . "\t" . $save{$pos}{$motif}{$k};

			if($first == 0) {
				$val = $save{$pos}{$motif}{$k};
			} else {
				if($val != $save{$pos}{$motif}{$k}) { $equal = 1; }
			}
			$first++;
		}
		if($equal == 1) {
			$cands{$pos}{$motif} = 1;
			print $print . "\n";
		}
		$first = 0;
	}
	print "\n";
}

my @pos;
my @chr;
my $middle;
my $start;
my @array;
my @name;
my $num = 0;
foreach my $c (keys %cands) {
	@chr = split(",", $c);
	print $c . "\n";
	print $chr[1] . "\t" . $chr[2] . "\n";
	for(my $i = $chr[1]; $i < $chr[2]; $i++) {
		print "|";
	}
	print "\n";
	$middle = $chr[1] + (($chr[2] - $chr[1])/2);
	print "middle: " . $middle . "\n";
	foreach my $blub (keys $cands{$c}) {
		print $blub . "\n";
		foreach my $s (keys $save_pos{$c}{$blub}) {
			print $s . "\t" . $save_pos{$c}{$blub}{$s} . "\n";

			my @k = split(",", $save_pos{$c}{$blub}{$s});
			my %k = ();
			for(my $j = 0; $j < @k; $j++) {
				$k{$k[$j]} = 1;
			}
			$start = $chr[1];
			foreach my $kkk (sort {$a cmp $b} keys %k) {
				for(my $l = $start; $l < ($middle + $kkk); $l++) {
					$array[$num] .= " ";
				}
				for(my $l = 0; $l < length($seq{$blub}); $l++) {
					$array[$num] .= "#";
				}
				$start = ($middle + $kkk) + length($seq{$blub});
			}
			for(my $l = $start; $l < $chr[1]; $l++) {
				$array[$num] .= " ";
			}
			$name[$num] = $blub . " " . $s;
			$num++;
		}
	}
	while(length($array[0]) > 0) {
		for(my $i = 0; $i < @array; $i++) {
			print $name[$i] . "\t";
			if(length($array[$i]) < 100) {
				print $array[$i] . "\n";
				$array[$i] = "";
			} else {
				print substr($array[$i], 0, 100) . "\n";
				$array[$i] = substr($array[$i], 100);
			}
		}
	}
}

