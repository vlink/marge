#!/usr/bin/perl -w

use strict;
require '../db_part/database_interaction.pm';

open FH, "<$ARGV[0]";

my @strains;
$strains[0] = "reference";
$strains[1] = "balbcj";
$strains[2] = "nodshiltj";
$strains[3] = "spreteij";
my $genome = "mm10";
my $homo = 1;
my $html = 0;
database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");


my @split;
my %mut;
my @m;
my @pos;
my $f = 0;
my @order;
my $count = 0;
my %seen = ();

foreach my $line (<FH>) {
	if($f == 0) { $f++; next; }
	chomp $line;
	@split = split('\t', $line);
	my $k = () = $line =~ /chr/g;
	my $exp = (@split - 1 - $k)/$k;
	%seen = ();
	for(my $i = 1; $i < @split; $i = $i + $exp + 1) {
		$split[$i] =~ s/ //g;
		if(exists $seen{$split[$i]}) {
			next;
		}
		$seen{$split[$i]} = 1;
	}
	foreach my $s (keys %seen) {
		my @a = split(":", $split[1]);
		my @b = split("-", $a[1]);
		my $seqs = database_interaction::get_genomic_seq($b[0], $b[1], substr($a[0], 3), "+", 0, 0, 0, 0);
		my @k = split('\n', $seqs);
		for(my $i = 0; $i < @k; $i++) {
			print ">" . $k[$i] . "_" . $s . "\n";;
			$i++;
			print $k[$i] . "\n";
		}
	}
	exit;
	print $split[0] . "\t" . $split[8] . "\t" . $split[9] . "\t" . $split[10] . "\t" . $split[11] . "\n"; 
	print $split[1] . "\n";
	for(my $i = @split - 3; $i < @split; $i++) {
		@m = split(";", $split[$i]);
		for(my $j  = 0; $j < @m; $j++) {
			@pos = split("\:", $m[$j]);
			$mut{$pos[0]} = 1;
		}
	}
	$count = 0;
	@order = ();
	foreach my $m (sort {$a <=> $b} keys %mut) { 
		$order[$count] = $m;
		$count++;
	}
	my $count2 = 0;
	for(my $i = @split - 3; $i < @split; $i++) {
		@m = split(";", $split[$i]);
		for(my $j = 0; $j < @m; $j++) {
			@pos = split(':', $m[$j]);
			while($order[$count2] != $pos[0]) {
				print " ";
				$count2++;
			}
			print "|";
			$count2++;
		}
		print "\n";
		$count2 = 0;
	}
	%mut = ();
}
