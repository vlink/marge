#!/usr/bin/perl -w

use strict;

if(@ARGV < 2) {
	print STDERR "Usage: perl $0 <snp file> <motif file>\n";
	exit;
}

open FH, "<$ARGV[0]";
my %snps;
my @split;
my @header;
my $f = 0;

foreach my $line (<FH>) {
	if($f == 0) {
		@header = split('\t', $line);
		$f++;
		next;
	}
	chomp $line;
	@split = split('\t', $line);
	$snps{$split[4]} = $line;
}
close FH;

open FH, "<$ARGV[1]";
my %output;
my @snps;
my $half;
my $motif_start;
my $motif_end;
my @mut;
my @pos;
my @motif_name;
my @strain;
my %all;

foreach my $line (<FH>) {
	chomp $line;
		@split = split('\t', $line);
#	print $line . "\n";
#	print $snps{$split[0]} . "\n";
	#Check if there are mutations at all
	@snps = split('\t', $snps{$split[0]});
#	print @snps . "\n";
	if(@snps > 7) {
		#keep going
		#Calculate half
		$half = ($snps[3] - $snps[2])/2;
		$motif_start = $half + $split[1];
		if($split[4] eq "-" && $split[1] < 0) {
			$motif_start = $motif_start - (length($split[2])-1);
		}
		if($split[4] eq "-" && $split[1] > 0) {
			$motif_start = $motif_start - (length($split[2])-1);
		}
		$motif_start = $snps[2] + $motif_start;
		$motif_end = $motif_start + length($split[2])-1;
	#	print $motif_start . " - " . $motif_end . "\n";
		for(my $i = 8; $i < @snps; $i++) {
		#	print $snps[$i] . "\n";
			@mut = split(";", $snps[$i]);		
			for(my $j = 0; $j < @mut; $j++) {
				@pos = split(":", $mut[$j]);
			#	print $pos[0] . "\n";
				if($pos[0] >= $motif_start && $pos[0] <= $motif_end) {
				#	print "OVERLAP!!!!\n";
					@motif_name = split('\(', $split[3]);
					my @a = split(":", $motif_name[0]);
					$motif_name[0] = $a[-1];
				#	print $motif_name[0] . "\n";
					@strain = split('\(', $header[$i]);
				#	print $strain[0] . "\n";
					$output{$split[0]}{$motif_name[0] . "_" . $strain[0]} = 1;
					$all{$motif_name[0] . "_" . $strain[0]} = 1;
				}
			}
		} 
#		$m_from_seq = substr($gapless{$s}, $motif_start, length($seq_motif{$motif . "_" . $t . "_" . $s}));
#		for(my $l = $start; $l < $motif_start; $l++) {
#			$array[$num] .= " ";
#		}
#		if($motif_start < $start) {
#			$array[$num] .= substr($seq_motif{$motif . "_" . $t . "_" . $s}, $start - $motif_start);
#			$start = $motif_start + length(substr($seq_motif{$motif . "_" . $t . "_" . $s}, $start - $motif_start));
#		} else {
#			$array[$num] .= $m_from_file;
#			$start = $motif_start + length($m_from_file) - 1;
#		}
#		$start = $motif_start + length($seq_motif{$motif . "_" . $t . "_" . $s});

#	exit;
	}
	$output{$split[0]}{'peak'} = $snps[0] . "\t" . $snps[1] . "\t" . $snps[2] . "\t" . $snps[3] . "\t" . $snps[4] . "\t" . $snps[5] . "\t" . $snps[6];
}

$f = 0;
foreach my $peak (keys %output) {
	if($f == 0) {
		print "ID\tchr\tstart\tend\torignalPeak\tmerge_num\tStrand";
		foreach my $ms (sort {$a cmp $b} keys %all) {
			print "\t" . $ms;
		}
		print "\n";
		$f++;
	}
	print $output{$peak}{'peak'};
	foreach my $ms (sort {$a cmp $b} keys %all) {
		print "\t";
		if(exists $output{$peak}{$ms}) { print "X"; }
	}
	print "\n";
}
