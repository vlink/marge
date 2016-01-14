#!/usr/bin/perl -w

use strict;

if(@ARGV < 1) {
	print STDERR "Please add more options!\n";
	print STDERR "Path!\n";
	exit;
}

my @files = `ls $ARGV[0]/*plots*R`;
my @split;

my $strain1;
my $strain2;
my $tmnt;
my $AB;
my $current_motif;
my @file_name;

open R, ">tmp.R";

foreach my $f (@files) {
	chomp $f;
	@file_name = split("/", $f);
	@split = split("_", $file_name[-1]);
	if($split[-1] ne "KLA.R" && $split[-1] ne "notx.R") {
		next;
	}
	$strain1 = $split[1];
	$strain2 = $split[2];
	$tmnt = substr($split[4], 0, length($split[4]) - 2);
	$AB = $split[3];
	open FH, "<$f";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 4) eq "plot") {
			@split = split('\\\n', $line);
			$current_motif = substr($split[-1], 0, length($split[-1]) - 2);
		}		
		if(substr($line, 0, 3) eq "one" || substr($line, 0,3 ) eq "two") {
			print R $line . "\n";
		}
		if(substr($line, 0, 6) eq "p_both") {
			print R $line . "\n";
			print R "print(paste(\"" . $strain1 . "_" . $strain2 . "_" . $tmnt . "_" . $AB . "_" . $current_motif . "\", \": \", p_both, sep=\"\"))\n";  
		}
	}	
	close FH;
}

`Rscript tmp.R > tmp_output_rfile.txt`;

open FH, "<tmp_output_rfile.txt";

my @name;
my %save;
my %comp;

my $print_line;
my $sig = 0;

foreach my $line (<FH>) {
	chomp $line;
	@split = split(": ", substr($line, 5, length($line) - 6));
	@name = split("_", $split[0]);
	$save{$name[4]}{$name[0] . "_" . $name[1] . "_" . $name[2]} = $split[1];
	$comp{$name[0] . "_" . $name[1] . "_" . $name[2]} = 1;
}
close FH;

my $first = 0;
foreach my $motif (keys %save) {
	if($first == 0) {
		print "TF";
		foreach my $comp (keys %comp) {
			print "\t" . $comp;
		}
		print "\n";
		$first++;
	}
	$print_line = $motif;
	$sig = 0;
	foreach my $comp (keys %comp) {
		if(exists $save{$motif}{$comp}) {
			if($save{$motif}{$comp} < 1.0e-40) {
				$save{$motif}{$comp} = 1.0e-40;
			}
		#	if($save{$motif}{$comp} == 0) {
		#		$save{$motif}{$comp} = 1.0e-300;
		#	}
			$print_line .= "\t" . $save{$motif}{$comp};
			if($save{$motif}{$comp} < 0.001) {
				$sig = 1;
			}
		} else {
			$print_line .= "\t1";
		}
	}
	if($sig == 1) {
		print $print_line . "\n";
	}
}
