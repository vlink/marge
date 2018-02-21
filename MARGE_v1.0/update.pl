#!/usr/bin/env perl
use strict;
use warnings;

# Copyright Verena M. Link <vlink@ucsd.edu>
# 
# This file is part of MARGE
#
# MARGE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MARGE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


use Getopt::Long;
use config;
use Data::Dumper;

$_ = () for my(@split, @MARGE_folders, @motif_bg, @motif_analysis_folders, %save_output, @tmp, $list);
$_ = 0 for my($MARGE_folders_list, $motif_bg_list, $motif_window, $motif_analysis_list, $software_update);
$_ = "" for my($tmp_output, $tmp_error, $command, $MARGE_genome, $file);

my $url = "http://homer.ucsd.edu/MARGE/";
 
sub printCMD{
	print STDERR "Script to download already pre-processed data for MARGE\n";
	print STDERR "The script offers downloads for MARGE folders and individualzed genomes, as well as background scans for motifs and pre-processed files for motif analysis per strain (200bp)\n";
	print STDERR "\n\n";
	print STDERR "Usage:\n";
	print STDERR "-genome <genome>: The genome for which the motifs were scanned for (e.g. mm10/hg19)\n";
	print STDERR "-MARGE_folders <list of folders>: Downloads the list of all folders specified in the list. If all: downloads all available folders\n";
	print STDERR "-MARGE_folders_list: Lists all MARGE folders that are availabe for downloadn\n";
	print STDERR "-motif_bg <list of motifs>: Downloads all the motif background scans for the motifs listed here. If all: downloads all availabile motifs\n";
	print STDERR "-motif_bg_list: Lists all motif backgrounds for genome specified in motif_genome\n";
	print STDERR "-motif_analysis_folders <list of folders> : Downloads all the preparsed motif analysis backgrounds specified in the list . If all: downloads all availabe backgrounds\n";
	print STDERR "-motif_window <bp>: Motif window for which files should be downloaded (default in HOMER motif analysis - 200bp)\n";
	print STDERR "-motif_analysis_list: Lists all motif analysis background scans\n";
	print STDERR "-software_update: Checks for updates\n";
	exit;
}
#CHECK FOR UPDATES

if(@ARGV < 1) {
        &printCMD();
}

GetOptions(	"MARGE_folders=s{,}" => \@MARGE_folders,
		"MARGE_folders_list" => \$MARGE_folders_list,
		"genome=s" => \$MARGE_genome,
		"motif_bg=s{,}" => \@motif_bg,
		"motif_bg_list" => \$motif_bg_list,	
		"motif_analysis_folders=s{,}" => \@motif_analysis_folders,
		"motif_window=s" => \$motif_window,
		"motif_analysis_list" => \$motif_analysis_list,
		"software_update" => \$software_update)
	or die(&printCMD());

my $config = config::read_config();

#First check for all the lists
if($MARGE_folders_list == 1) {
	&list_MARGE_folders(1);
	exit;
}
if($motif_bg_list == 1) {
	&list_motif_bg(1);
	exit;
}
if($motif_analysis_list == 1) {
	&list_motif_analysis(1);
	exit;
}
if($software_update == 1) {
	print STDERR "We need to add this routine in here\n";
	exit;
}

if((@MARGE_folders >= 1 || @motif_bg >= 1 || @motif_analysis_folders >= 1) && $MARGE_genome eq "") {
	print STDERR "Please specify genome!\n\n";
	exit;
}

if(@MARGE_folders == 1) {
	@tmp = split(",", $MARGE_folders[0]);
	if(@tmp > @MARGE_folders) {
		@MARGE_folders = @tmp;
	}
}

for(my $i = 0; $i < @MARGE_folders; $i++) {
	$MARGE_folders[$i] =~ s/,//g;
	$MARGE_folders[$i] = uc($MARGE_folders[$i]);
}

if(@MARGE_folders >= 1) {
	$list = &list_MARGE_folders(0);
	if(@MARGE_folders == 1 && uc($MARGE_folders[0]) eq "ALL") {
		@MARGE_folders = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @MARGE_folders, $s;
		}
	}
	for(my $i = 0; $i < @MARGE_folders; $i++) {
		if(!exists $list->{$MARGE_genome}->{$MARGE_folders[$i]}) {
			print STDERR "Could not find " . $MARGE_folders[$i] . " for " . $MARGE_genome . "\n";
		} else {
			&download_file_MARGE($MARGE_genome, $MARGE_folders[$i]);
		}
	}
}

if(@motif_analysis_folders == 1) {
	@tmp = split(",", $motif_analysis_folders[0]);
	if(@tmp > @motif_analysis_folders) {
		@motif_analysis_folders = @tmp;
	}
}

for(my $i = 0; $i < @motif_analysis_folders; $i++) {
	$motif_analysis_folders[$i] =~ s/,//g;
	$motif_analysis_folders[$i] = uc($motif_analysis_folders[$i]);
}

if(@motif_analysis_folders >= 1) {
	if($motif_window == 0) {
		print STDERR "Motif window was not specified!\n";
		print STDERR "Set it to default - 200bp\n";
		$motif_window = 200;
	} 
	$motif_window =~ s/bp//g;
	$motif_window =~ s/BP//g;
	$motif_window = $motif_window . "bp";

	$list = &list_motif_analysis(0);
	if(@motif_analysis_folders == 1 && uc($motif_analysis_folders[0]) eq "ALL") {
		@motif_analysis_folders = ();
		foreach my $s (keys %{$list->{$MARGE_genome}->{$motif_window}}) {
			push @motif_analysis_folders, $s;
		}
	}
	for(my $i = 0; $i < @motif_analysis_folders; $i++) {
		if(!exists $list->{$MARGE_genome}->{$motif_window}->{$motif_analysis_folders[$i]}) {
			print STDERR "Could not find " . $motif_analysis_folders[$i] . " for " . $MARGE_genome . " with window " . $motif_window . "\n";
			print STDERR "Run this script with parameter -motif_analysis_list to see what is available\n";  
		} else {
			&download_motif_analysis($MARGE_genome, $motif_analysis_folders[$i], $motif_window);
		}
	}
}

if(@motif_bg == 1) {
	@tmp = split(",", $motif_bg[0]);
	if(@tmp > @motif_bg) {
		@motif_bg = @tmp;
	}
}

for(my $i = 0; $i < @motif_bg; $i++) {
	$motif_bg[$i] =~ s/,//g;
	$motif_bg[$i] = uc($motif_bg[$i]);
}

if(@motif_bg >= 1) {
	print STDERR "Waiting for data to come in!\n";
	exit;
}

sub download_file_MARGE{
	my $g = $_[0];
	my $f = $_[1];
	print STDERR "Downloading " . $f . "\n";
	$command = "wget " . $url . "/MARGE_folders/" . $g . "_" . $f . ".tar.gz";
	`$command`;
	#Unzip file
	print STDERR "Extracting and moving " . $f . "\n";
	$command = "tar xfvz " . $g . "_" . $f . ".tar.gz";
	`$command`;
	$command = "mv " . uc($f) . " " . $config->{'data_folder'};
	#Move file
	`$command`;
	$command = "rm " . $g . "_" . $f . ".tar.gz";
	if(-e $g . "_" . $f . ".tar.gz") {
		`$command`;
	}
}

sub download_motif_analysis{
	my $g = $_[0];
	my $f = $_[1];
	my $w = $_[2];
	if(!-e $config->{'data_folder'} . "/" . $f) {
		print STDERR "There is no data folder for $f!\n";
		print STDERR "Please download first the MARGE folder for this individual genome!\n";
		exit;
	}

	print STDERR "Downloading background for findMotifs for " . $f . " for " . $w . "\n";
	$command = "wget " . $url . "/motif_analysis/" . $g . "_preparsed_" . $f . "_" . $w . ".tar.gz";
	print $command . "\n";
	`$command`; 
	#Extract
	print STDERR "Extracting and moving " . $f . "\n";
	$command = "tar xfvz " . $g . "_preparsed_" . $f . "_" . $w . ".tar.gz";
	print $command . "\n";
	`$command`;
	if(!-e $config->{'data_folder'} . "/" . $f . "/preparsed") {
		print STDERR "Generating folder for preparsed data!\n";
		$command = "mkdir " . $config->{'data_folder'} . "/" . $f . "/preparsed";
		`$command`;
	}
	$command = "cp " . $g . "_preparsed_" . $f . "_" . $w . "/* " . $config->{'data_folder'} . "/" . $f . "/preparsed/";	
	print $command . "\n";
	`$command`;
	$command = "rm -rf " . $g . "_preparsed_" . $f . "_" . $w;
	if(-e $g . "_preparsed_" . $f . "_" . $w) {
		print $command . "\n";
		`$command`;
	}
	if(-e $g . "_preparsed_" . $f . "_" . $w . ".tar.gz") {
		$command = "rm " . $g . "_preparsed_" . $f . "_" . $w . ".tar.gz";
		print $command . "\n";
		`$command`;
	}
}

sub list_MARGE_folders {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/MARGE_folders/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", substr($line, 0, length($line) - 7));
		$save_output{$split[0]}{$split[1]} = 1;
	}
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print "Folders for genome " . $genome . "\n";
			foreach my $strain (sort {$a cmp $b} keys %{$save_output{$genome}}) {
				print "\t" . $strain . "\n";
			}
		}
	}
	return \%save_output;
}

sub list_motif_bg {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/motif_bg/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", substr($line, 0, length($line) - 4));
		$save_output{$split[0]}{$split[1]} = 1;
	}
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print "Motifs that were scanned for genome " . $genome . "\n";
			foreach my $strain (sort {$a cmp $b} keys %{$save_output{$genome}}) {
				print "\t" . $strain . "\n";
			}
		}
	}
	return \%save_output;
}

sub list_motif_analysis {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/motif_analysis/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", substr($line, 0, length($line) - 7));
		$save_output{$split[0]}{$split[3]}{$split[2]} = 1;
	}
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print "Motifs that were scanned for genome " . $genome . "\n";
			foreach my $window (keys %{$save_output{$genome}}) {
				print "Window that was used for scan: " . $window . "\n";
				foreach my $strain (sort {$a cmp $b} keys %{$save_output{$genome}{$window}}) {
					print "\t\t" . $strain . "\n";
				}
			}
		}	
	}
	return \%save_output;
}
