#!/usr/bin/env perl
BEGIN {push @INC, '/gpfs/data01/glasslab/home/vlink/code/marge/bin'}
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
use Storable;
use Set::IntervalTree;
use config;
use general;
use analysis_tree;
use Data::Dumper; 

$_ = "" for my($file, $output, $data_dir, $genome_dir);
$_ = () for my(@strains, %peaks, %strand, @split, %tree, %lookup_strain, %last_strain, @tmp_split, %save_id);
$_ = 0 for my($hetero, $allele, $line_number, $id);

sub printCMD {
        print STDERR "\nUsage:\n";
        print STDERR "\t-ind <individuals>: Comma-separated list of individuals\n";
        print STDERR "\t-file <file>: File with genomic coordinates to pull the sequences\n";
	print STDERR "\t-output <file>: Name of the output files (Default: sequences.txt)\n";
	print STDERR "\t-id: Uses peak ID as identifier for sequences - can only be used when only one individual was specified\n";
	print STDERR "\t-hetero: Data is heterozygous (Default: homozygous)\n";
	print STDERR "\nAdditional parameters:\n";
	print STDERR "\t-data_dir <directory>: default defined in config\n";
	print STDERR "\t-genome_dir <directory>: default defined in config\n\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my %mandatory = ('-file' => 1, '-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(   	"file=s" => \$file,
		"ind=s{,}" => \@strains,
		"genome_dir=s" => \$genome_dir,
		"data_dir=s" => \$data_dir,
		"hetero" => \$hetero,
		"id" => \$id,
		"output=s" => \$output)
	or die(&printCMD());
#First step: Get the sequences for the peaks

#Set variables
if($hetero == 1) {
	$allele = 2;
} else {
	$allele = 1;
}
if($data_dir eq "") {
	$data_dir = config::read_config()->{'data_folder'};
}
if($genome_dir eq "") {
	$genome_dir = $data_dir;
}
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}
if($output eq "") {
	$output = "sequences.txt";
}

if(@strains > 1 && $id == 1) {
	$id = 0;
}

print STDERR "Saving peaks from " . $file . "\n";
open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#" || substr($line, 0, 6) eq "PeakID") {
		next;
	}
	@split = split('\t', $line);
	@tmp_split = split("_", $split[1]);
	$peaks{substr($tmp_split[0], 3)}{$split[2]} = $split[3];
	$strand{substr($tmp_split[0], 3)}{$split[2]} = $split[4];
	$save_id{$tmp_split[0] . "_" . $split[2] . "_" . $split[3]} = $split[0];
	#print STDERR $line . "\n";
	#print STDERR $tmp_split[0] . "\t" . $split[2] . "\t" . $split[3] . "\n\n";
	$line_number++;
}
close FH;
#for(my $a = 1; $a <= $allele; $a++) {
	#Read in strains data
	print STDERR "Loading shift vectors\n";
	for(my $i = 0; $i < @strains; $i++) {
		my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data_dir, "ref_to_strain");
		$tree{$strains[$i]} = $tree_ref;
		$lookup_strain{$strains[$i]} = $lookup;
		$last_strain{$strains[$i]} = $last;
	}
	#Get sequences for every peak
	if($id == 0) {
		analysis::get_seq_for_peaks($output, \%peaks, \@strains, $genome_dir, $allele, $line_number, 0, 0, \%tree, \%lookup_strain, \%last_strain, \%strand);
	} else {
		my $tmp = "tmp_" . rand(10) . ".txt";
		my ($seq, $long, $mut) = analysis::get_seq_for_peaks($tmp, \%peaks, \@strains, $genome_dir, $allele, $line_number, 0, 0, \%tree, \%lookup_strain, \%last_strain, \%strand);
		#	`rm $tmp`;
		open OUT, ">$output";
		#use Data::Dumper;
		#my $a = Dumper %save_id;
		#print STDERR $a . "\n\n\n\n";
		#$a = Dumper $seq;
		#print STDERR $a . "\n\n\n\n";
		foreach my $header (keys %{$seq}) {
			@split = split("_", $header);
			if(length($seq->{$header}) == 0) { next; }
			#print STDERR $header . "\n";
			#print STDERR $split[0] . "_" . $split[1] . "_" . $split[2] . "\n";
			#print STDERR $save_id{$split[0] . "_" . $split[1] . "_" . $split[2]} . "\t" . $seq->{$header} . "\n";
			print OUT $save_id{$split[0] . "_" . $split[1] . "_" . $split[2]} . "\t" . $seq->{$header} . "\n";
		}
		close OUT;
	}
#}
