#!/usr/bin/env perl
use warnings;
use strict;
package general;

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


use Storable;
use Set::IntervalTree;
use Data::Dumper;

#Read in the lookup and last files for every strain
sub read_strains_data {
	my $strain = $_[0];
	my $data_dir = $_[1];	
	my $shift_direction = $_[2];
	$_ = () for my (%last, %lookup, %tree_strain, @split, @shift_files, @split_tree);
	$_ = "" for my($chr, $allele);
	my $file;
	#Check which shfit direction is specified
	if($shift_direction eq "ref_to_strain") {
		$file = $data_dir . "/" . $strain . "/last_shift_ref.txt";
		if(-e $file) {
			%last = %{retrieve($file)};
		}
	} else {
		$file = $data_dir . "/" . $strain . "/last_shift_strain.txt";
		if(-e $file) {
			%last = %{retrieve($file)};
		}
	}
	$file = $data_dir . "/" . $strain . "/lookup_table_chr.txt";
	if(-e $file) {
        	%lookup = %{retrieve($file)};
	}
	#Read in all shift files and store in Interval Tree
       	@shift_files = `ls $data_dir/$strain/*$shift_direction* 2> /dev/null`;
        foreach my $shift (@shift_files) {
                chomp $shift;
		@split = split("chr", $shift);
		@split = split('\.', $split[-1]);
		@split = split("_", $split[0]);
		$chr = $split[0];
		$allele = $split[2];
		if(exists $lookup{$chr}) {
			$chr = $lookup{$chr};
		}
		my $tree = Set::IntervalTree->new;
		open SHIFT, "<$shift";
		foreach my $line (<SHIFT>) {
			chomp $line;
			@split_tree = split('\t', $line);
			$tree->insert( {'shift'=>$split_tree[0] }, $split_tree[1], $split_tree[2]+1);
		}
		close SHIFT;
		$tree_strain{$chr}{$allele} = $tree;
        }
	return (\%tree_strain, \%last, \%lookup);	
}

sub read_strains_mut {
	my $strain = $_[0];
	my $data_dir = $_[1];
	$_ = () for my (%tree_mut, @mut_files, @split, @name, %last, %tree_pos);
	$_ = "" for my ($chr, $allele);
	@mut_files = `ls $data_dir/$strain/*mut 2> /dev/null`;
	foreach my $file (@mut_files) {
		chomp $file;
		@split = split("/", $file);
		@split = split('\.', $split[-1]);
		@split = split("_", $split[0]);
		$chr = substr($split[0], 3);
		if(@split < 2) {
			$allele = 1;
		} else {
			$allele = $split[2];
		}
		print STDERR "\tfor chromosome" . $chr . " - allele " . $allele . "\n";
		open FH, "<$file";
		my $tree = Set::IntervalTree->new;
		my $tree_pos = Set::IntervalTree->new;
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$tree->insert( {'mut'=> $split[0] }, $split[0], $split[0] + length($split[2]));
		       	$tree_pos->insert( {'pos'=> $split[1] . "->" . $split[2] }, $split[0], $split[0] + length($split[2]));	
		}
		close FH;
		$tree_mut{$chr}{$allele} = $tree;
		$tree_pos{$chr}{$allele} = $tree_pos;
	}
	$file = $data_dir . "/" . $strain . "/last_shift_ref.txt";
	if(-e $file) {
		%last = %{retrieve($file)};
	}
	return (\%tree_mut, \%tree_pos, \%last);
}

sub read_mutations_from_two_strains {
	my $strain_one = $_[0];
	my $strain_two = $_[1];
	my $data_dir = $_[2];
	$_ = () for my (%tree_mut, @mut_files, @split, @name, %last, %save);
	$_ = "" for my ($chr, $allele);
	#Read in all mutations from strain1 and save in hash
	@mut_files = `ls $data_dir/$strain_one/*mut 2> /dev/null`;

	foreach my $file (@mut_files) {
		chomp $file;
		@split = split("/", $file);
		@split = split('\.', $split[-1]);
		@split = split("_", $split[0]);
		$chr = substr($split[0], 3);
		if(@split < 2) {
			$allele = 1;
		} else {
			$allele = $split[2];
		}
		print STDERR "\tfor chromosome" . $chr . " - allele " . $allel . "\n";
		open FH, "<$file";
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$save{$chr}{$allele}{$split[0]} = length($split[2]);
		}
		close FH;
	}

	$file = $data_dir . "/" . $strain_one . "/last_shift_ref.txt";
	if(-e $file) {
		%last = %{retrieve($file)};
	}

	@mut_files = `ls $data_dir/$strain_two/*mut 2> /dev/null`;
	
	foreach my $file (@mut_files) {
		chomp $file;
		@split = split("/", $file);
		@split = split('\.', $split[-1]);
		@split = split("_", $split[0]);
		$chr = substr($split[0], 3);
		if(@split < 2) {
			$allele = 1;
		} else {
			$allele = $split[2];
		}
		print STDERR "\tFile two for chromosome " . $chr . "\n";
		open FH, "<$file";
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			if(exists $save{$chr}{$allele}{$split[0]}) {
				delete $save{$chr}{$allele}{$split[0]};
			} else {
				$save{$chr}{$allele}{$split[0]} = length($split[2]);
			}
		}
		close FH;
	}
	print STDERR "Save in tree\n";
	foreach my $c (keys %save) {
		foreach my $allele (keys %{$save{$c}}) {
			my $tree = Set::IntervalTree->new;
			foreach my $pos (sort {$a <=> $b} keys %{$save{$c}{$allele}}) {
				$tree->insert( {'mut'=>1 }, $pos, $pos + $save{$c}{$allele}{$pos});
			}
			$tree_mut{$c}{$allele} = $tree;
		}
	}
	return (\%tree_mut, \%last);
}


sub wait_10_secs{
	my $file = $_[0];
	print STDERR "$file exists already!\n";
	print STDERR "Waiting for 10 seconds before overwriting!\n";
	print STDERR "Press Ctrl + C to abort\n";
	print STDERR "Wait";
	for(my $i = 0; $i < 10; $i++) {
		print STDERR ".";
		sleep(1);
	}
	print STDERR "\n";
}
1;

