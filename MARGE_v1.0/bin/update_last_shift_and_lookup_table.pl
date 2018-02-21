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
use Storable;
use config;
use processing;
use Set::IntervalTree;
use Data::Dumper;

$_ = 0 for my($hetero, $num_strains, $diff, $allele);
$_ = () for my(@strains, @split, %last_shift_pos_strain, %last_shift_pos_ref, %lookup_table, @tmp, @split_last, @tmp2);
$_ = "" for my($data, $last_line_strains, $last_line_ref, $command, $out);


sub printCMD{
	print STDERR "\nThis script updates the last_shift and lookup tables\n";
	print STDERR "It should be used when mutation files for each individual were generated using a different file per chromosome\n";
	print STDERR "The -add parameter must have been used\n";
	print STDERR "The last_shift and lookup tables only contain the last entry, so it needs to be updated and all chromosomes need to be added\n\n\n";
	print STDERR "Usage:\n";
	print STDERR "\t-ind <list of individuals>\n";
	print STDERR "\t-hetero: Strain is heterozygous (Default: homozygous)\n";
	print STDERR "\t-data_dir: Directory where mutation files are located - default: folder specified in config file\n\n\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

#Read in command line arguments
GetOptions(	"-ind=s{,}" => \@strains,
		"-hetero" => \$hetero,
		"-data_dir=s" => \$data)
or die(&printCMD());

if($data eq "") {
        $data = config::read_config()->{'data_folder'};
}

$num_strains = @strains;
for(my $i = 0; $i < @strains; $i++) {
        $strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

for(my $i = 0; $i < @strains; $i++) {
	my @files = `ls $data/$strains[$i]/*mut`;

	foreach my $f (@files) {
		print $f;
		chomp $f;
		@tmp2 = split("/", $f);
		@tmp = split("_", $tmp2[-1]);
		$tmp[0] = substr($tmp[0], 3);
		#First step: get the last diff
		open FH, "<$f";
		$allele = substr($f, length($f) - 5, 1);
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			if(length($split[1]) != length($split[2])) {
				$diff = length($split[2]) - length($split[1]);
			}
		}
		#Second step: Get the last line in ref to strain
		$command = "tail -n1 " . substr($f, 0, length($f) - 4) . ".ref_to_strain.vector";
		$last_line_ref = `$command`;
		chomp $last_line_ref;
		@split_last = split('\t', $last_line_ref);
		$last_shift_pos_strain{$tmp[0]}{$allele}{'pos'} = ($split_last[-1]+1);
		$last_shift_pos_strain{$tmp[0]}{$allele}{'shift'} = $split_last[0] + $diff;
		#Third step: Get the last line in strain to ref
		$command = "tail -n1 " . substr($f, 0, length($f) - 4) . ".strain_to_ref.vector";
		$last_line_ref = `$command`;
		chomp $last_line_ref;
		@split_last = split('\t', $last_line_ref);
		$last_shift_pos_ref{$tmp[0]}{$allele}{'pos'} = $split_last[1];
		$last_shift_pos_ref{$tmp[0]}{$allele}{'shift'} = $split_last[0];
	}
	$out = $data . "/" . $strains[$i] . "/last_shift_strain.txt";
        store \%last_shift_pos_ref, "$out";
	$out = $data . "/" . $strains[$i] . "/last_shift_ref.txt";
        store \%last_shift_pos_strain, "$out";
	%last_shift_pos_strain = ();
	%last_shift_pos_ref = ();

}
