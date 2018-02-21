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
use Set::IntervalTree;
use config;
use general;
use analysis_tree;
use Data::Dumper;

$_ = "" for my($output, $data_dir, $genome_dir, $strain, $chr, $method, $output_file, $print_line);
$_ = () for my(%peaks, %strand, @split, %tree, %tree_pos, %last_strain, @tmp_split, %save_id, $tree_pos, $tree, @files, $tree_tmp, @strains, %last, $last, @chr_split, $tree_tmp_1, $tree_tmp_2, %tree_detail, $tree_detail, $tree_tmp_detail_1, $tree_tmp_detail_2, @detail);
$_ = 0 for my($hetero, $line_number, $id, $verbose);

sub printCMD {
        print STDERR "\n\nUsage:\n";
        print STDERR "\t-ind <individual>: Individual we look for muts versus reference\n";
	print STDERR "\t-inds <individuals>: Two individuals we look for muts against each other\n";
        print STDERR "\t-files <files>: Comma separated list of files\n";
	print STDERR "\t-method <bowtie|star>: define mapping method (default: bowtie)\n";
	print STDERR "\t-hetero: Data is heterozygous (Default: homozygous)\n";
	print STDERR "\nAdditional parameters:\n";
	print STDERR "\t-v: verbose mode - progress monitoring\n";
	print STDERR "\t-data_dir <directory>: default defined in config\n";
	print STDERR "\t-genome_dir <directory>: default defined in config\n\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my %mandatory = ('-files' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(   	"files=s{,}" => \@files,
		"ind=s" => \$strain,
		"inds=s{,}" => \@strains,
		"method=s" => \$method,
		"genome_dir=s" => \$genome_dir,
		"data_dir=s" => \$data_dir,
		"hetero" => \$hetero, 
		"v" => $verbose)
	or die(&printCMD());
#First step: Get the sequences for the peaks

#Set variables
if($data_dir eq "") {
	$data_dir = config::read_config()->{'data_folder'};
}
if($genome_dir eq "") {
	$genome_dir = $data_dir;
}
if($strain eq "" && @strains < 1) {
	&printCMD();
}
if($method ne "bowtie" && $method ne "star") {
	$method = "bowtie";
}
if($strain ne "") {
	$strain = uc($strain);
}

if(@strains == 1) {
	my @a = split(",", $strains[0]);
	for(my $i = 0; $i < @a; $i++) {
		$a[$i] =~ s/,//g;
		$strains[$i] = uc($a[$i]);
	}
} elsif(@strains > 1) {
	for(my $i = 0; $i < @strains; $i++) {
		$strains[$i] =~ s/,//g;
		$strains[$i] = uc($strains[$i]);
	}
}

for(my $i = 0; $i < @files; $i++) {
	$files[$i] =~ s/,//g;
}
print STDERR "Read in mutations\n";
if(@strains < 1) {
#	($tree, $tree_detail, $last) = general::read_strains_mut($strain, $data_dir);
	($tree, $tree_pos, $last) = general::read_strains_mut($strain, $data_dir);
} else {
	($tree, $last) = general::read_mutations_from_two_strains($strains[0], $strains[1], $data_dir);
}

my $next = 0;
my $no_mut = 0;
my $all_lines = 0;
my $printed_lines = 0;
print STDERR "Processing sam file\n";
my @name;
my %save;
my @base;
my $complete_lines = 0;

foreach my $file (@files) {
	print STDERR $file . "\n";
	if($verbose == 1) {
		$complete_lines = (split('\s', `wc -l $file`))[0];
	}
	open my $fh, "<", $file;
#	open FH, "<$file";
	@split = split("/", $file);
	$output = "";
	for(my $i = 0; $i < @split - 1; $i++) {
		$output .= $split[$i] . "/";
	}
	$output_file = $output .  "only_muts_" . $split[-1];
	open my $mut_reads, ">", $output_file;
	$output_file = $output . "perfect_reads_" . $split[-1];
	open my $perfect, ">", $output_file;
	$_ = 0 for($next, $no_mut, $all_lines, $printed_lines);
	while(my $line = <$fh>) {
		$all_lines++;
		if($verbose == 1) {
			print STDERR "Processing " . $file . ": " . $all_lines . " of " . $complete_lines . "(" . sprintf("%.2f", ($all_lines/$complete_lines)) . "%)" . "\r"; 
		}
		chomp $line;
		if(substr($line, 0, 1) eq "@") {
			if(substr($line, 0, 3) eq "\@SQ") {
				@split = split('\t', $line);
				@detail = split(":", $split[1]);
				@chr_split = split("_", $detail[1]);
				#	if(@chr_split > 1) {
				#	print STDERR "Your input file was not shifted!\n";
				#	print STDERR "Abort\n";
				#	exit;
				#}
				print $mut_reads $split[0] . "\tSN:" . $chr_split[0] . "\t" . $split[2] . "\n";
				print $perfect $split[0] . "\tSN:" . $chr_split[0] . "\t" . $split[2] . "\n";
			} else {
				print $mut_reads $line . "\n";
				print $perfect $line . "\n";
			}	
		} else {
			@split = split('\t', $line);
			if($split[2] eq "*") { next; }
			@chr_split = split("_", $split[2]);
			if(@chr_split > 1) {
				$split[2] = $chr_split[0];
			}
			$chr = substr($split[2], 3);
			if(!exists $tree->{$chr}) { next; }
			$print_line = $split[0] . "\t" . $split[1] . "\t" . $split[2];
			for(my $i = 3; $i < @split; $i++) {
				$print_line .= "\t" . $split[$i];
			}
			if($method eq "bowtie") {
				 if($line =~ m/XM:i:0/) {
					 print $perfect $print_line . "\n";
				 } else {
					$next++;
					next;
				 }
			} else {
				if($line =~ m/nM:i:0/) {
					print $perfect $print_line . "\n";
				} else {
					$next++;
					next;
				}
			}
			if($hetero == 0) {
				if($split[3] + length($split[9]) > $last->{$chr}->{'1'}->{'pos'}) { $next++; next; }
				$tree_tmp = $tree->{$chr}->{1}->fetch($split[3], $split[3] + length($split[9]));
				if(exists $tree_tmp->[0]->{'mut'}) {
					print $mut_reads $print_line . "\n";
					$printed_lines++;
				} else {
					$no_mut++;
				}
			} else {
				if($split[3] + length($split[9]) > $last->{$chr}->{'1'}->{'pos'} || $split[3] + length($split[9]) > $last->{$chr}->{'2'}->{'pos'}) { $next++; next; }
				$tree_tmp_1 = $tree->{$chr}->{'1'}->fetch($split[3], $split[3] + length($split[9]));
				$tree_tmp_2 = $tree->{$chr}->{'2'}->fetch($split[3], $split[3] + length($split[9]));
				$tree_tmp_detail_1 = $tree_pos->{$chr}->{'1'}->fetch($split[3], $split[3] + length($split[9]));
				$tree_tmp_detail_2 = $tree_pos->{$chr}->{'2'}->fetch($split[3], $split[3] + length($split[9]));
				if(exists $tree_tmp_1->[0]->{'mut'} && !exists $tree_tmp_2->[0]->{'mut'}) {
					print $mut_reads $print_line . "\n";
					$printed_lines++;
				} elsif(!exists $tree_tmp_1->[0]->{'mut'} && exists $tree_tmp_2->[0]->{'mut'}) {
					print $mut_reads $print_line . "\n";
					$printed_lines++;
				} elsif(exists $tree_tmp_1->[0]->{'mut'} && exists $tree_tmp_2->[0]->{'mut'}) {
					if((scalar @{$tree_tmp_1}) != (scalar @{$tree_tmp_2})) {
						print $mut_reads $print_line . "\n";
						$printed_lines++;
					} else {
						for(my $i = 0; $i < (scalar @{$tree_tmp_1}); $i++) {
							if($tree_tmp_1->[$i]->{'mut'} ne $tree_tmp_2->[$i]->{'mut'}) {
								print $mut_reads $print_line . "\n";
								$printed_lines++;
								last;
							} else {
								if($tree_tmp_detail_1->[$i]->{'pos'} ne $tree_tmp_detail_2->[$i]->{'pos'}) {
									print $mut_reads $print_line . "\n";
									$printed_lines++;
									last;
							 	}	
							}
						}
					}
				} else {
					$no_mut++;
				}
			}
		}
	}
	close $fh;
	close $mut_reads;
	close $perfect;
	if($verbose == 1) {
		print STDERR "\n\n";
	}
	open LOG, ">$output.log";
	print LOG "Get only sequences spanning mutations for $file\n";
	print LOG "All lines looked at:\t\t" . $all_lines. "\n";
	print LOG "All lines spanning mutations:\t" . $printed_lines . "\t(" . ($printed_lines/$all_lines) . ")\n";
	print LOG "All lines not spanning mutations:\t" . $no_mut . "\t(" . ($no_mut/$all_lines) . ")\n";
	print LOG "All lines skipped:\t\t" . $next . "\t(" . ($next/$all_lines) . ")\n";
	close LOG;
}


