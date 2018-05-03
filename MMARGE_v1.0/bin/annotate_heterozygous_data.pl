#!/usr/bin/env perl
use strict;
use warnings;

# Copyright Verena M. Link <vlink@ucsd.edu>
# 
# This file is part of MMARGE
#
# MMARGE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MMARGE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


use Getopt::Long;
use Storable;
use config;
use general;
use analysis_tree;
use Set::IntervalTree;
use Data::Dumper;


#Check if files are empty before starting to calcualte ratio
#
#
$_ = "" for my($output, $command, $data, $chr, $output_resize, $merged_file);
$_ = () for my(@tag_perfect, @tag_mut, @ind, @input, @peak_files, @name, @split, $tree_tmp, $mut, $tree, %lookup_strain, %last_strain, %remove, %mut_in_peak, @muts, @ann_perfect, @ann_mut, %total_reads, %all_reads, $fh, @allele, @align_ind, %mut_peak, @tag_perfect_name, @tag_mut_name);
$_ = 0 for my($hetero, $bp, $resize, $allele, $mut_exists, $min, $max, $center, $first, $no_resize, $no_ann, $F1, $bed_format);

sub printCMD {
	print STDERR "\n\nUsage:\n";
	print STDERR "\t-tag_perfect <list of tag directories from perfectly aligned reads>\n";
	print STDERR "\t-tag_mut <list of tag directories from reads overlapping mutations only>\n";
	print STDERR "\t-ind <list of individuals>\n";
	print STDERR "\t-hetero: individuals are heterozygous\n";
	print STDERR "\t-bp <basepairs>: read lengths (default: 50bp)\n";
	print STDERR "\t-input_perfect <list of tag directores from input - either input per individual, or input per allele>\n";
	print STDERR "\t-peak_files <list of peak files>: Input peak files - if not provided script calls peaks with HOMER\n";
	print STDERR "\t-ann_perfect <list of peak files annotated with perfect mapped reads>\n";
	print STDERR "\t\tFiles have to be annotated with HOMER's annotatePeaks or be in the exact same format. Bed format does not work\n";
	print STDERR "\t-ann_mut <list of peak files annotated with mut mapped reads>\n";
	print STDERR "\t\tFiles have to be annotated with HOMER's annotatePeaks or be in the exact same format. Bed format does not work!\n";
	print STDERR "\t-merged_file <merged peak file>\n";
	print STDERR "\t\tFile can be HOMER's peak file format or bed format. If file has bed format, it is converted into HOMER's peak format.\n";
	print STDERR "\t-no_resize: does not resize peak file\n";
	print STDERR "\t-no_ann: do not annotate individual peak files\n";
	print STDERR "\t-resize: only resizes peaks without any further annotation\n";
	print STDERR "\t-output_resize <output file name>: output file for resized coordinates (default: merged_peaks_resized.txt)\n"; 
	print STDERR "\t-output <output file>: default annotated_heterozygous_peaks.txt\n";
	print STDERR "\t-data_dir <folder to data directory>: default defined in config\n"; 
	print STDERR "\t-F1: Annotate allele-specific based on a F1 mouse strain. In this case, please specify only 2 genomes - both parental mouse strains. These genomes are then used as different alleles. This does only work for two genomes, not more!\n";
	print STDERR "\n\n\n";
	exit;
}

if(@ARGV < 1) { 
	&printCMD();
}

#Check mandatory command line arguments
my %mandatory = ('-tag_perfect' => 1, '-tag_mut' => 1, '-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(     "tag_perfect=s{,}" => \@tag_perfect,
		"tag_mut=s{,}" => \@tag_mut, 
		"ind=s{,}" => \@ind,
		"hetero" => \$hetero,
		"bg=s" => \$bp,
		"input_perfect=s{,}" => \@input,
		"peak_files=s{,}" => \@peak_files,
		"ann_perfect=s{,}" => \@ann_perfect,
		"ann_mut=s{,}" => \@ann_mut,
		"merged_file=s{,}" => \$merged_file,
		"no_resize" => \$no_resize,
		"resize" => \$resize,
		"no_ann" => \$no_ann,
		"output_resize=s" => \$output_resize,
		"data_dir=s" => \$data,
		"F1" => \$F1,
		"output=s" => \$output)
        or die("Error in command line options!\n");




if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}

if($hetero == 0) {
	$allele = 1;
} else { 
	$allele = 2;
}

if(@tag_perfect == 1) {
	@tag_perfect = split(",", $tag_perfect[0]);
}
if(@tag_mut == 1) {
	@tag_mut = split(",", $tag_mut[0]);
}
if(@ind == 1) {
	@ind = split(",", $ind[0]);
}
for(my $i = 0; $i < @ind; $i++) {
	$ind[$i] =~ s/,//g;
	$ind[$i] = uc($ind[$i]);
}
#Check for length
if(@tag_perfect != @tag_mut) { 
	print STDERR "ERROR: Different number of perfect tag directories and mutated tag directores\n\n\n";
	&printCMD();
}
if(@input > 0 && (@input != @tag_perfect && @input != (@tag_perfect/2))) {
	print STDERR "ERROR: Number of input directories does not overlap with perfect tag directories\n\n\n";
	&printCMD();
}
if(@tag_perfect != @ind && $hetero == 0) {
	print STDERR "ERROR: Number of individuals and tag directories does not match!\n\n\n";
	&printCMD();
}

for(my $i = 0; $i < @tag_perfect; $i++) {
	$tag_perfect[$i] =~ s/,//g;
	$tag_mut[$i] =~ s/,//g;
	@split = split("/", $tag_perfect[$i]);
	if(@split > 1) {
		$tag_perfect_name[$i] = $split[-1];
	} else {
		$tag_perfect_name[$i] = $tag_perfect[$i];
	}
	@split = split("/", $tag_mut[$i]);
	if(@split > 1) {
		$tag_mut_name[$i] = $split[-1];
	} else {
		$tag_mut_name[$i] = $tag_mut[$i];
	}
}

for(my $i = 0; $i < @ind; $i++) {
	$ind[$i] =~ s/,//g;
	my($tree_ref, $test, $last) = general::read_strains_mut($ind[$i], $data);
	$tree->{$ind[$i]} = $tree_ref;
}

if(@ann_mut == 1) {
	@ann_mut = split(",", $ann_mut[0]);
}
if(@ann_mut > 0) {
	for(my $i = 0; $i < @ann_mut; $i++) {
		$ann_mut[$i] =~ s/,//g;
	}
}

if(@ann_perfect == 1) {
	@ann_perfect = split(",", $ann_perfect[0]);
}
if(@ann_perfect > 0) {
	for(my $i = 0; $i < @ann_perfect; $i++) {
		$ann_perfect[$i] =~ s/,//g;
	}
}

if(@input == 1) {
	@input = split(",", $input[0]);
}
if(@peak_files == 1) {
	@peak_files = split(",", $peak_files[0]);
}

if($F1 == 1) {
	$allele = 1;
	if(@ind != 2) {
		print STDERR "F1 parameter set: Please specify exactly two parental genomes!\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if(@tag_perfect != 2) {
		print STDERR "F1 parameter set: Please specify exactly two perfectly mapped tag directories\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if(@tag_mut != 2) {
		print STDERR "F1 parameter set: Please specify excatly two tag directories with mutations only!\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if(@input > 0 && @input != 2) {
		print STDERR "F1 parameter set: Please specify exactly two input tag directories!\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if(@ann_perfect > 0 && @ann_perfect != 2) {
		print STDERR "F1 parameter set: Please specify exactly two files with annotations from perfectly mapped tag directories\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if(@ann_mut > 0 && @ann_mut != 2) {
		print STDERR "F1 parameter set: Please specify exactly two files with annotations from tag directories from mutations only!\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if(@peak_files > 1 && @peak_files != 2) {
		print STDERR "F1 parameter set: Please specify excatly two peak files!\n";
		print STDERR "\n\n";
		&printCMD();
	}
	if($hetero == 1) {
		print STDERR "F1 parameter set: Change parameter to homozygous!\n";
		$hetero = 0;
	}
}


if($output eq "") { $output = "annotated_heterozygous_peaks.txt"; }
if($output_resize eq "") { 
	if($merged_file eq "") {
		$output_resize = "merged_peaks_resized.txt"; 
	} else {
		$output_resize = $merged_file . "_resized.txt";
	}
}


if($bp == 0) { $bp = 50; }


#Remove last / from directories
for(my $i = 0; $i < @tag_perfect; $i++) {
	$tag_perfect[$i] =~ s/,//g;
	$tag_mut[$i] =~ s/,//g;
	while(substr($tag_perfect[$i], length($tag_perfect[$i]) - 1) eq "/") {
		chop $tag_perfect[$i];
	}
	while(substr($tag_mut[$i], length($tag_mut[$i]) - 1) eq "/") {
		chop $tag_mut[$i];
	}
	if($i < @input) {
		$input[$i] =~ s/,//g;
		while(substr($input[$i], length($input[$i]) - 1) eq "/") {
			chop $input[$i];
		}
	}
	if($i < @peak_files) {
		$peak_files[$i] =~ s/,//g;
	}
}


#First call peaks if peak files are not inserted
if(@peak_files == 0 && $merged_file eq "") {
	print STDERR "\tCalling peaks for all different tag directories\n";
	if(@input == 0) {
		for(my $i = 0; $i < @tag_perfect; $i++) {
			print STDERR "\t" . $tag_perfect[$i] . "\n";
			@name = split("/", $tag_perfect[$i]);
			$command = "findPeaks " . $tag_perfect[$i] . " > peak_" . $name[-1] . ".txt 2> tmp";
			#print $command . "\n";			
			`$command`;
			$peak_files[$i] = "peak_" . $name[-1] . ".txt";
		}
	} elsif(@input == @tag_perfect) {
		for(my $i = 0; $i < @tag_perfect; $i++) {
			print STDERR "\t" . $tag_perfect[$i] . "\n";
			@name = split("/", $tag_perfect[$i]);
			$command = "findPeaks " . $tag_perfect[$i] . " -i " . $input[$i] . " > peak_" . $name[-1] . ".txt 2>> tmp";
			#print $command . "\n";
			`$command`;
			$peak_files[$i] = "peak_" . $name[-1] . ".txt";
		}
	} else {
		my $input_count = 0;
		for(my $i = 0; $i < @tag_perfect; $i++) {
			print STDERR "\t" . $tag_perfect[$i] . "\n";
			@name = split("/", $tag_perfect[$i]);
			$command = "findPeaks " . $tag_perfect[$i] . " -i " . $input[$input_count] . " > peak_" . $name[-1] . ".txt 2>> tmp";
			if($i > 0 && $i % 2 == 0) {
				$input_count++;
			}
			#print $command . "\n";
			`$command`;
			$peak_files[$i] = "peak_" . $name[-1] . ".txt";
		}
	}
}

if($merged_file eq "") {
	#Merge peak file
	$merged_file = "merged_peaks.txt";
	$command = "mergePeaks "; 
	for(my $i = 0; $i < @peak_files; $i++) {
		$command .= $peak_files[$i] . " ";
	}
	$command .= " > merged_peaks.txt 2>> tmp";
	`$command`;
} else {
	#Check if merged file is bed format
	$first = 0;
	open my $fh, "<", $merged_file;
	while(my $line = <$fh>) {
		if($first == 0) { $first++; next; }
		@split = split('\t', $line);
		if(@split == 3) {
			$bed_format = 1;
		}
	}
	close $fh;
}
$first = 0;

my $tmp = "tmp_" . rand(20);
if($bed_format == 1) {
	print STDERR "Merged input file is bed format - convert to HOMER peak format\n";
	$command = "annotatePeaks.pl " . $merged_file . " mm10 -noann -nogene > " . $tmp;
	`$command`;
	$command = "mv " . $tmp . " " . $merged_file;
	`$command`;
}

#Read in all peaks
open FH, "<", $merged_file;
$merged_file = $output_resize;
open OUT, ">", $output_resize;
#Resize peaks to
#	no mutation: center +/- bp
#	one mutation: mut +/- bp
#	several muts: center muts furtherst apart +/- bp
#This needs to be discussed
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#") { next; }
	@split = split('\t', $line);
	$chr = substr($split[1], 3);
	$mut_exists = 0;
	@muts = ();
	for(my $i = 0; $i < @ind; $i++) {
		for(my $a = 1; $a <= $allele; $a++) {
			if(exists $tree->{$ind[$i]}->{$chr}) {
				$tree_tmp = $tree->{$ind[$i]}->{$chr}->{$a}->fetch($split[2], $split[3]);
				if(exists $tree_tmp->[0]->{'mut'}) {
					$mut_peak{$split[0]}{$ind[$i]} = 1;
					$mut = $tree_tmp->[0]->{'mut'};
					for(my $j = 0; $j < @{$tree_tmp}; $j++) {
						push @muts, $mut;
					}
					$mut_exists = 1;
				}
			}
		}
	}
	if($no_resize == 0) {
		#print STDERR "Resize\n";
	#	if(exists $remove{$split[0]}) { next; }
		if($mut_exists == 0) {
			#print STDERR "mot exists\n";
			$center = $split[2] + ($split[3] - $split[2])/2;
			print OUT $split[0] . "\t" . $split[1] . "\t" . int($center - $bp) . "\t" . int($center + $bp) . "\t+\n"; 
		} else {
			if(@muts == 1) {
				print OUT $split[0] . "\t" . $split[1] . "\t" . int($muts[0] - $bp) . "\t" . int($muts[0] + $bp) . "\t+\n";
			} else {
				$min = $muts[0];
				$max = $muts[0];
				for(my $i = 1; $i < @muts; $i++) {
					if($min > $muts[$i]) { $min = $muts[$i]; }
					if($max < $muts[$i]) { $max = $muts[$i]; }
				}
				$center = $min + ($max - $min)/2;
				print OUT $split[0] . "\t" . $split[1] . "\t" . int($center - $bp) . "\t" . int($center + $bp) . "\t+\n";
			}
		}
	} else {
		print OUT $line . "\n";
	}
}
close OUT;

if($resize == 1) {
	print STDERR "File is successfully resized\n";
	print STDERR "Resized coordinates are written to file: " . $output_resize . "\n";
	exit;
}

if($no_ann == 0) {
	#Annotate peak file
	if(@ann_perfect == 0) {
		for(my $i = 0; $i < @tag_perfect; $i++) {
			$command = "annotatePeaks.pl " . $merged_file . " none -strand + -d " . $tag_perfect[$i] . " -noadj > merged_peaks_ann_perfect_" . $tag_perfect_name[$i] . ".txt";
			print STDERR $command . "\n";
			`$command`;
		}
	}
	if(@ann_mut == 0) {
		for(my $i = 0; $i < @tag_mut; $i++) {
			$command = "annotatePeaks.pl " . $merged_file . " none -strand + -d " . $tag_mut[$i] . " -noadj > merged_peaks_ann_muts_" . $tag_mut_name[$i] . ".txt";
			print STDERR $command . "\n";
			`$command`;
		}
	}
}

#map individuals to tag directories
if(@tag_perfect == @ind) { 
	@allele = (1)x@ind;
       	@align_ind = @ind;
} else {
	for(my $i = 0; $i < @ind; $i++) {
		push @align_ind, $ind[$i];
		push @align_ind, $ind[$i];
		push @allele, 1;
		push @allele, 2;
	}
}



#Read in files and normalize
#	That needs to be discuessed
for(my $i = 0; $i < @tag_perfect; $i++) {
	print STDERR "open merged_peaks_ann_perfect_" . $tag_perfect_name[$i] . ".txt\n";
	open $fh, "<", "merged_peaks_ann_perfect_" . $tag_perfect_name[$i] . ".txt";
	$first = 0;
	while(my $line = <$fh>) {
		chomp $line;
		@split = split('\t', $line);
		if($first == 0) {
		#	@name = split('\(', $split[19]);
		#      	@split = split(' ', $name[1]);
		#	$total_reads{$align_ind[$i]}{$allele[$i]}{'p'} = $split[0];
			$first++;
		} else {
			$all_reads{$split[0]}{$align_ind[$i]}{$allele[$i]}{'p'} = $split[19];
			$total_reads{$align_ind[$i]}{$allele[$i]}{'p'} += $split[19];

		}
	}
	close $fh;
}

for(my $i = 0; $i < @tag_mut ;$i++) {
	print STDERR "open merged_peaks_ann_muts_" . $tag_mut_name[$i] . ".txt\n";
	open $fh, "<", "merged_peaks_ann_muts_" . $tag_mut_name[$i] . ".txt";
	$first = 0;
	while(my $line = <$fh>) {
		chomp $line;
		@split = split('\t', $line);
		if($first == 0) {
		#	@name = split('\(', $split[19]);
		#	@split = split(' ', $name[1]);
		#	$total_reads{$align_ind[$i]}{$allele[$i]}{'m'} = $split[0];
			$first++;
		} else {
			$all_reads{$split[0]}{$align_ind[$i]}{$allele[$i]}{'m'} = $split[19];
			$total_reads{$align_ind[$i]}{$allele[$i]}{'m'} += $split[19];
		}
	}
	close $fh;
}

my %ratio;
my $complete_reads = 0;
#Now annotate the peaks
open $fh, "<", $merged_file;
print STDERR "Opening $merged_file\n";
open my $out, ">", $output;
print STDERR "Output all data into " . $output . "\n";
$first = 0;
my $max_perfect = 0;

for(my $i = 0; $i < @ind; $i++) {
	for(my $a = 1; $a <= $allele; $a++) {
		if($max_perfect == 0) {
			$max_perfect = $total_reads{$align_ind[$i]}{$allele[$i]}{'p'};
		}
		if($total_reads{$align_ind[$i]}{$allele[$i]}{'p'} > $max_perfect) {
			$max_perfect = $total_reads{$align_ind[$i]}{$allele[$i]}{'p'};
		}
	}
}
my $value;
foreach my $line (<$fh>) {
	if($first == 0) { $first++; next; }
	chomp $line;
	@split = split('\t', $line);
	print $out $line;
	if(substr($line, 0, 1) eq "#") { print $out "\n"; next; }
	if($F1 == 0) {
		for(my $i = 0; $i < @ind; $i++) {
			if(exists $mut_peak{$split[0]}{$ind[$i]}) {
				$complete_reads = 0;
				%ratio = ();
				for(my $a = 1; $a <= $allele; $a++) {
					$complete_reads += $all_reads{$split[0]}{$ind[$i]}{$a}{'m'};
				}
				if($complete_reads == 0) { $complete_reads = 1; }
				for(my $a = 1; $a <= $allele; $a++) {
					$ratio{$a} = ($all_reads{$split[0]}{$ind[$i]}{$a}{'m'})/$complete_reads;
				}
			}
			for(my $a = 1; $a <= $allele; $a++) {
				#if(!exists $ratio{$a}) { print STDERR "Could not find any ratio!\n"; print STDERR "Make sure that mutation directory is correct!\n"; exit;}
				if(exists $ratio{$a} && $ratio{$a} == 0) { $ratio{$a} = 1; }
				if(exists $mut_peak{$split[0]}{$ind[$i]}) {
					print $out "\t" . ((($all_reads{$split[0]}{$ind[$i]}{$a}{'p'}/$total_reads{$ind[$i]}{$a}{'p'})/2) * $ratio{$a}) * ($max_perfect/$total_reads{$ind[$i]}{$a}{'p'}) * 1000000;
				} else {
					$value = (($all_reads{$split[0]}{$ind[$i]}{$a}{'p'}/$total_reads{$ind[$i]}{$a}{'p'})/2) * 1000000 * ($max_perfect/$total_reads{$ind[$i]}{$a}{'p'});
					print $out "\t" . ($value + int(rand($value * 0.10)));
				}	
			}
		}
		print $out "\n";
	} else {
		if(exists $mut_peak{$split[0]} == 1) {
			$complete_reads = 0;
			%ratio = ();
			for(my $i = 0; $i < @ind; $i++) {
				$complete_reads += $all_reads{$split[0]}{$ind[$i]}{'1'}{'m'};
			}
			if($complete_reads == 0) { $complete_reads = 1; }
			for(my $i = 0; $i < @ind; $i++) {
				$ratio{$i} = ($all_reads{$split[0]}{$ind[$i]}{'1'}{'m'})/$complete_reads;
			}
		} else {
			for(my $i = 0; $i < @ind; $i++) {
				$ratio{$i} = 0;
			}
		}
		for(my $i = 0; $i < @ind; $i++) {
			#if(!exists $ratio{$i}) { print STDERR "Could not find any ratio!\n"; print STDERR "Make sure that mutation directory is correct2!\n"; exit;}
			if(exists $ratio{$i} && $ratio{$i} == 0) {
				print $out "\t" . 0;
			} else {
				if(exists $mut_peak{$split[0]}) {
					print $out "\t" . ((($all_reads{$split[0]}{$ind[$i]}{'1'}{'p'}/$total_reads{$ind[$i]}{'1'}{'p'})/2) * $ratio{$i}) * 1000000;
				} else {
					print $out "\t" . (($all_reads{$split[0]}{$ind[$i]}{'1'}{'p'}/$total_reads{$ind[$i]}{'1'}{'p'})/2) * 1000000;
				}
			}
		}
		print $out "\n";
	}
}
close $fh;
close $out;
