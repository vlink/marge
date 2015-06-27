#!/usr/bin/perl -w

#Merge peaks with different parameters
#Use mergePeaks 
#	just keep unique ones
#	keep all combinations
#Use tag counts and foldchange
#Use getDiffPeaks

use strict;
use Getopt::Long;
#use DBI;
require '../general/config.pm';
require '../general/system_interaction.pm';
use Data::Dumper;
use File::Basename;
use Cwd;

$|=1;

if(@ARGV < 1) {
	&printCmd();
}

sub printCmd {
        print STDERR "Usage:\n";
	print STDERR "\n";
	print STDERR "If you want to use follow up programs, please make sure that the filename contains the strain in the correct spelling!\n";
	print STDERR "\n\n";
        print STDERR "General commands:\n";
        print STDERR "\t-method <merge|foldchange|diffPeaks>\n";
	print STDERR "\t-files <list of files, comma separated>\n";
	print STDERR "\t-out <output file>\n";
	print STDERR "\nMethod merge:\n";
	print STDERR "\t-unique: <Just keeps unique peaks in output file> (default)\n";
	print STDERR "\t-all: <Keeps all peaks in output file>\n";
	print STDERR "\nMethod foldchange:\n";
	print STDERR "\tInput one annotated peak file\n";
	print STDERR "\t-order <order of strains, comma separated>\n";
	print STDERR "\t-fc <foldchange>: default: 2\n";
	print STDERR "\t-filter <threshold>: filters out all peaks where every strain has less than threshold peaks (recommended!) (default: 0)\n";
	print STDERR "\t-diff: Will only output peaks that are different between at least 2 strains\n";
	print STDERR "\nMethod diffPeaks\n";
	print STDERR "\t-dirs <list of tag directories>: comma-seperated list of tag directores (same order as -order)\n";
	print STDERR "\t-dirlist <file with tag directores>\n";
	print STDERR "\t\tFormat: <path to tag directory>,strain\n";
	print STDERR "\t-order <order of strains, comma separated>\n";
	print STDERR "\t-diff: Will only output peaks that are different between at least 2 strains\n";
	print STDERR "\nRest TODO\n";
        exit;
}

$_ = "" for my($method, $output, $command, $file_list, $ref_pos, $up, $down, $dirlist);
$_ = 0 for my($unique, $all, $is_installed, $first, $fc, $filter, $over_filter, $diff);
$_ = () for my (@files, %del, @split, @order, @dirs, %tagdirs, %save_diff_results_up, %save_diff_results_down);

my %mandatory = ('-method' => 1, '-output' => 1, '-files' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "method=s" => \$method,
                "output=s" => \$output,
		"files=s{,}" => \@files,
		"diff" => \$diff,
                "unique" => \$unique,
                "all" => \$all, 
		"order=s{,}" => \@order,
		"fc=s" => \$fc, 
		"filter=s" => \$filter,
		"dirs=s{,}" => \@dirs,
		"dirlist=s" => \$dirlist)
        or die ("Error in command line arguments!\n");

my $tmp = "100000";

if($method eq "merge") {
	$is_installed = system_interaction::check_if_homer_is_installed();
	if($is_installed == 1) {
		print STDERR "Please install homer or choose another merging method!\n";
		exit;
	}
	$command = "rm " . $tmp . ".tmp";
	`$command`;
	my $pwd = cwd();
	$pwd .= "/";
	my ($filename, $path_to_file) = fileparse($files[0]);
	if($pwd ne $path_to_file) {
		print STDERR "copy files to tmp file\n";
		for(my $i = 0; $i < @files; $i++) {
			($filename, $path_to_file) = fileparse($files[$i]);
			$command = "cp " . $files[$i] . " .";
			`$command`;
			$files[$i] = $filename;
			$del{$filename} = 1;
		}
	}
	for(my $i = 0; $i < @files; $i++) {
		$file_list .= $files[$i] . " ";
	}

	if($unique == 1) {
		print STDERR "merging peaks!\n";
		$command = "mergePeaks " . $file_list . " -prefix t 2> /dev/null";
		$del{"t_*"} = 1;
		`$command`;
		for(my $i = 0; $i < @files; $i++) {
			$command = "cat t_" . $files[$i] . " >> " . $tmp . ".tmp";
			`$command`; 
		}
	} else {
		print STDERR "merging peaks!\n";
		$command = "mergePeaks " . $file_list . " > " . $tmp . ".tmp 2> /dev/null";
		`$command`;
	}
	$del{$tmp . ".tmp"} = 1;
	$command = "cat " . $tmp . ".tmp | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7}' > final_" . $tmp . ".tmp";
	`$command`;
	$command = "mv final_" . $tmp . ".tmp " . $output;
	`$command`;
} elsif($method eq "foldchange") {
	open OUT, ">$output";
	#Check file format first
	#Needs to be a peak file with annotation
	if($fc == 0) {
		$fc = 2;
	}
	if(@order == 0) {
		print STDERR "Please define order of strains!\n";
		exit;
	}
	for(my $i = 0; $i < @order; $i++) {
		$order[$i] =~ s/,//g;
		if($order[$i] eq "reference") {
			$ref_pos = (@split - @order) + $i;
		}
	}
	if($ref_pos eq "") {
		print STDERR "Please specify the position of the reference in order!\n";
		exit;
	}

	open FH, "<$files[0]";
	if(@files > 1) {
		print STDERR "Only first file is considered! Ignore the rest!\n";
	}
	#First line check format
	foreach my $line (<FH>) {
		if($first == 0) {
			if(substr($line, 0, 6) ne "PeakID") {
				print STDERR "Please check file format!\n";
				print STDERR "Does not look like an original peak file!\n";
			}
			@split = split('\t', $line);
			if(@split < 6) {
				print STDERR "Please check file format!\n";
				print STDERR "Is this peak file annotated?\n";
			}
			$first++;
			print OUT "#id\tchr\tstart\tend\tstrand\tStat\tParent files\n";
			next;
		}
		#Now let's check the pattern
		chomp $line;
		@split = split('\t', $line);
		$up = "up_";
		$down = "down_";
		$over_filter = 0;
		print $line . "\n";
		for(my $i = @split - @order; $i < @split; $i++) {
			if($i == $ref_pos) {
				if($split[$ref_pos] > $filter) {
					$over_filter = 1;
				}
				next;
			}
			if($split[$i] > $filter) {
				$over_filter = 1;
			}
			if(($split[$i] + 1)/($split[$ref_pos]+1) > $fc) {
				$up .= $order[$i - (@split - @order)] . "_";
			}
			if(($split[$i] + 1)/($split[$ref_pos]+1) < 1/$fc) {
				$down .= $order[$i - (@split - @order)] . "_";
			}
		}
		if($over_filter == 0) {
			next;
		}
		chop $up;
		chop $down;
		if($diff == 1 && $up eq "" && $down eq "") {
			next;
		}
		print OUT $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . $split[3] . "\t" . $split[4] . "\t" . $split[5] . "\t";
		if(length($up) > 3) {
			print OUT $up;
		}
		if(length($down) > 5) {
			if(length($up) > 3) {
				print OUT "_" . $down;
			} else {
				print OUT $down;
			}
		}
		print OUT "\n";
	}
} elsif($method eq "diffPeaks") {
	$is_installed = system_interaction::check_if_homer_is_installed();
	if($is_installed == 1) {
		print STDERR "Please install homer or choose another merging method!\n";
		exit;
	}

	#Check if all parameters are given
	if(@order == 0) {
		print STDERR "Please define order of strains!\n";
		exit;
	}
	for(my $i = 0; $i < @order; $i++) {
		$order[$i] =~ s/,//g;
		if($order[$i] eq "reference") {
			$ref_pos = $i;
		}
	}
	if($ref_pos eq "") {
		print STDERR "Please specify the position of the reference in order!\n";
		exit;
	}
	#tag directores
	if(@dirs == 0 && $dirlist eq "") {
		print STDERR "Please specify either a comma separated list of tag directires (-dirs) or a file containing tag directories and strains (-dirlist)\n";
		exit;
	}
	#Check dirlist
	if($dirlist ne "") {
		open FH, "<$dirlist";
		foreach my $line (<FH>) {
			chomp $line;
			@split = split(",", $line);
			$tagdirs{$split[1]} = $split[0];
		}
		close FH;
	}
	if(@dirs > 0) {
		#Check if we have as many tag dirs as we have strains
		if(@dirs != @order) {
			print STDERR "The number of inserted strains and tag directories does not overlap!\n";
			print STDERR "Please check that every strain has a corresponding tag directory!\n";
			exit;
		} else {
			for(my $i = 0; $i < @dirs; $i++) {
				$dirs[$i] =~ s/,//g;
				$tagdirs{$order[$i]} = $dirs[$i];
			}
		}
	}
	#Start getting differential peaks, then open the file and read it in
	print STDERR "Running getDifferentialPeaks for all strains versus reference!\n";
	for(my $i = 0; $i < @order; $i++) {
		if($i == $ref_pos) {
			next;
		}
		print "compare " . $order[$ref_pos] . " vs " . $order[$i] . "\n";
		$command =  "getDifferentialPeaks " . $files[0] . " " . $tagdirs{$order[$ref_pos]} . " " . $tagdirs{$order[$i]} . " > compare_" . $order[$ref_pos] . "_" . $order[$i] . "_spec_" . $order[$ref_pos] . " 2> /dev/null";
	#	print $command . "\n";
		`$command`;
		$del{"compare_" . $order[$ref_pos] . "_" . $order[$i] . "_spec_" . $order[$ref_pos]} = 1;
		open FH, "<compare_" . $order[$ref_pos] . "_" . $order[$i] . "_spec_" . $order[$ref_pos];
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			if(substr($line, 0, 1) eq "#") {
				next;
			}
			@split = split('\t', $line);
			$save_diff_results_down{$split[0]} .= $order[$i] . "_";
		}
		close FH;
		$command =  "getDifferentialPeaks " . $files[0] . " " . $tagdirs{$order[$ref_pos]} . " " . $tagdirs{$order[$i]} . " -rev > compare_" . $order[$ref_pos] . "_" . $order[$i] . "_spec_" . $order[$i] . " 2> /dev/null";
	#	print $command . "\n\n\n";
		$del{"compare_" . $order[$ref_pos] . "_" . $order[$i] . "_spec_" . $order[$i]} = 1;
		`$command`;
		open FH, "<compare_" . $order[$ref_pos] . "_" . $order[$i] . "_spec_" . $order[$i];
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			if(substr($line, 0, 1) eq "#") {
				next;
			}
			$save_diff_results_up{$split[0]} .= $order[$i] . "_";
		}
		close FH;
	}
	#Write output file
	open OUT, ">$output";
	open FH, "<$files[0]";
	print OUT "#id\tchr\tstart\tend\tstrand\tStat\tParent files\n";
	foreach my $line (<FH>) {
		if(substr($line, 0, 1) eq "#") {
			next;
		}
		chomp $line;
		@split = split('\t', $line);
		if($diff == 1 && !exists $save_diff_results_up{$split[0]} && !exists $save_diff_results_down{$split[0]}) {
			next;
		}
		print OUT $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . $split[3] . "\t" . $split[4] . "\t" . $split[5] . "\t";
		if(exists $save_diff_results_up{$split[0]}) {
			chop $save_diff_results_up{$split[0]};
			print OUT "up_" . $save_diff_results_up{$split[0]};
		}
		if(exists $save_diff_results_down{$split[0]}) {
			chop $save_diff_results_down{$split[0]};
			if(exists $save_diff_results_up{$split[0]}) {
				print OUT "_";
			}
			print OUT "down_" . $save_diff_results_down{$split[0]};
		}
		print OUT "\n";
	}	
	close FH;
	close OUT;
} else {
	print STDERR "Unknown method!\n";
	&printCmd();
}

foreach my $f (keys %del) {
	`rm $f`;
}
