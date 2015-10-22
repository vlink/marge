#!/usr/bin/perl -w

use strict;

use FileHandle;
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
	print STDERR "General commands:\n";
	print STDERR "\t-file <input>\n";
	print STDERR "Statistics:\n";
	print STDERR "\t-stats: outputs only the statistics on STDOUT (default)\n";
	print STDERR "\t-files <prefix>: outputs the different groups into different files startin with prefix\n";
	exit;
}

$_ = "" for my($file, $fh, $files, $tmp_name);
$_ = () for my(%header, @split, @name, %mut_exists, %exact, %strain_overlap, %mut_no_overlap, %mut_in_others, %no_mut_at_all, %peak_in_strain, @filehandle, %save_filehandle_number, %all_combis);
$_ = 0 for my($add_strain, $right_strain, $number_of_muts, $number_of_strains, $number_of_overlap, $mut_overlap_with_peak, $stats);
my %mandatory = ('-file' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "file=s" => \$file, 
		"stats" => \$stats, 
		"files=s" => \$files)
        or die ("Error in command line arguments!\n");

#Start with only unique peaks!
open FH, "<$file";

print STDERR "So far just tested for unique peaks\n";

if($files eq "") {
	$stats = 1;
}

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if(substr($line, 0, 1) eq "#") { 
		for(my $i = 8; $i < @split; $i++) {
			@name = split('\s+', $split[$i]);
			$header{$i} = $name[0];
		}
		next; 
	}
	#Create all the output files`
	if($files ne "") {
		if(!exists $save_filehandle_number{$files . "_exact_" . $split[6]}) {
			$tmp_name = $files . "_exact_" . $split[6];
			open my $fh, ">", $tmp_name or die "Can't open $tmp_name: $!";
			push @filehandle, $fh;
			$save_filehandle_number{$tmp_name} = @filehandle - 1;
		
			$tmp_name = $files . "_add_muts_" . $split[6];
			open my $fh4, ">", $tmp_name or die "Can't open $tmp_name: $!";
			push @filehandle, $fh4;
			$save_filehandle_number{$tmp_name} = @filehandle - 1;

	
			$tmp_name = $files . "_no_mut_in_peak_" . $split[6];
			open my $fh1, ">", $tmp_name or die "Can't open $tmp_name: $!";
			push @filehandle, $fh1;
			$save_filehandle_number{$tmp_name} = @filehandle - 1;

			$tmp_name = $files . "_mut_in_others_" . $split[6];
			open my $fh2, ">", $tmp_name or die "Can't open $tmp_name: $!";
			push @filehandle, $fh2;
			$save_filehandle_number{$tmp_name} = @filehandle - 1;	

			$tmp_name = $files . "_no_mut_at_all_" . $split[6];
			open my $fh3, ">", $tmp_name or die "Can't open $tmp_name: $!";
			push @filehandle, $fh3;
			$save_filehandle_number{$tmp_name} = @filehandle - 1;
		}					
	} else {
		if(!exists $all_combis{$split[6]}) {
			$all_combis{$split[6]} = 1;
		}
	}
	%mut_exists = ();
	%peak_in_strain = ();
	$mut_overlap_with_peak = $number_of_muts = $number_of_strains = $number_of_overlap = 0;
	for(my $i = 8; $i < @split; $i++) {
		if(length($split[$i]) > 2) {
			$mut_exists{$header{$i}} = 1;
			$number_of_muts++;
		}
	}
	foreach my $s (keys %header) {
		if(index($split[6], $header{$s}) != -1) {
			$number_of_strains++;
			#We want to know if all strains with this peak also have a mutation
			$peak_in_strain{$header{$s}} = 1;
		}
		if(exists $mut_exists{$header{$s}}) {
			$number_of_overlap++;
			if(exists $peak_in_strain{$header{$s}}) {
				$mut_overlap_with_peak++;
			}
		}
	}

#	print $line . "\n";
#	print Dumper @split;
#	print "\nmut exists\n";
#	print Dumper %mut_exists;
#	print "\npeak in strain:\n";
#	print Dumper %peak_in_strain;
#	print "all important numbers:\n";
#	print "number of overlap: " . $number_of_overlap . "\n";
#	print "number of strains: " . $number_of_strains . "\n";
#	print "mut overlap with peak: " . $mut_overlap_with_peak . "\n";
#	print "number of muts: " . $number_of_muts . "\n";
#	print "keys peak_in_strain: " . (keys %peak_in_strain) . "\n";

	if($mut_overlap_with_peak == keys %peak_in_strain) {
		#The mutations in the strains overlap with the peak
		#Now figure out if there are any other mutations in any other strains
		if($number_of_overlap == $mut_overlap_with_peak) {
			#Exactly the number of mutations
		#	print "DECISION: exact\n";
			$exact{$split[6]}++;
			if($files ne "") {
				$filehandle[$save_filehandle_number{$files . "_exact_" . $split[6]}]->print($line . "\n");
			}
		} else {
			#There are additional mutations
			$strain_overlap{$split[6]}++;
		#	print "DECISION: add_muts\n";
			if($files ne "") {
				$filehandle[$save_filehandle_number{$files . "_add_muts_" . $split[6]}]->print($line . "\n");
			}
		}
	} elsif($mut_overlap_with_peak > 0 && $mut_overlap_with_peak < keys %peak_in_strain) {
		print STDERR "That is the part for the not unique merging!\n";
		print STDERR "Still todo\n";
		exit;
	} else {
		#The mutations do not overlap with the peak
		#Are there any mutations at all?
		if($number_of_muts == 0) {
			$no_mut_at_all{$split[6]}++;
		#	print "DECISION: no mut at all\n";
			if($files ne "") {
				$filehandle[$save_filehandle_number{$files . "_no_mut_at_all_" . $split[6]}]->print($line . "\n");
			}
		}
		#All strains except for the one with the peaks are mutated
		if($number_of_muts == (keys %header) - $number_of_strains) {
			$mut_in_others{$split[6]}++;
		#	print "DECISION: mut in others and just not in the strain with the peak\n";
			if($files ne "") {
				$filehandle[$save_filehandle_number{$files . "_mut_in_others_" . $split[6]}]->print($line . "\n");
			}
		}
		#Some strains have mutations
		if($number_of_muts > 0) {
			$mut_no_overlap{$split[6]}++;
		#	print "DECISION: some strains have muts\n";
			if($files ne "") {
				$filehandle[$save_filehandle_number{$files . "_no_mut_in_peak_" . $split[6]}]->print($line . "\n");

			}
		}
	}
}

if($files ne "") {
	foreach my $number (keys %save_filehandle_number) {
		close $filehandle[$save_filehandle_number{$number}];
	}
} else {
	print "strains\t#exact mut\t#additional mut\t#mut only in others\t#mut in all others\t#no mut at all\n";
	foreach my $s (keys %all_combis) {
		print $s . "\t";
		if(exists $exact{$s}) {
			print $exact{$s} . "\t";
		} else {
			print "0\t";
		}
		if(exists $strain_overlap{$s}) {
			print $strain_overlap{$s} . "\t";
		} else {
			print "0\t";
		}
		if(exists $mut_no_overlap{$s}) {
			print $mut_no_overlap{$s} . "\t";
		} else {
			print "0\t";
		}
		if(exists $mut_in_others{$s}) {
			print $mut_in_others{$s} . "\t";
		} else {
			print "0\t";
		}
		if(exists $no_mut_at_all{$s}) {
			print $no_mut_at_all{$s} . "\n"; 
		} else {
			print "0\n";
		}
	}
}
#Format for insert file
#ID	chr	start	stop	strain	Stat	unique to which strain	Total subpeaks	strain1	strain2	strain3	etc
