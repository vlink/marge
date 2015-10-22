#!/usr/bin/perl -w

use strict;
use Getopt::Long;
require '../general/config.pm';
use Data::Dumper;
require '../db_part/database_interaction.pm';

$_ = "" for my($genome, $file, $tf, $filename, $candidates);
$_ = () for my(@strains, %peaks, @split, %mutation_pos, %shift, %current_pos);
$_ = 0 for my($homo, $allele, $region);
my $data = config::read_config()->{'data_folder'};
sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - if not defined only the reference is used - if all is specified all strains are used\n";
        print STDERR "\t-file <file>: Peaks\n";
	print STDERR "\t-TF <transcription factor motif matrix>: Matrix of the TF that was chipped for\n";
	print STDERR "\t-candidates <TF considered as collaborative candidates>\n";
        print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-region: Size of the region used to look for other motifs (Default: 200)\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $param = config::read_config();

$region = 200;

GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
                "strains=s{,}" => \@strains,
                "homo" => \$homo, 
		"-TF=s" => \$tf, 
		"-candidates=s" => \$candidates,
		"-region" => \$region)
        or die("Error in command line options!\n");

#First step: Get the sequences for the peaks
$allele = 1;

open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#") {
		next;
	}
	@split = split('\t', $line);
	$peaks{substr($split[1], 3)}{$split[2]} = $split[3];
}
close FH;
my $tmp_out = "tmp" . rand(15);
&get_seq_for_peaks($tmp_out);
print $tmp_out . "\n";
my $tmp_out_main_motif = "tmp" . rand(15);

#Now scan for motifs
print STDERR "Scanning for chipped motif\n";
my $command = "homer2 find -i " . $tmp_out . " -m " . $tf . " > " . $tmp_out_main_motif; 
print $command . "\n";
`$command`;

sub get_seq_for_peaks {
	open OUT, ">$_[0]";
	my $general_offset = 0;
	my $seq = "";
	my $length;
	my $newlines = 0;
	my $newlines_seq = 0;
	my $start_pos;
	my @fileHandles;
	my $shift_start = 0;
	my $shift_stop = 0;
	my $strain_spec_start;
	my $strain_spec_stop;
	foreach my $p (keys %peaks) {
		#Start with offset for firest line
		$general_offset = 5 + length($p);
		for(my $i = 0; $i < @strains; $i++) {
			#Read in the offsets for the different strains
			$strains[$i] =~ s/,//g;
			$strains[$i] = uc($strains[$i]);
			$current_pos{$strains[$i]} = 0;
			&create_offset($p, $allele, uc($strains[$i]));
			#Open the chromosome sequences of the different strains
			$filename = $data . "/" . uc($strains[$i]) . "/chr" . $p . ".fa";
			open my $fh, "<", $filename or die "Can't open $filename: $!\n";
			$fileHandles[$i] = $fh;
		}
		foreach my $start (sort {$a <=> $b} keys %{$peaks{$p}}) {
			$length = $region + $peaks{$p}{$start} - $start;
			$start_pos = $start + $general_offset + $newlines - int($region/2);
			for(my $i = 0; $i < @strains; $i++) {
				#Check the offset for the strains
				#We run through a sorted list, so we can just remember the last offset per strains so we save time
				$shift_start = 0;
				$shift_stop = 0;
				$strain_spec_start = $start;
				$strain_spec_stop = $start + $length;
				#Start at the last position we used for shifting in the previous peak
				for(my $run = $current_pos{$strains[$i]}; $run < @{$shift{$strains[$i]}}; $run++) {
					#Position of our peak is greater than the last mutation, so we just use the last shifting value
					if($start > $mutation_pos{$strains[$i]}[-1]) {
						$strain_spec_start = $start + $shift{$strains[$i]}[-1];
						$strain_spec_stop = $peaks{$p}{$start} + $shift{$strains[$i]}[-1];
						$shift_start = 1;
						$shift_stop = 1;
					}
					#No shifting of position yet, and mutation pos greater than our current pos - so we use one position earlier and use this shifting vector
					if($shift_start == 0 && $mutation_pos{$strains[$i]}[$run] > $start) {
						$current_pos{$strains[$i]} = $run - 1;
						$shift_start = 1;
						$strain_spec_start = $start + $shift{$strains[$i]}[$run-1];
					}
					#No shifting of position yet, and mutation pos greater than our current pos - so we use one position earlier and use this shifting vector
					if($shift_stop == 0 && $mutation_pos{$strains[$i]}[$run] > $peaks{$p}{$start}) {
						$strain_spec_stop = $peaks{$p}{$start} + $shift{$strains[$i]}[$run - 1] + 1;
						$shift_stop = 1;
						last;
					}				
				}
				print OUT ">chr" . $p . "_" . $start . "_" . $peaks{$p}{$start} . "_" . $strains[$i] . "\n";
				$length = $strain_spec_stop - $strain_spec_start + $region;
				$newlines = int($strain_spec_stop/50) - int($strain_spec_start/50);
				$length = $length + $newlines;
				$newlines = int($strain_spec_start/50);
				$strain_spec_start = $strain_spec_start + $general_offset + $newlines - int($region/2);
				seek($fileHandles[$i], $strain_spec_start, 0);
				read $fileHandles[$i], $seq, $length;
				$seq =~ s/\n//g;
				print OUT uc($seq) . "\n";
			}
		}
	}
	close OUT;
}

sub create_offset{
	my $chr = $_[0];
	my $a = $_[1];
	my $strain = $_[2];
	$filename = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".shift";
	my @array_mut = ();
	my @array_shift = ();
	if(-e $filename) {
		open FH, "<$filename";
		my $run = 0;
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$array_mut[$run] = $split[0];
			$array_shift[$run] = $split[1] - $split[0];
			$run++;
		}
	}
	$shift{$strain} = [@array_shift];
	$mutation_pos{$strain} = [@array_mut];
}
