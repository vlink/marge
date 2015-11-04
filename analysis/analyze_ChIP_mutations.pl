#!/usr/bin/perl


use strict;
use Getopt::Long;
#require '/Users/verenalink/workspace/strains/general/config.pm';
require '../general/config.pm';
use Data::Dumper;
#require '../db_part/database_interaction.pm';

$_ = "" for my($genome, $file, $tf, $filename, $candidates);
$_ = () for my(@strains, %peaks, @split, %mutation_pos, %shift, %current_pos, %save_local_shift, %seq, %PWM, @fileHandlesMotif, %index_motifs);
$_ = 0 for my($homo, $allele, $region, $motif_score);
my $data = config::read_config()->{'data_folder'};
#$data = "/Users/verenalink/workspace/strains/data/";

my %comp;
$comp{'A'} = 'T';
$comp{'C'} = 'G';
$comp{'G'} = 'C';
$comp{'T'} = 'A';

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
		"-region=s" => \$region)
        or die("Error in command line options!\n");
#First step: Get the sequences for the peaks
$allele = 1;
my $ref_save = 0;

#Save motif files
&read_motifs($tf);

#Make sure the reference is in the strains array
my $reference = config::reference_genome($genome);
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
	if($strains[$i] eq $reference) {
		$ref_save = 1;
	}
}
if($ref_save == 0) {
	push @strains, $reference;
}


print STDERR "Saving peaks\n";
open FH, "<$file";
my $line_number = 0;
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#") {
		next;
	}
	@split = split('\t', $line);
	if(length($split[1]) > 10) { next; }
	$peaks{substr($split[1], 3)}{$split[2]} = $split[3];
	$line_number++;
}
close FH;


my $tmp_out = "tmp" . rand(15);
print STDERR "Extracting sequences from strain genomes\n";
print STDERR $tmp_out . "\n";
&get_seq_for_peaks($tmp_out);

my $tmp_out_main_motif = "tmp" . rand(15);

#Now scan for motifs
print STDERR "Scanning for chipped motif\n";
my $command = "homer2 find -i " . $tmp_out . " -m " . $tf . " > " . $tmp_out_main_motif; 
print STDERR $command . "\n";
`$command`;

print STDERR "Read in the motif file!\n";
for(my $files = 0; $files < keys %index_motifs; $files++) {
	$fileHandlesMotif[$files]->print("motif\tpos\t");
	for(my $i = 0; $i < @strains; $i++) {
		$fileHandlesMotif[$files]->print("tg " . $strains[$i] . "\t");
	}
	for(my $i = 0; $i < @strains; $i++) {
		$fileHandlesMotif[$files]->print("motif_score " . $strains[$i] . "\t");
	} 
	for(my $i = 0; $i < @strains; $i++) {
		$fileHandlesMotif[$files]->print("exist " . $strains[$i] . "\t");
	}
	$fileHandlesMotif[$files]->print("#motifs\n");
}
&analyze_motifs();
&print_results();

for(my $i = 0; $i < @fileHandlesMotif; $i++) {
	close $fileHandlesMotif[$i];
}

my $last_line = "";
my @header;
my $motif_start = 0;
my %block;
my %analysis_result = ();
my $more_motifs = 0;
my $save_pos = 0;
my %existance;
my %diff;


sub analyze_motifs{
	open FH, "<$tmp_out_main_motif";
	open FH, "<$tmp_out_main_motif";
	$last_line = "";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@header = split('_', $split[0]);
		#Read in blocks - stop after gathering all information about one position
		if($last_line ne "" && $header[0] . "_" . $header[1] . "_" . $header[2] ne $last_line) {
			&print_results();
		}
		$last_line = $header[0] . "_" . $header[1] . "_" . $header[2];
		my $half = length($seq{$split[0]})/2;
		$half = int($half);
		$motif_start = $header[1] + $half;
		#First get start pos of the motif
		if($split[4] eq "-" && $split[1] < 0) {
		#	$motif_start = $motif_start - ((length($split[2])/2)) + $split[1];	
			$motif_start = $motif_start - (length($split[2])) + $split[1] + 1;	
		}
		if($split[4] eq "-" && $split[1] >= 0) {
		#	$motif_start = $motif_start - ((length($split[2])/2)) + $split[1];
			$motif_start = $motif_start - (length($split[2])) + $split[1] + 1;
		}
		if($split[4] eq "+") {
			$motif_start = $motif_start + $split[1];
		}
		$motif_start = $motif_start - $header[1];
	#	print $line . "\n";
	#	print "without shifting: " . $motif_start . "\n";
	#	print "shifting vector: " . $save_local_shift{$split[0]}{$motif_start} . "\n";
	#	foreach my $key (sort {$a <=> $b} keys %{$save_local_shift{$split[0]}}) {
	#		print $key . "\t" . $save_local_shift{$split[0]}{$key} . "\n";
	#	}
		$block{$split[3]}{$motif_start + $save_local_shift{$split[0]}{$motif_start}}{$header[-1]} = $split[-1];
		$block{$split[3]}{$motif_start + $save_local_shift{$split[0]}{$motif_start}}{'length'} = length($split[2]);
	}
}

sub print_results {
	%analysis_result = ();
	#Now run through all motifs and see if they are in all the other strains
#	print "\n\n\n";
#	print Dumper %block;
	my %ignore = ();
	foreach my $motif (keys %block) {
		%existance = ();
		%diff = ();
		foreach my $motif_pos (keys %{$block{$motif}}) {
			if(exists $ignore{$motif_pos}) { 
				delete $block{$motif}{$motif_pos};
				next;
			}
			$save_pos = $motif_pos;
			$motif_score = 0;
			for(my $i = 0; $i < @strains; $i++) {
				if(!exists $diff{$strains[$i]}) {
					$diff{$strains[$i]} = 0;
				}
				if(!exists $block{$motif}{$motif_pos}{$strains[$i]}) {
					#Check in vicinity (+/- motif_lenght-1) - if there is a motif within motif/2 -> save the current motif under this position, to not count them as separate motifs
					$existance{$strains[$i]} = 1;
				#	print "Check motif pos: " . $motif_pos . "\n";
					for(my $run_motif = $motif_pos - ($block{$motif}{$motif_pos}{'length'} - 1); $run_motif < $motif_pos + ($block{$motif}{$motif_pos}{'length'} - 1); $run_motif++) {
				#		print "run_motif: " . $run_motif. "\n";
						if(exists $block{$motif}{$run_motif} && exists $block{$motif}{$run_motif}{$strains[$i]} && $run_motif != $motif_pos) {
						
							$block{$motif}{$motif_pos}{$strains[$i]} = $block{$motif}{$run_motif}{$strains[$i]};
							delete $existance{$strains[$i]};
							delete $block{$motif}{$run_motif}{$strains[$i]};
							$ignore{$run_motif} = 1;
							if(exists $block{$motif}{$run_motif}{'length'}) {
								delete $block{$motif}{$run_motif}{'length'};
							}
						}
					}
				} 
				$diff{$strains[$i]} += $block{$motif}{$motif_pos}{$strains[$i]};
			}
		}
		$fileHandlesMotif[$index_motifs{$motif}]->print($motif . "\t" . $last_line . "\t");
	#	print "\n";
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$index_motifs{$motif}]->print("tg " . $strains[$i] . "\t");
		}
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$index_motifs{$motif}]->print($diff{$strains[$i]} . "\t");
		}
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$index_motifs{$motif}]->print($strains[$i]. " ");
			if(exists $existance{$strains[$i]}) {
				$fileHandlesMotif[$index_motifs{$motif}]->print("1\t");
			} else {
				$fileHandlesMotif[$index_motifs{$motif}]->print("0\t");
			}
		}
		$fileHandlesMotif[$index_motifs{$motif}]->print("" . (keys %{$block{$motif}}) . "\n");
#		print Dumper %block;
#		print "\n\n";
	}
	%block = ();
}

sub rev_comp{
	my $local_seq = $_[0];
	my $comp_seq = "";
	my @local_split = split("", $local_seq);
	for(my $i = 0; $i < @local_split; $i++) {
		$comp_seq .= $comp{$local_split[$i]};
	}	
	return reverse($comp_seq);
}


sub calculate_motif_score{
	my $score_forward = 0;
	my $score_reverse = 0;
	my $local_seq = $_[1];
	my @base = split('', $local_seq);
	for(my $b = 0; $b < @base; $b++) {
		$score_forward += log($PWM{$_[0]}{$b}{$base[$b]}/0.25); 
	}

	$local_seq = &rev_comp($local_seq);	
	@base = split('', $local_seq);
	for(my $b = 0; $b < @base; $b++) {
		$score_reverse += log($PWM{$_[0]}{$b}{$base[$b]}/0.25); 
	}
	if($score_reverse < 0) { $score_reverse = 0;}
	if($score_forward < 0) { $score_forward = 0;}
	if($score_reverse > $score_forward) {
		return sprintf("%.6f", $score_reverse);	
	} else {
		return sprintf("%.6f", $score_forward);
	}
}


sub read_motifs{
	open FH, "<$_[0]";
	my $pos = 0;
	my $motif = "";
	my $index = 0;
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			$motif = (split('\t', $line))[1];
			$pos = 0;
			#Open file per motif
			my @short = split /[\(,\/]/, $motif;
			my $filename = "output_mutation_" . $short[0] . ".txt";
			open my $fh, ">", $filename or die "Can't open $filename: $!\n";
			$index_motifs{$motif} = $index;
			$fileHandlesMotif[$index] = $fh;
		#	print Dumper %index_motifs;
		#	print "\n\n";
		#	print Dumper @fileHandlesMotif;
			$index++;
		} else {
			@split = split('\t', $line);
			$PWM{$motif}{$pos}{'A'} = $split[0];
			$PWM{$motif}{$pos}{'C'} = $split[1];
			$PWM{$motif}{$pos}{'G'} = $split[2];
			$PWM{$motif}{$pos}{'T'} = $split[3];
			$pos++;
		}
	}
}

sub get_seq_for_peaks {
	open OUT, ">$_[0]";
	my $general_offset = 0;
	my $seq_number = 0;
	my $seq = "";
	my $length;
	my $newlines = 0;
	my $newlines_seq = 0;
	my @fileHandles;
	my $shift_start = 0;
	my $shift_stop = 0;
	my $strain_spec_start;
	my $strain_spec_stop;
	my $mut_number_start = 0;
	my $mut_number_stop = 0;

	foreach my $p (keys %peaks) {
		#Start with offset for first line
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
		my $working_start;
		my $working_stop;
		my $byte_offset;
		#File needs to be sorted, because the shift vector is saved that way
		foreach my $start (sort {$a <=> $b} keys %{$peaks{$p}}) {
			for(my $i = 0; $i < @strains; $i++) {
				$mut_number_start = $mut_number_stop = 0;
				$length = $peaks{$p}{$start} - $start;
				#Check the offset for the strains
				#We run through a sorted list, so we can just remember the last offset per strains so we save time
				$shift_start = 0;
				$shift_stop = 0;
				$working_start = $start - int($region/2);
				$working_stop = $peaks{$p}{$start} + int($region/2);
				#set strain spec start and stop the working start and stop (including the subtraction/addition of the region) for strains without any shiftings
				$strain_spec_start = $working_start;
				$strain_spec_stop = $working_stop;
				#Start at the last position we used for shifting in the previous peak
			#	print "current pos: " . $current_pos{$strains[$i]} . "\n";
			#	print "working start: " . $working_start . "\n";
				for(my $run = $current_pos{$strains[$i]}; $run < @{$shift{$strains[$i]}}; $run++) {
				#	print "run: " . $run . "\n";
				#	print "strain: " . $strains[$i] . "\tmutation pos: " . $mutation_pos{$strains[$i]}[$run] . "\n";
					#Position of our peak is greater than the last mutation, so we just use the last shifting value
					if($working_start > $mutation_pos{$strains[$i]}[-1]) {
						$strain_spec_start = $working_start + $shift{$strains[$i]}[-1];
						$strain_spec_stop = $working_stop + $shift{$strains[$i]}[-1];
						$shift_start = 1;
						$shift_stop = 1;
					}
					#No shifting of position yet, and mutation pos greater than our current pos - so we use one position earlier and use this shifting vector
					if($shift_start == 0 && $mutation_pos{$strains[$i]}[$run] > $working_start) {
					#	print "mut pos: " . $mutation_pos{$strains[$i]}[$run] . "\tworking start: " . $working_start . "\n";
						if($run == 0) {
							$current_pos{$strains[$i]} = 0;
							$strain_spec_start = $working_start;
							$mut_number_start = $run;
						} else {
						#	print "here!\n";
							$current_pos{$strains[$i]} = $run;
							$strain_spec_start = $working_start + $shift{$strains[$i]}[$run-1];
							$mut_number_start = $run - 1;
						}
						$shift_start = 1;
					#	print "mut number start: " . $mut_number_start . "\n";
					#	print "mut pos: " . $mutation_pos{$strains[$i]}[$run] . "\n";
					}
					#No shifting of position yet, and mutation pos greater than our current pos - so we use one position earlier and use this shifting vector
					if($shift_stop == 0 && $mutation_pos{$strains[$i]}[$run] > $working_stop) {
						if($run == 0) {
							$strain_spec_stop = $working_stop;
							$mut_number_stop = 0;
						} else {
							$strain_spec_stop = $working_stop + $shift{$strains[$i]}[$run - 1];
							$mut_number_stop = $run - 1;
						}
						$shift_stop = 1;
					}				
				}
				my $header = "chr" . $p . "_" . $start . "_" . $peaks{$p}{$start} . "_" . $strains[$i];
				print OUT ">chr" . $p . "_" . $start . "_" . $peaks{$p}{$start} . "_" . $strains[$i] . "\n";
				$length = $strain_spec_stop - $strain_spec_start;
				$newlines = int($strain_spec_stop/50) - int($strain_spec_start/50);
				$length = $length + $newlines;
				$newlines = int(($strain_spec_start)/50);
				$byte_offset = $strain_spec_start + $general_offset + $newlines;
				seek($fileHandles[$i], $byte_offset, 0);
				read $fileHandles[$i], $seq, $length;
				$seq =~ s/\n//g;
				print OUT uc($seq) . "\n";
				$seq{$header} = uc($seq);
				$seq = "";
				my $local_shift = 0;
				my $insert = 0;
			#	print $header . "\n";
				if($mut_number_start == $mut_number_stop) {
					for(my $j = 0; $j < ($strain_spec_stop - $strain_spec_start); $j++) {
						$save_local_shift{$header}{$j} = 0; 
					}
				} else {
					my $tmp_mut = $mut_number_start + 1;
					my $before_vector;
					if(defined $mutation_pos{$strains[$i]}[$tmp_mut]) {
						$before_vector = $mutation_pos{$strains[$i]}[$tmp_mut];
					} else {
						$before_vector = $mutation_pos{$strains[$i]}[-1];
					}
					my $last_shift = 0;
					my $index_shift = 0;
					for(my $j = 0; $j < ($strain_spec_stop - $strain_spec_start); $j++) {
						if($working_start + $j + $index_shift == $before_vector) {
						#	print "go in before vector\n";
						#	print "current pos: " . ($working_start + $j + $index_shift) . "\n";
						#	print "before vector: " . $before_vector . "\n";
							$local_shift = ($shift{$strains[$i]}[$tmp_mut - 1] - $shift{$strains[$i]}[$tmp_mut]);
							$last_shift = $last_shift + $local_shift;
							$tmp_mut++;
							$before_vector = $mutation_pos{$strains[$i]}[$tmp_mut];
						#	print "new before vector: " . $before_vector . "\n";
							$insert = 0;
						}
						if($last_shift >= 0 && $insert == 0) {
							for(my $k = 0; $k <= $last_shift + 1; $k++) {
							#	print "k: " . $k . "\tj: " . $j . "\n";
							#	print "index shift: " . $index_shift . "\n";
								$save_local_shift{$header}{$j - $k} = $last_shift;
							#	print "header: " . $header . "\tpos: " . ($j - $k) . "\tshift: " . $last_shift . "\n";
							}
							$insert = 1;
						} else {
							if($local_shift < 0) {
							#	print "relative position: " . ($j + $index_shift) . "\n";
								for(my $k = 0; $k > $last_shift; $k--) {
							#		print "local shift k < 0: " . $k . "\n";
							#		print "before vector: " . $mutation_pos{$strains[$i]}[$tmp_mut - 1] . "\tand current pos: " . ($working_start + $j + $index_shift) . "\n";
							#		print $shift{$strains[$i]}[$tmp_mut] . "\n";
							#		print "save local shift{" . $header . "}{" . ($j + $index_shift - $k) . "} = " . $last_shift . "\n";
									$save_local_shift{$header}{$j + $index_shift - $k} = $last_shift;
							#		print "header: " . $header . "\tpos: " . ($j + $index_shift - $k) . "\tshift: " . ($last_shift - $local_shift) . "\n";
								}
								$save_local_shift{$header}{$j + $index_shift - $last_shift} = $last_shift;
							#	print "after save local shift{" . $header . "}{" . ($j + $index_shift - $last_shift) . "} = " . ($last_shift - $local_shift) . "\n";
							#	$index_shift = $index_shift + ($last_shift * (-1) - 1);
								$index_shift = 0;
								$local_shift = 0;
							} else {
								$save_local_shift{$header}{$j + $index_shift} = $last_shift;
							#	print "else header: " . $header . "\tpos: " . ($j + $index_shift) . "\tshift: " . $last_shift . "\n";
							}
						}
					#	print "j: " . $j . "\n";
					#	print "index shift: " . $index_shift . "\n";	
					#	print "new final position: " . ($working_start + $j + $index_shift) . "\n";
					#	print "before vector: " . $before_vector . "\n";
					}
				}
			}
			$seq_number++;
			print STDERR "" . int(($seq_number/$line_number) * 100) . "% of sequences are gathered\r";
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
	my @tmp;
	my $i = 0;
	if(-e $filename) {
		close FH;
		open FH, "<$filename";
		my $run = 0;
		my @a;
		$i = 0;
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
