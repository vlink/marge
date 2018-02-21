#!/usr/bin/env perl
package analysis;
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



use config;
use POSIX;
#use Statistics::Basic qw(:all);
#use Statistics::Distributions;
#use Statistics::RankCorrelation;;
#use Statistics::TTest;
use Data::Dumper;
use Storable;
#use Statistics::R;
#use Memory::Usage;

$_ = () for my(@split, %PWM, @tmp_split, %comp, $seed);
$comp{'A'} = 'T';
$comp{'T'} = 'A';
$comp{'C'} = 'G';
$comp{'G'} = 'C';
$comp{'-'} = '-';
$comp{'N'} = 'N';

sub write_header{
	my $filehandle = $_[0];
	my @strains = @{$_[1]};
	my $fc = $_[2];
	my $delta = $_[3];
	my $allele = $_[4];
	print $filehandle "motif\tpos\t";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			print $filehandle "tg " . $strains[$i] . " - " . $al . "\t";
		}
	}
	#Check if foldchange is set 
	if($fc == 0) {
		for(my $i = 0; $i < @strains - 1; $i++) {
			for(my $j = $i + 1; $j < @strains; $j++) {
				for(my $a1 = 1; $a1 <= $allele; $a1++) {
					for(my $a2 = 1; $a2 <= $allele; $a2++) {
						if($delta == 1) {
							print $filehandle "delta of tg" . $strains[$i] . "-" . $a1 . " vs " . $strains[$j] . "-" . $a2 . "\t";
						} else {
							print $filehandle "log2 fc ". $strains[$i] . "-" . $a1 . " vs " . $strains[$j] . "-" . $a2 . "\t";
						}
					}
				}
			}
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			print $filehandle "motif_score " . $strains[$i] . "-" . $al . "\t";
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			print $filehandle "exist " . $strains[$i] . "-" . $al . "\t";
		}
	}
	print $filehandle "#motifs\n";
}




#This method summarizes the output of HOMER's motif scan for every peak. It generates a hash that contains all motifs for all strains per peak
#The method reports the motif positions from the beginning of the sequence and takes care of shifting through indels within this sequence
#All motif positions are realtive to the reference motif position
sub analyze_motifs{
	$_ = () for my(@split, @header, %block, $tree_tmp, $last_line);
	$_ = 0 for my($pos_beginning, $shift_beginning, $pos_current, $pos_end, $shift_current, $shift_diff, $shift_end, $motif_pos, $chr_num, $motif_start, $allele_num);
	my $motif_file = $_[0];
	my @strains = @{$_[1]};
	my $tree = $_[2];
	my $lookup = $_[3];
	my $last = $_[4];
	my $allele = $_[5];
	my $region = $_[6];
	my $original_chr_num;
	open(my $fh, "<", $motif_file);
	while(my $line = <$fh>) {
                chomp $line;
                @split = split('\t', $line);
		$split[3] = &get_motif_name($split[3]);
                @header = split('_', $split[0]);
		$last_line = $header[0] . "_" . $header[1] . "_" . $header[2];
		$motif_start = $split[1];
		$motif_pos = $split[1];
		#define motif start
		$motif_start = $motif_start + $header[1];
		#Motif is on negative strand - adjust start position
		if($split[4] eq "-") {
			$motif_start = $motif_start - length($split[2]) + 1;
			$motif_pos = $motif_pos - length($split[2]) + 1;
		}
		$pos_beginning = $header[1] - ($region/2);
		@tmp_split = split("_", $header[0]);
		$chr_num = substr($tmp_split[0], 3);
		$allele_num = $header[4];
		
		if($chr_num !~ /\d+/ && !exists $lookup->{$header[3]}->{$chr_num}) {
			print STDERR "Skip analysis of chromosome " . $chr_num . " allele " . $allele_num . " in analyze_motifs\n";
		}
		$original_chr_num = $chr_num;
		if(exists $lookup->{$header[3]}->{$chr_num}) {
			$chr_num = $lookup->{$header[3]}->{$chr_num};
		}
		#Check if there is a interval tree for this strain and chromosome
		if(defined $tree->{$header[3]}->{$chr_num}->{$allele_num}) {
			#Get shifting vector for the beginning of the sequence that was used to scan for the motif
		#	$tree_tmp = $tree->{$header[3]}->{$chr_num}->fetch($pos_beginning, $pos_beginning + 1);
			#$tree_tmp = $tree->{$header[3]}->{$chr_num}->{$allele_num}->fetch($pos_beginning, $pos_beginning + 1);
			$tree_tmp = $tree->{$header[3]}->{$chr_num}->{$allele_num}->fetch($pos_beginning, $pos_beginning);
			if(!exists $tree_tmp->[0]->{'shift'}) {
				$shift_beginning = $last->{$header[3]}->{$original_chr_num}->{$allele_num}->{'shift'};
			} else {
				$shift_beginning = $tree_tmp->[0]->{'shift'};
			}
			#Get shifting vector for the position of the motif
			#$tree_tmp = $tree->{$header[3]}->{$chr_num}->{$allele_num}->fetch($motif_start, $motif_start + 1);
			$tree_tmp = $tree->{$header[3]}->{$chr_num}->{$allele_num}->fetch($motif_start, $motif_start);
			if(!exists $tree_tmp->[0]->{'shift'}) {
				$shift_current = $last->{$header[3]}->{$original_chr_num}->{$allele_num}->{'shift'};
			} else {
				$shift_current = $tree_tmp->[0]->{'shift'};
			}
			#Calculate the difference between the shfiting vector at the beginning of the sequence and the shifting vector at the beginning of the motif - if different then there is a InDel in this sequence and shifting is necessary
			$shift_diff = $shift_beginning - $shift_current;
			$motif_pos = $motif_pos + $shift_diff;
		}
		#Save motif, length and orientation for this strain at the position in a hash
		$block{$last_line}{$split[3]}{$motif_pos}{$header[-2]}{$header[-1]} = $split[-1];
		$block{$last_line}{$split[3]}{$motif_pos}{'length'} = length($split[2]);
		$block{$last_line}{$split[3]}{$motif_pos}{'orientation'} = $split[4];
        }
	close $fh;
	return \%block;
}

#The block hash is merged when necessary when overlap is set and motifs are overlapping
#The method also determines whether or not there is a mutation in one of the strains
sub merge_block {
	$_ = () for my($tree_tmp, $shift_vector, %ignore, %existance, %diff);
	$_ = 0 for my($num_of_bp, $chr_num, $current_pos, $prev_shift, $current_shift, $save_pos, $motif_score);
	my $block_ref = $_[0];
	my %block = %$block_ref;
	my @strains = @{$_[2]};
	my $overlap = $_[1];
	my $seq = $_[3];
	my $tree = $_[4];
	my $lookup = $_[5];
	my $last = $_[6];
	my $allele = $_[7];
	my $chr_num_original;
        #Now run through all motifs and see if they are in all the other strains
	foreach my $chr_pos (keys %block) {
		$current_pos = (split("_", $chr_pos))[1];
		foreach my $motif (keys %{$block{$chr_pos}}) {
			%existance = ();
			%diff = ();
			#When overlap is set, define the number of basepairs the motif can overlap in order to be merged
			foreach my $motif_pos (keys %{$block{$chr_pos}{$motif}}) {
				if($overlap eq "half") {
					$num_of_bp = int($block{$chr_pos}{$motif}{$motif_pos}{'length'}/2)
				} elsif($overlap eq "complete") {
					$num_of_bp = $block{$chr_pos}{$motif}{$motif_pos}{'length'};
				} elsif($overlap eq "") {
					$num_of_bp = 0;
				} else {
					$num_of_bp = ($overlap*1);
				}
				#Motif was handled already in another part of this method (e.g. motif overlapped and was merged into another peak, so this one has to be deleted)
				if(exists $ignore{$motif_pos}) {
					delete $block{$chr_pos}{$motif}{$motif_pos};
					next;
				}
				$save_pos = $motif_pos;
				$motif_score = 0;
				for(my $i = 0; $i < @strains; $i++) {
					for(my $al = 1; $al <= $allele; $al++) {
						$shift_vector = 0;
						$chr_num = substr((split('_', $chr_pos))[0], 3);
						$chr_num_original = $chr_num;
						if($chr_num !~ /\d+/ && !exists $lookup->{$strains[$i]}->{$chr_num}) {
							print STDERR "Skip analysis of chromosome " . $chr_num . " allele " . $al . " in merge_block\n";
							next;
						}
						if(exists $lookup->{$strains[$i]}->{$chr_num}) {
							$chr_num = $lookup->{$strains[$i]}->{$chr_num};
						}
						if(!exists $diff{$strains[$i]}->{$al}) {
							$diff{$strains[$i]}->{$al} = 0;
						}
						#Motif does not exist at this position in this strain
						if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al}) {
							#Check in vicinity (+/- motif_lenght-1) - if there is a motif within motif/2 -> save the current motif under this position, to not count them as separate motifs
							$existance{$strains[$i]}{$al} = 1;
							#Check if there is an overlap of a motif nearby in this strain (e.g. in this strain there was a indel that shifted the motif by some basepairs)
							#If overlap is not set, $number_of_bp is 0
							for(my $run_motif = $motif_pos - $num_of_bp; $run_motif < $motif_pos + $num_of_bp; $run_motif++) {
								#There is a motif nearby and the current position we are looking it is not the motif position 
								if(exists $block{$chr_pos}{$motif}{$run_motif} && $run_motif != $motif_pos) {
									#Run through all strains now and shift this motif to the current motif position and then delete the other motif position
									for(my $run_strains = 0; $run_strains < @strains; $run_strains++) {
										for(my $run_allele = 1; $run_allele <= $allele; $run_allele++) {
											if(!exists $block{$chr_pos}{$motif}{$run_motif}{$strains[$run_strains]}{$run_allele}) { next; }
											$block{$chr_pos}{$motif}{$motif_pos}{$strains[$run_strains]}{$run_allele} += $block{$chr_pos}{$motif}{$run_motif}{$strains[$run_strains]}{$run_allele};
											$diff{$strains[$run_strains]}{$run_allele} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$run_strains]}{$run_allele};
											delete $existance{$strains[$run_strains]}{$run_allele};
											delete $block{$chr_pos}{$motif}{$run_motif}{$strains[$run_strains]}{$run_allele};
										}
									}
									#Add motif to ignore motif
									$ignore{$run_motif} = 1;
									if(exists $block{$chr_pos}{$motif}{$run_motif}{'length'} || $block{$chr_pos}{$motif}{$run_motif}{'orientation'}) {
										delete $block{$chr_pos}{$motif}{$run_motif};
									}
								}
							}
						}
						#Motif either has no score in this strain or does not exist at all and also was not merged
						if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} || $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} eq "") {
							#Define the position of this motif in the strain that is currently investigated
							#Then pull this sequence from the sequence hash and calculate the motif score
						#	print STDERR "strain: " . $strains[$i] . "\t" . $chr_num . "\t" . $al . "\n";
						#	print STDERR "current pos: " . $current_pos . "\tmotif pos: " . $motif_pos ."\n";
						#	print STDERR "chr original: " . $chr_num_original . "\n";
						#	print STDERR "last: " . $last->{$strains[$i]}->{$chr_num_original}->{$al}->{'pos'} . " vs " . ($current_pos + $motif_pos) . "\n";
						#	print STDERR $seq->{$chr_pos . "_" . $strains[$i] . "_" . $al} . "\n";
							if(defined $tree->{$strains[$i]}->{$chr_num}->{$al}) {
								if($current_pos + $motif_pos >= $last->{$strains[$i]}->{$chr_num_original}->{$al}->{'pos'}) {
						#			print STDERR "Last\n";
									$current_shift = $last->{$strains[$i]}->{$chr_num_original}->{$al}->{'pos'};
								} else {
									#$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$al}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos + 1);
									$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$al}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos);
								#	print STDERR "From sequence\n";
									$current_shift = $tree_tmp->[0]->{'shift'};
								}
								if($current_pos >= $last->{$strains[$i]}->{$chr_num_original}->{$al}->{'pos'}) {
									$prev_shift = $last->{$strains[$i]}->{$chr_num_original}->{$al}->{'pos'};
								} else { 
									#$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$al}->fetch($current_pos, $current_pos + 1);
									$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$al}->fetch($current_pos, $current_pos);
									$prev_shift = $tree_tmp->[0]->{'shift'};
								}
							#	print STDERR "motif: " . $motif . "\t";
							#	print STDERR $current_shift . "\t" . $prev_shift . "\n";
								$shift_vector = $current_shift - $prev_shift;	
							}
						#	print STDERR $seq->{$chr_pos . "_" . $strains[$i] . "_" . $al} . "\n";
						#	print STDERR "motif pos + shift vector: " . ($motif_pos + $shift_vector) . "\n";
						#	print STDERR "length of seq: " . length($seq->{$chr_pos . "_" . $strains[$i] . "_" . $al}) . "\n";
						#	print STDERR "length of motif: " . $block{$chr_pos}{$motif}{$motif_pos}{'length'} . "\n";
							if(length($seq->{$chr_pos . "_" . $strains[$i] . "_" . $al}) < $motif_pos + $shift_vector || length(substr($seq->{$chr_pos . "_" . $strains[$i] . "_" . $al}, $motif_pos + $shift_vector)) < $block{$chr_pos}{$motif}{$motif_pos}{'length'} || $motif_pos + $shift_vector < 0) {
						#		print STDERR "Why are we not setting it to zero?\n";
								$block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} = 0;
							} else {
							#	print STDERR "chr pos: " . $chr_pos .  "\tstrain: " . $strains[$i] . "\tallele: " . $al . "\tmotif pos" . ($motif_pos + $shift_vector) . "\tlength: " . $block{$chr_pos}{$motif}{$motif_pos}{'length'} . "\n";
							#	print STDERR $seq->{$chr_pos . "_" . $strains[$i] . "_" . $al} . "\n";
								if(substr($seq->{$chr_pos . "_" . $strains[$i] . "_" . $al}, $motif_pos + $shift_vector, $block{$chr_pos}{$motif}{$motif_pos}{'length'}) eq "") {
							#		print STDERR "Set to zero\n";
									$block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} = 0;
								} else {
									#	print STDERR substr($seq->{$chr_pos . "_" . $strains[$i] . "_" . $al}, $motif_pos + $shift_vector, $block{$chr_pos}{$motif}{$motif_pos}{'length'}) . "\n";
									$block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} = &calculate_motif_score($motif, substr($seq->{$chr_pos . "_" . $strains[$i] . "_" . $al}, $motif_pos + $shift_vector, $block{$chr_pos}{$motif}{$motif_pos}{'length'}), $block{$chr_pos}{$motif}{$motif_pos}{'orientation'});
								}
							}
							$diff{$strains[$i]}{$al} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al};
						} else {
							$diff{$strains[$i]}{$al} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al};
						}
					}
				}
			}
		}
	}
	return \%block;
}

#Method outputs the analyze motifs for further analysis in other methods 
#Files are per motif and summarize the number of mutation in every peak per strain
sub output_motifs{
	$_ = 0 for my($fc_exists, $max_motif);
	$_ = () for my($curr_chr, $curr_pos, %diff, %existance);
	my $block_ref = $_[0];
	my %block = %$block_ref;
	my $output_file = $_[1];
#	my @fileHandlesMotif = @{$_[1]};
	my %tag_counts = %{$_[2]};
	my @strains = @{$_[3]};
	my %index_motifs = %{$_[4]};
	my %fc = %{$_[5]};
	if((keys %fc) > 1) {
		$fc_exists = 1;
	}
	my %recenter_conversion = %{$_[6]};
	my $motif_diff = $_[7];
	my $motif_diff_percentage = $_[8];
	my $allele = $_[9];
	open OUT, ">$output_file";
        #Now run through all motifs and see if they are in all the other strains
	foreach my $chr_pos (keys %block) {
		@split = split("_", $chr_pos);
		$curr_chr = substr($split[0], 3);
		$curr_pos = $split[1];
		if(exists $recenter_conversion{$curr_chr}{$curr_pos}{'start'}) {
			$curr_pos = $recenter_conversion{$curr_chr}{$curr_pos}{'start'};
		}
		foreach my $motif (keys %{$block{$chr_pos}}) {
			%diff = ();
			%existance = ();
		#	$fileHandlesMotif[$index_motifs{$motif}]->print($motif . "\t" . $chr_pos . "\t");
		#	$fileHandlesMotif[$index_motifs{$motif}]->print($tag_counts{$curr_chr}{$curr_pos});
			print OUT $motif . "\t" . $chr_pos . "\t";
			print OUT $tag_counts{$curr_chr}{$curr_pos};
			foreach my $motif_pos (keys %{$block{$chr_pos}{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					for(my $al = 1; $al <= $allele; $al++) {
						if(!exists $diff{$strains[$i]}{$al}) { $diff{$strains[$i]}{$al} = 0; }
						#Filter by motif existance - either motif is not there at all (default) or diff_motif is set
						if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} || ($motif_diff == 0 && $motif_diff_percentage == 0 && $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} == 0)) {
							$existance{$strains[$i]}{$al}++;
						#	$existance{$strains[$i]} = 1;
						}
						if($block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} > $max_motif) {
							$max_motif = $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al};
						}
					#	$diff{$strains[$i]}{$al} = $diff{$strains[$i]} + $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al};
						$diff{$strains[$i]}{$al} = $diff{$strains[$i]}{$al} + $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al};
					}
				}
				#not binary motif existance, but difference in motif score is used to determine if there is a mutation within the peak
				if($motif_diff != 0) {
					for(my $i = 0; $i < @strains; $i++) {
						for(my $al = 1; $al <= $allele; $al++) {
							if($max_motif - $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} > $motif_diff) {
								$existance{$strains[$i]}{$al}++;
							#	$existance{$strains[$i]} = 1;
							}
						}
					}
				}
				#not binary motif existance or absolute difference in motif score is used, but difference percentage between the two motif scores
				if($motif_diff_percentage != 0) {
					for(my $i = 0; $i < @strains; $i++) {
						for(my $al = 1; $al <= $allele; $al++) {
							if($max_motif - $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}{$al} > ($max_motif * $motif_diff_percentage)) {
							#	$existance{$strains[$i]} = 1;
								$existance{$strains[$i]}{$al}++;
							}
						}
					}
				}
				$max_motif = 0;
			}
			if($fc_exists == 1) {
			#	$fileHandlesMotif[$index_motifs{$motif}]->print($fc{$curr_chr}{$curr_pos});
				print OUT $fc{$curr_chr}{$curr_pos};
			}
			for(my $i = 0; $i < @strains; $i++) {
				for(my $al = 1; $al <= $allele; $al++) {
				#	$fileHandlesMotif[$index_motifs{$motif}]->print($diff{$strains[$i]}{$al} . "\t");
					print OUT $diff{$strains[$i]}{$al} . "\t";
				}
			}
			#output number of mutations in peak for this motif
			for(my $i = 0; $i < @strains; $i++) {
				for(my $al = 1; $al <= $allele; $al++) {
					if(exists $existance{$strains[$i]}{$al}) {
					#	$fileHandlesMotif[$index_motifs{$motif}]->print("1\t");
						print OUT "1\t";
					} else {
					#	$fileHandlesMotif[$index_motifs{$motif}]->print("0\t");
						print OUT "0\t";
					}
				}
			}
		#	$fileHandlesMotif[$index_motifs{$motif}]->print("" . (keys %{$block{$chr_pos}{$motif}}) . "\n");
			print OUT "" . (keys %{$block{$chr_pos}{$motif}}) . "\n";
		}
	}
	close OUT;
}


#Method prepares output hash for distance plots
sub distance_plot{
	my %block = %{$_[0]};
	my @strains = @{$_[1]};
	my %fc = %{$_[2]};
	my $fc_significant = $_[3];
	my $effect = $_[4];
	my $longest_seq = $_[5];
	my $seq = $_[6];
	my $delta_tag = $_[7];
	my $delta_threshold = $_[8];
	my $allele = $_[9];
	my %dist_plot;
	my $offset;
	foreach my $last_line (keys %block) {
		@split = split("_", $last_line);
		#Only keeps peaks where there is an effect, so when the foldchange between the two peaks is not significant peak is ignored - only if fold change of strains, not absoulte difference of tag counts per strain is used
		if($delta_tag == 0 && $effect == 1 && ($fc{substr($split[0], 3)}{$split[1]} < log($fc_significant)/log(2) && $fc{substr($split[0], 3)}{$split[1]} > log(1/$fc_significant)/log(2))) {
			next;
		}
		#Only keeps peaks where there is an effect, so when the delta between the tag counts of the two strains is smaller then the threshold set for significance (absolute value) the peak is ignored
		if($delta_tag == 1 && $effect == 1 && abs($fc{substr($split[0], 3)}{$split[1]}) < $delta_threshold) {
			next;
		}
		#Define the offset, because sequences are differently long but so the center of the sequence is different - the longest sequence is used so all sequences will be able to fit
		$offset = int(($longest_seq - length($seq->{$last_line . "_" . $strains[0] . "_1"}))/2);
		foreach my $motif (keys %{$block{$last_line}}) {
			foreach my $pos (keys %{$block{$last_line}{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					for(my $al = 1; $al <= $allele; $al++) {
						if(!exists $block{$last_line}{$motif}{$pos}{$strains[$i]}{$al}) {
							next;
						}
						#add 1 to counter whenever there is a motif at this position in a peak
						for(my $l = 0; $l < $block{$last_line}{$motif}{$pos}{'length'}; $l++) {
							$dist_plot{$strains[$i]}{($pos + $l + $offset)}++;
						#	$dist_plot{$motif}{$strains[$i]}{($pos + $l + $offset)}++;
						}
					}
				}
			}	
		}
	}
	return \%dist_plot;
}

#Method prepares output hash for the background distribution in the distance plots
sub background_dist_plot{
	$_ = () for my (%scan_candidates, %motifs, @fileHandlesBackground, @name, %background, %background_saved, %dist_plot_background, %missed_bg);
	$_ = "" for my ($motif_genomewide, $command, $tmp_motif, $tmp_motif2, $filename, $current_chr, $file);
	$_ = 0 for my ($i, $first, $count, $current_dist, $next_dist, $smallest_dist, $main_tf_start, $main_tf_end, $length);
	my $background_folder = $_[0];	
	my %index_motifs = %{$_[1]};
	my $delete = $_[2];
	my $motif_scan_score = $_[3];
	my $genome = $_[4];
	my $ab = $_[5];
	my $region = $_[6];
	my $longest_seq = $_[7];
	my $core = $_[8];
	#Once the motif was scanned in this genome it is saved - so check if the data is already there otherwise add it to the list of motifs to scan for
	print STDERR "Checking for genome wide scan of motifs\n";
	foreach my $motif (keys %index_motifs) {
		$motif_genomewide = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		if(!-e $motif_genomewide) {
			$scan_candidates{$motif} = 1;
		}
	}
	my $split_motifs = ceil((keys %scan_candidates)/$core);
	$count = 0;
	my @split_array;
	foreach my $key (sort {$a cmp $b} keys %scan_candidates) {
		$split_array[$count] = $key;
		$count++;
	}
	my @save_motifs;
	my $run = 0;
	my $k = 0;
	for(my $i = 0; $i < @split_array; $i = $i + $split_motifs) {
		$k = 0;
		for(my $j = $run * $split_motifs; $j < ($run + 1) * $split_motifs; $j++) {
			$save_motifs[$run][$k] = $split_array[$j];
			$k++;
		}
		$run++;
	}
	my $count_fork = 0;
	for(1 .. $core) {
		my $pid = fork;
		if(not $pid) {
			&process_background_scanning(\@{$save_motifs[$count_fork]}, $count_fork, $motif_scan_score, $genome, $background_folder);
			exit();
		}
		$count_fork++;
	}
	my $pid = fork;
	if(not $pid) {
		&monitor_scanning_process($count);
		exit();
	}
	for(1 .. ($core+1)) {
		wait();
	}
	#}
	#Read in all files with background distribution of motifs that are in the list - start with transcription factor that was use for the chip
	$motif_genomewide = $background_folder . "/" . $genome . "_" . $ab . ".txt";
	open(my $fh, "<", $motif_genomewide); 
	print STDERR "Saving chipped TF " . $motif_genomewide . "\n";
	#Save the positions of the TF used for the chip, which will be the reference for the distance of all the other TF
	while(my $line = <$fh>) {
		chomp $line;
		@split = split('\t', $line);
		if($first == 0) {
			$current_chr = $split[1];
			$first++;
		}
		if($split[1] ne $current_chr) {
			$count = 0;
		}
		$current_chr = $split[1];
		$background_saved{$current_chr}[$count] = $split[2];
		$background{$split[1]}{$split[2]} = $split[3];
		$count++;
	}
	close $fh;
	$count = 0;
	$current_chr = "";
	#After knowing where the main TF is, calculate the distance distribution of all the other TF
	foreach my $motif (keys %index_motifs) {
		print STDERR "Processing " . $motif . "\n";
		$file = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		open(my $fh, "<", $file);
		$count = 0;
		while(my $line = <$fh>) {
			chomp $line;
			#Start with first binding motif of the reference transcription factor and run through second TF till the distance of count is smaller than the distance of count + 1 (we start running away from the motif again - so this is the closest instance) 
			@split = split('\t', $line);
			if($split[1] ne $current_chr) {
				$count = 0;
			}
			if(!defined $background_saved{$split[1]}) {
				$missed_bg{$split[1]} = 1;
				next;
			}
			$current_chr = $split[1];
			$current_dist = $split[2] - $background_saved{$split[1]}[$count];
			if(defined $background_saved{$split[1]}[$count + 1]) {
			$next_dist = $split[2] - $background_saved{$split[1]}[$count + 1];
			} else {
				$next_dist = 100000000;
			}
			while(defined $background_saved{$split[1]}[$count + 1] && abs($split[2] - $background_saved{$split[1]}[$count]) > abs($split[2] - $background_saved{$split[1]}[$count + 1])) {
				$count++;
			}
			$current_dist = $split[2] - $background_saved{$split[1]}[$count];
			if(defined $background_saved{$split[1]}[$count + 1]) {
				$next_dist = $split[2] - $background_saved{$split[1]}[$count + 1];
			} else {
				$next_dist = 100000000;
			}
			#Determine which motif is closer to the main TF (the one to the left (current_dist) or to the right (next_dist)
			if(abs($current_dist) < abs($next_dist)) {
				$smallest_dist = $current_dist;
				$main_tf_start = $background_saved{$split[1]}[$count];
				$main_tf_end = $background{$split[1]}{$main_tf_start};
			} else {
				$smallest_dist = $next_dist;
				$main_tf_start = $background_saved{$split[1]}[$count + 1];
				$main_tf_end = $background{$split[1]}{$main_tf_start};
			}
			#Test if smallest distance is within the length of sequences that was defined by the longest sequence in the peak file
			if(abs($smallest_dist) < int($longest_seq/2) + $region) {
				$length = ($main_tf_end - $main_tf_start);
				#If this is the case add motif to background distribution plot
				for(my $i = int($length/2) * -1; $i < int($length/2); $i++) {
					$dist_plot_background{$motif}{$smallest_dist + $i}++;
				}
			}
		}
		close $fh;
		if(keys %missed_bg > 0) {
			print STDERR "Could not calculate a background for these chromosomes:\n";
			foreach my $k (keys %missed_bg) {
				print STDERR "\t" . $k . "\n";
			}
		}
		%missed_bg = ();
	}
	return ($delete, \%dist_plot_background);
}

#Method that analyzes the different bases within the motif to find which bases are mutated the most and have the most impact on the TF binding
sub analyze_motif_pos{
	$_ = 0 for my($current_shift, $prev_shift, $shift_vector, $chr_num, $current_pos, $max_motif, $num_of_muts, $max);
	$_ = () for my($tree_tmp, @fc_split, @char_one, @char_two, @header, %matrix_pos_muts, $last_line_ref, @last_line_split);
	my %block = %{$_[0]};
	my $fc = $_[1];
	my @strains = @{$_[2]};
	my $fc_significant = $_[3];
	my $seq = $_[4];
	my $tree = $_[5];
	my $last = $_[6];
	my $lookup = $_[7];
	my $motif_diff = $_[8];
	my $motif_diff_percentage = $_[9];
	my $delta_tag = $_[10];
	my $delta_threshold = $_[11];
	my $allele = $_[12];
	my $region = $_[13];
	my $chr_num_original;
	foreach my $last_line (keys %block) {
		chomp $last_line;
		@header = split("_", $last_line);
		$chr_num = substr($header[0], 3);
		$chr_num_original = $chr_num;
		@fc_split = split('\t', $fc->{$chr_num}->{$header[1]});
		$current_pos = $header[1] - ($region/2);
		for(my $i = 0; $i < @strains; $i++) {
			for(my $al = 1; $al <= $allele; $al++) {
				if($chr_num !~ /\d+/ && !exists $lookup->{$strains[$i]}->{$chr_num}) {
					print STDERR "Skip analysis of chromosome " . $chr_num . " allele " . $al . " in analyze_motifs_pos\n";
				}
				if(exists $lookup->{$strains[$i]}->{$chr_num}) {
					$chr_num = $lookup->{$strains[$i]}->{$chr_num};
				}
			}
		}
		foreach my $motif (keys %{$block{$last_line}}) {
			#Define indels, snps and multiple snps (N == not significant, S == significant) 
			if(!exists $matrix_pos_muts{'indel'}) {
				$matrix_pos_muts{'indel'}{'N'} = 0;
				$matrix_pos_muts{'indel'}{'S'} = 0;
			}
			if(!exists $matrix_pos_muts{'multi'}) {
				$matrix_pos_muts{'multi'}{'N'} = 0;
				$matrix_pos_muts{'multi'}{'S'} = 0;
			}

			#Run through all the motif occurrences in the peak file
			foreach my $motif_pos (keys %{$block{$last_line}{$motif}}) {
				#If not binary decisions but motif score difference define significane, the strain motif with the highest score needs to be found first
				if($motif_diff > 0 || $motif_diff_percentage > 0) {
					$max_motif = 0;
					for(my $i = 0; $i < @strains; $i++) {
						for(my $al = 1; $al <= $allele; $al++) {
							if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$al} > $max_motif) {
								$max_motif = $block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$al};
							}
						}
					}
				}
				#Run through all pairwise strain comparisons
				for(my $i = 0; $i < @strains - 1; $i++) {
					for(my $j = $i + 1; $j < @strains; $j++) {
						for(my $a1 = 1; $a1 <= $allele; $a1++) {
							for(my $a2 = 1; $a2 <= $allele; $a2++) {
								#Motif is different between the strains
								#motif differences and percentage of motif difference is not used to call different motifs - the score of the motif is 0 in one of the strains
								if(($motif_diff == 0 && $motif_diff_percentage == 0 && ($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} == 0 || $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2} == 0)) ||
								#Motif difference is used to call different motifs and the difference between the two motif scores is greater than what is defined as minmum different to be significant
								($motif_diff > 0 && abs($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} - $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2}) > $motif_diff) ||
								#Percentage of motif score difference is used to call different motifs and the percentual difference is bigger than the defined minimum to call it as significant
								($motif_diff_percentage > 0 && abs($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} - $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2}) > ($max_motif * $motif_diff_percentage))) {
									$shift_vector = 0;
									$current_shift = 0;
									$prev_shift = 0;
									#Check if there was an indel between the start of the sequence and the position of the motif in the strain because original sequences needs to be used in this method
									if(defined $tree->{$strains[$i]}->{$chr_num}->{$a1}) {
										if($current_pos + $motif_pos >= $last->{$strains[$i]}->{$chr_num_original}->{$a1}->{'pos'}) {
											$current_shift = $last->{$strains[$i]}->{$chr_num_original}->{$a1}->{'shift'};
										} else {
											#$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$a1}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos + 1);
											$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$a1}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos);
											$current_shift = $tree_tmp->[0]->{'shift'};
										}
										if($current_pos >= $last->{$strains[$i]}->{$chr_num_original}->{$a1}->{'pos'}) {
											$prev_shift = $last->{$strains[$i]}->{$chr_num_original}->{$a1}->{'shift'};
										} else {
											#$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$a1}->fetch($current_pos, $current_pos + 1);
											$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->{$a1}->fetch($current_pos, $current_pos);
											$prev_shift = $tree_tmp->[0]->{'shift'};
										}
										$shift_vector = $current_shift - $prev_shift;	
									}
									#Split the original strain motif sequence into an array, so a base by base comparison is possible
								#	print STDERR "Get char one\n";	
								#	print STDERR "last line: " . $last_line . "\tstrain: " . $strains[$i] . "\tallele: " . $a1. "\n";
								#	print STDERR "this seq: " . $seq->{$last_line . "_" . $strains[$i] . "_" . $a1} . "\n";
								#	print STDERR "motif pos: " . $motif_pos . "\tshift vector: " . $shift_vector . "\tlength: " . $block{$last_line}{$motif}{$motif_pos}{'length'} . "\n";
									if(length($seq->{$last_line . "_" . $strains[$i] . "_" . $a1}) < $motif_pos + $shift_vector +  $block{$last_line}{$motif}{$motif_pos}{'length'}){
										@char_one = ();
									} else {
										if($block{$last_line}{$motif}{$motif_pos}{'orientation'} eq "+") {
											@char_one = split("", substr($seq->{$last_line . "_" . $strains[$i] . "_" . $a1}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'}));
										} else {
											@char_one = split("", &rev_comp(substr($seq->{$last_line . "_" . $strains[$i] . "_" . $a1}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'})));
										}
									}
									$shift_vector = 0;
									$current_shift = 0;
									$prev_shift = 0;
									#Check if there was an indel between the start of the sequence and the position of the motif in the second strain
									if(defined $tree->{$strains[$j]}->{$chr_num}->{$a2}) {
										if($current_pos + $motif_pos >= $last->{$strains[$j]}->{$chr_num_original}->{$a2}->{'pos'}) {
											$current_shift = $last->{$strains[$j]}->{$chr_num_original}->{$a2}->{'shift'};
										} else {
											#$tree_tmp = $tree->{$strains[$j]}->{$chr_num}->{$a2}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos + 1);
											$tree_tmp = $tree->{$strains[$j]}->{$chr_num}->{$a2}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos);
											$current_shift = $tree_tmp->[0]->{'shift'};
										}
										if($current_pos >= $last->{$strains[$j]}->{$chr_num_original}->{$a2}->{'pos'}) {
											$prev_shift = $last->{$strains[$j]}->{$chr_num_original}->{$a2}->{'shift'};
										} else {
											#$tree_tmp = $tree->{$strains[$j]}->{$chr_num}->{$a2}->fetch($current_pos, $current_pos + 1);
											$tree_tmp = $tree->{$strains[$j]}->{$chr_num}->{$a2}->fetch($current_pos, $current_pos);
											$prev_shift = $tree_tmp->[0]->{'shift'};
										}
										$shift_vector = $current_shift - $prev_shift;	
									}
								#	print STDERR "Get char two\n";	
								#	print STDERR "last line: " . $last_line . "\tstrain: " . $strains[$j] . "\tallele: " . $a2. "\n";
								#	print STDERR "this seq: " . $seq->{$last_line . "_" . $strains[$j] . "_" . $a2} . "\n";
								#	print STDERR "motif pos: " . $motif_pos . "\tshift vector: " . $shift_vector . "\tlength: " . $block{$last_line}{$motif}{$motif_pos}{'length'} . "\n";
									#Split the original sequence of the second strain into an array for base by base comparison
									if(length($seq->{$last_line . "_" . $strains[$j] . "_" . $a2}) < $motif_pos + $shift_vector + $block{$last_line}{$motif}{$motif_pos}{'length'}) {
										@char_two = ();
									} else {
										if($block{$last_line}{$motif}{$motif_pos}{'orientation'} eq "+") {
											@char_two = split("", substr($seq->{$last_line . "_" . $strains[$j] . "_" . $a2}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'}));
										} else {
											@char_two = split("", &rev_comp(substr($seq->{$last_line . "_" . $strains[$j] . "_" . $a2}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'})));
										}
									}

									#Count number of differences - if more than 2 dismiss this comparison because motifs are not the same and count it as indel
									$num_of_muts = 0;
								#	my $t = Dumper @char_one;
								#	print STDERR $t . "\n";
								#	my $tt = Dumper @char_two;
								#	print STDERR "second:\n$tt\n";
									if(@char_one != @char_two) {
										$num_of_muts = (@char_two - @char_one);
									} else {
										for(my $c = 0; $c < @char_one; $c++) {
											if($char_one[$c] ne $char_two[$c]) {
												$num_of_muts++;
											}
										}
									}
								#	print "motif: " . $motif . "\n";
								#	print "strain1: " . join("", @char_one) . "\n";
								#	print "strain2: " . join("", @char_two) . "\n";
								#	print "num: " . $num_of_muts . "\n";
									if($num_of_muts > 2) {
										#Number of mutation is greater than 2 -> mutation is counted as indel - check FC in order to define if mutation changed the binding of TF
										if($delta_tag == 0) {
											if($fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1] > log($fc_significant)/log(2) || $fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1] < log(1/$fc_significant)/log(2)) {
												$matrix_pos_muts{'indel'}{'S'}++;
											} else {
												$matrix_pos_muts{'indel'}{'N'}++;
											}
										} else {
											if(abs($fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1]) > $delta_threshold) {
												$matrix_pos_muts{'indel'}{'S'}++;
											} else {
												$matrix_pos_muts{'indel'}{'N'}++;
											}
										}

									} elsif($num_of_muts > 1) {
									#Mnumber of mutation is greather than 1 (and smaller than 2) -> count it as multiple mutations within the same motif - check FC in order to define if mutation changed the binding of TF
										if($delta_tag == 0) {
											if($fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1] > log($fc_significant)/log(2) || $fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1] < log(1/$fc_significant)/log(2)) {
												$matrix_pos_muts{'multi'}{'S'}++;
											} else {
												$matrix_pos_muts{'multi'}{'N'}++;
											}
										} else {
											if(abs($fc_split[$i] + ($a1 - 1) + $j + ($a2 - 1) - 1) > $delta_threshold) {
												$matrix_pos_muts{'multi'}{'S'}++;
											} else {
												$matrix_pos_muts{'multi'}{'N'}++;
											}
										}


									} elsif($num_of_muts == 1) {
										#Number of mutation is exactly 1 - check the kind of subsitution and determine if it was significant or not
										for(my $c = 0; $c < @char_one; $c++) {
											if($char_one[$c] ne $char_two[$c]) {
												if($delta_tag == 0) {
													if($fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1] > log($fc_significant)/log(2) || $fc_split[$i + ($a1 - 1) + $j + ($a2 - 1) - 1] < log(1/$fc_significant)/log(2)) {
														#S meaning significant
														if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2}) {
															$matrix_pos_muts{$char_two[$c]}{$c}{'S'}++;
														} else {
															$matrix_pos_muts{$char_one[$c]}{$c}{'S'}++;
														}
													} else {
														if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2}) {
															$matrix_pos_muts{$char_two[$c]}{$c}{'N'}++;
														} else {
															$matrix_pos_muts{$char_one[$c]}{$c}{'N'}++;
														}
													}
												} else {
													if(abs($fc_split[$i + $j - 1]) > $delta_threshold) {
														#S meaning significant
														if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2}) {
															$matrix_pos_muts{$char_two[$c]}{$c}{'S'}++;
														} else {
															$matrix_pos_muts{$char_one[$c]}{$c}{'S'}++;
														}
													} else {
														if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]}{$a1} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}{$a2}) {
															$matrix_pos_muts{$char_two[$c]}{$c}{'N'}++;
														} else {
															$matrix_pos_muts{$char_one[$c]}{$c}{'N'}++;
														}
													}

												}

											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return \%matrix_pos_muts;
}

#Calculate the score of the motif
sub calculate_motif_score{
	my $score = 0;
	my $motif = $_[0];
	my $local_seq = $_[1];
	my $orientation = $_[2];
	if($local_seq =~ m/N/) { return 0; }
	#Sometimes motif is missing because of indel
	if(!defined $local_seq || length($local_seq) == 0) {
		$score = 0;
		print STDERR "Return\n";
		return sprintf("%.6f", $score);	
	}
	my @base = split('', $local_seq);
#	print STDERR "we keep going\n";
	if($orientation eq "+") {
		for(my $b = 0; $b < @base; $b++) {
			$score += log($PWM{$motif}{$b}{$base[$b]}/0.25); 
		}
	} else {
		$local_seq = &rev_comp($local_seq);	
		@base = split('', $local_seq);
		for(my $b = 0; $b < @base; $b++) {
			$score += log($PWM{$motif}{$b}{$base[$b]}/0.25); 
		}
	}
	if($score < 0) { $score = 0;}
	return sprintf("%.6f", $score);	
}

#Generate reverse complementary sequence
sub rev_comp{
	my $local_seq = $_[0];
	my $rev = reverse($local_seq);
	my @split = split('', $rev);
	my $comp = "";
	foreach my $c (@split) {
		$comp .= $comp{$c};
	}	
	return $comp;
}

#Print motif name from HOMER motif file in are more human readable format - maybe skip for other motif files
sub get_motif_name{
        my $name = $_[0];
	$_ = () for my (@clean1, @clean2, @underscore, @paren);
	#First check if it comes from homer motif analysis or database
	if(index($name, "BestGuess") != -1) {
		#Get rid of consensus sequence at the beginning
		@clean1 = split("BestGuess:", $name);
		#Get rid of best guess
		#Get rif of everything in paraenthesis 
		@clean2 = split/[\(,\/,\:]+/, $clean1[-1];
		$name = $clean2[0];
		@underscore = split('_', $clean2[0]);
		for(my $i = @underscore - 1; $i > 0; $i--) {
			if(length($underscore[$i]) > 2 && $underscore[$i] =~ /[a-zA-Z]+/) {
				$name = $underscore[$i];
				return $name;
				last;
			}
		}
		@paren = split('\(', $name);
		$name = $paren[0];
	} else {
		@clean1 = split(':', $name);
		if(@clean1 > 1) {
			$name = $clean1[0];
		} else {
			@clean1 = split/[\,,\:]+/, $name;
			if(@clean1 > 1) {
				$name = $clean1[1];
			} else {
				$name = $clean1[0];
			}
			@clean2 = split/[\(,\/,\:]+/, $name;
			$name = $clean2[0];
		}
	}
	$name =~ s/\-/_/g;
	$name =~ s/ //g;
	$name =~ s/://g;
	$name =~ s/\)//g;
	$name =~ s/\(//g;
	$name =~ s/\|/_/g;
	$name =~ s/\?//g;
	$name =~ s/\+//g;
	return $name;
}

#Open filehandles for writing
sub open_filehandles{
	my %motifs = %{$_[0]};
	my %del = %{$_[1]};
	my $output_mut = $_[2];
	my @fileHandlesMotif;
	my $index = 0;
	foreach my $motif (keys %motifs) {
		my $filename = $output_mut . "_" . $motif . ".txt";
		$del{$filename} = 1;
	#	my $fh = cacheout ">", $filename or die "Can'topen $filename: $!\n"; 
		open my $fh, ">", $filename or die "Can't open $filename: $!\n";
		$fileHandlesMotif[$motifs{$motif}] = $fh;
	}
	return (\@fileHandlesMotif, \%del);
}

#Read motifs from motif file and save PWM in a global variable that every method in this package can access - should be changed at some point
sub read_motifs{
	open FH, "<$_[0]" or die "Could not open motif file $_[0]: $!\n";
	my $pos = 0;
	my $motif = "";
	my $index = 0;
	my $name;
	my %index_motifs;
	my %score;
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			$motif = (split('\t', $line))[1];
			$pos = 0;
			#Open file per motif
			$name = &get_motif_name($motif);
			#Read in when not already in index_motifs (different motifs with same name are overwritten - maybe change at some point)
			if(!exists $index_motifs{$name}) {
				$index_motifs{$name} = $index;
				$index++;
				@split = split('\t', $line);
				$score{$name} = $split[2];
			}
		} else {
			#Change positions with probabailty 0 to 0.001 (easier when calculating motif score)
			@split = split('\t', $line);
			if($split[0] == 0) { $split[0] += 0.001; }
			if($split[1] == 0) { $split[1] += 0.001; }
			if($split[2] == 0) { $split[2] += 0.001; }
			if($split[3] == 0) { $split[3] += 0.001; }
			$PWM{$name}{$pos}{'A'} = $split[0];
			$PWM{$name}{$pos}{'C'} = $split[1];
			$PWM{$name}{$pos}{'G'} = $split[2];
			$PWM{$name}{$pos}{'T'} = $split[3];
			$pos++;
		}
	}
	return (\%index_motifs, \%PWM, \%score);
}

#Scan the sequecnes with HOMER to find motifs
sub scan_motif_with_homer{
	my $seq_file = $_[0];
	my $out_file = $_[1];
	my $tf = $_[2];
	my $rfolder = $_[3];
	my $wd = $_[4];
#	my $command = "homer2 find -i " . $seq_file . " -m " . $tf . " -offset 0 > " . $out_file . " 2> /dev/null";
	my $command = "cd " . $rfolder . " && homer2 find -i " . $wd . "/" . $seq_file . " -m " . $wd . "/" . $tf . " -offset 0 -o " . $wd . "/" . $out_file . " 2> /dev/null; cd -";
	`$command`;	
}

#Get strain sequences for every peak
sub get_seq_for_peaks {
	$_ = 0 for my ($general_offset, $seq_number, $newlines, $longest_seq, $no_shift, $filter_no_mut, $length, $ref_start, $working_start, $working_stop, $byte_offset, $chr_num, $shift_vector, $mut, $allele_num); 
	$_ = "" for my($seq, $h, $header, $seq_first, $filename);
	$_ = () for my(@fileHandles, %current_pos, @tmp_split, %seen_part, @split_part, %seq); 
	my $current_tree = Set::IntervalTree->new;
	open OUT, ">$_[0]";
	my %peaks = %{$_[1]};
	my @strains = @{$_[2]};
	my $data = $_[3];
	my $allele = $_[4];
	my $line_number = $_[5];
	my $mut_only = $_[6];
	my $region = $_[7];
	my $tree = $_[8];
	my $lookup = $_[9];
	my $last = $_[10];
	my %strand;
	if(@_ > 11) {
		%strand = %{$_[11]};
	}
	my $first_line;
	$line_number = $line_number * @strains;
	for(my $i = 0; $i < @strains; $i++) {
		foreach my $chr (keys %peaks) {
			for(my $al = 1; $al <= $allele; $al++) {
				$chr_num = $chr;
				$allele_num = $al;
				$no_shift = 0;
				#Start with offset for first line
			#	$filename = $data . "/" . uc($strains[$i]) . "/chr" . $chr . "_allele_1.fa";
				$filename = $data . "/" . uc($strains[$i]) . "/chr" . $chr . "_allele_" . $allele_num . ".fa";
				if(!-e $filename) {
					print STDERR "Can't open $filename\n";
					print STDERR "Exclude chromosome $chr allele $allele_num from analysis\n";
					delete $peaks{$chr};
					#Also remove it from sequence 
					foreach my $s (%seq) {
						if($s eq "") { next; }
						@tmp_split = split("_", $s);
						$h = $tmp_split[0];
						for(my $k = 1; $k < @tmp_split - 3; $k++) {
							$h .= "_" . $tmp_split[$k];
						}
						if($h eq "chr" . $chr) {
							delete $seq{$s};
						}
					}
					next;	
				}	
				$first_line = `head -n1 $filename`;
			#	$general_offset = 5 + length($chr);
				$general_offset = length($first_line);
				if(exists $lookup->{$strains[$i]}->{$chr_num}) {
					$chr_num = $lookup->{$strains[0]}->{$chr};
				}
				if($chr !~ /\d+/ && !exists $lookup->{$strains[$i]}->{$chr}) {
				#	print STDERR "Skip peaks for chromosome " . $chr . " in get_seq_for_peaks\n";
					print STDERR "Skip peaks for chromosome " . $chr . " allele " . $allele_num . " in get_seq_for_peaks\n";
					next;
				}
			#	$current_tree = $tree->{$strains[$i]}->{$chr_num};
				$current_tree = $tree->{$strains[$i]}->{$chr_num}->{$allele_num};
				if(!defined $current_tree) {
					$no_shift = 1;
				}
				open my $fh, "<", $filename or die "Can't open $filename: $!\n";
				#Save file in filehandle to access it via byteOffsets
				$fileHandles[0] = $fh;
				open FH, "<$filename";
				foreach my $start_pos (sort {$a <=> $b} keys %{$peaks{$chr}}) {
					#Define start and stop pos from sequence with region
					$working_start = $start_pos - int($region/2);
					$working_stop = $peaks{$chr}{$start_pos} + int($region/2);
					#Calculate strain specific offset from shifting vector
					if($no_shift == 0) {
						#$ref_start = $current_tree->fetch($working_start, $working_start + 1);
						$ref_start = $current_tree->fetch($working_start, $working_start);
						if(!exists $ref_start->[0]->{'shift'}) {
							$shift_vector = $last->{$strains[$i]}->{$chr}->{$allele_num}->{'shift'};
						} else {
							$shift_vector = $ref_start->[0]->{'shift'};
						}
						$working_start = $working_start + $shift_vector;
						#$ref_start = $current_tree->fetch($working_stop, $working_stop + 1);
						$ref_start = $current_tree->fetch($working_stop, $working_stop);
						if(!exists $ref_start->[0]->{'shift'}) {
							$shift_vector = $last->{$strains[$i]}->{$chr}->{$allele_num}->{'shift'};
						} else {
							$shift_vector = $ref_start->[0]->{'shift'};
						}
						$working_stop = $working_stop + $shift_vector;
					}
				#	$header = "chr" . $chr . "_" . $start_pos . "_" . $peaks{$chr}{$start_pos} . "_" . $strains[$i];
					$header = "chr" . $chr . "_" . $start_pos . "_" . $peaks{$chr}{$start_pos} . "_" . $strains[$i] . "_" . $allele_num;
					#Calculate length and newlines (fastq file saves seq in 50bp lines) 
					$length = $working_stop - $working_start;
					$newlines = int($working_stop/50) - int($working_start/50);
					#Add number of newlines to length (because working on bytes)
					$length = $length + $newlines;
					#Calculate the number of newlines before before working start in sequence
					$newlines = int(($working_start)/50);
					#Calculate offset - general offset is first line
					$byte_offset = $working_start + $general_offset + $newlines;
					#When we are looking at a sequences that is totally gone in the current strain then length can be negative
					if($length < 0) { $length = 1;} 
					seek($fileHandles[0], $byte_offset, 0);
					read $fileHandles[0], $seq, $length;
					$seq =~ s/\n//g;
					if(length($seq) > $longest_seq) {
						$longest_seq = length($seq);
					}
					$seq = uc($seq);
					if(exists $strand{$chr} && ($strand{$chr}{$start_pos} && $strand{$chr}{$start_pos} =~ /^[0-9,.E]+$/ && $strand{$chr}{$start_pos} == 1) || ($strand{$chr}{$start_pos} && $strand{$chr}{$start_pos} eq "-")) {
						$seq = &rev_comp($seq);
					}
					$current_pos{$strains[$i]} = $seq;
					$seq{$header} = uc($seq);
					$seq = "";
					$seq_number++;
				}
			}
		}
	}
	#Only peaks with mutations are considered in this analysis
	if($mut_only == 1) {
		print STDERR "Filter out peaks without mutations\n";
		foreach my $key (keys %seq) {
			@split_part = split("_", $key);
			$h = $split_part[0] . "_" . $split_part[1] . "_" . $split_part[2];
			if(exists $seen_part{$h}) { next; }
			for(my $al = 1; $al <= $allele; $al++) {
				$seq_first = $seq{$h . "_" . $strains[0] . "_" . $al};
				#Compare all strains to first strain - if one sequence is different keep sequences
				for(my $i = 1; $i < @strains; $i++) {
					if($seq_first ne $seq{$h . "_" . $strains[$i] . "_" . $al}) {
						$mut = 1;
					} 
				}
			}
			if($mut == 0) {
				$filter_no_mut++;
				#No mutation was found - delete this peak
				for(my $i = 0; $i < @strains; $i++) {
					for(my $al = 1; $al <= $allele; $al++) {
						delete $seq{$h . "_" . $strains[$i] . "_" . $al};
					}
				}
			}
			$mut = 0;
			$seen_part{$h} = 1;
		}
	}
	if($filter_no_mut > 0) {
		print STDERR "" . $filter_no_mut . " sequences were filtered out because of no variance between the strains\n";
	}
	#Save sequences in a file so HOMER can scan sequences for motifs
	foreach my $key (sort {$a cmp $b} keys %seq) {
		print OUT ">" . $key . "\n";
		print OUT $seq{$key}  . "\n";
	}
	close OUT;
	return(\%seq, $longest_seq, $filter_no_mut);
}


sub all_vs_all_comparison {
	my $block = $_[0];
	my $recenter_conversion = $_[1];
	my $tag_counts = $_[2];
	my @strains = @{$_[3]};
	my $allele = $_[4];
	my $score = $_[5];
	my $output = $_[6];
	my $same = 0;
	my $value = 0;
	my $all = keys %{$block};
	my $count = 0;
	my $filename;
	$_ = () for my(@a, %motif_matrix, %motif_matrix_number, %tag_matrix, @tag_sum);
	foreach my $pos (sort {$a cmp $b} keys %{$block}) {
		@a = split("_", $pos);
		@tag_sum = split('\s+', $tag_counts->{substr($a[0], 3)}->{$a[1]});
		for(my $i = 0; $i < @strains; $i++) {
			for(my $al = 1; $al <= $allele; $al++) {
				$tag_matrix{$pos}{$strains[$i]}{$al} = $tag_sum[$i];
			}
		}
		foreach my $motif (sort {$a cmp $b} keys %{$block->{$pos}}) {
			foreach my $motif_pos (sort {$a <=> $b} keys %{$block->{$pos}->{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					for(my $al = 1; $al <= $allele; $al++) {
						if(!exists $motif_matrix{$motif}{$pos}{$strains[$i]}{$al}) {
							$motif_matrix{$motif}{$pos}{$strains[$i]}{$al} = 0;
							$motif_matrix_number{$motif}{$pos}{$strains[$i]}{$al} = 0;
						}
						if($block->{$pos}->{$motif}->{$motif_pos}->{$strains[$i]}->{$al} > $score->{$motif}) {
					#	if($block->{$pos}->{$motif}->{$motif_pos}->{$strains[$i]}->{$al} > 5) {
							$motif_matrix{$motif}{$pos}{$strains[$i]}{$al}++;
						}
						$motif_matrix_number{$motif}{$pos}{$strains[$i]}{$al} += $block->{$pos}->{$motif}->{$motif_pos}->{$strains[$i]}->{$al};
					}
				}
			}
		}
	}
	foreach my $motif (keys %motif_matrix) {
		$filename = $output . "_" . $motif . ".txt";
		open OUT, ">$filename";
		print OUT "Strain\tLocus\tMotif\tMotif_score\tbinding\n";
		foreach my $pos (sort {$a cmp $b} keys %tag_matrix) {
			for(my $i = 0; $i < @strains; $i++) {
				for(my $al = 1; $al <= $allele; $al++) {
					print OUT $strains[$i] . "_" . $al . "\t" . $pos . "\t";
					if(!exists $motif_matrix{$motif}{$pos}{$strains[$i]}{$al}) {
						print OUT "0\t0";
					} else {
						print OUT $motif_matrix{$motif}{$pos}{$strains[$i]}{$al} . "\t" . $motif_matrix_number{$motif}{$pos}{$strains[$i]}{$al};
					}	
					print OUT "\t" . $tag_matrix{$pos}{$strains[$i]}{$al} . "\n";
				}
			}
		}
		close OUT;
	}
}

sub split_into_single_files{
	my $file = $_[0];
	open (my $fh, "<$file");
	my $current_motif = "";
	my $short_name;
	my $open = 0;
	while (my $line = <$fh>) {
		chomp $line;
		@split = split('\t', $line);
		$short_name = &get_motif_name($split[3]);		
		if($short_name ne $current_motif) {
			if($open == 1) {
				close OUT;
			}
			open OUT, ">" . $file . "_" . $short_name;
			$current_motif = $short_name;
			$open = 1;
		}
		print OUT $line . "\n";
	}
	close OUT;
	close($fh);
}

sub process_background_scanning {
	my @fork_motifs = @{$_[0]};
	my $c = $_[1];
	my $motif_scan_score = $_[2];
	my $genome = $_[3];
	my $background_folder = $_[4];
	$_ = () for my (%delete, $motif, $tmp_motif, $tmp_motif2, $command, $filename, $delete, @name);
	for(my $i = 0; $i < @fork_motifs; $i++) {
		if(!defined $fork_motifs[$i]) { next; }
		$motif = $fork_motifs[$i];
		$tmp_motif = "tmp" . rand() . ".txt";
		$delete->{$tmp_motif} = 1;
		open TMP, ">$tmp_motif";
	#	foreach my $motif (keys %scan_candidates) {
			print TMP ">consensus_" . $motif . "\t$motif\t" . $motif_scan_score->{$motif} . "\n";
			foreach my $pos (sort {$a <=> $b} keys %{$PWM{$motif}}) {
				print TMP $PWM{$motif}{$pos}{'A'} . "\t" . $PWM{$motif}{$pos}{'C'} . "\t" . $PWM{$motif}{$pos}{'G'} . "\t" . $PWM{$motif}{$pos}{'T'} . "\n";
			}
	#	}
		close TMP;

		print STDERR "Scan motifs genome wide\n";
		$tmp_motif2 = "tmp" . rand();
		$delete->{$tmp_motif2} = 1;
		$command = "scanMotifGenomeWide.pl " . $tmp_motif . " " . $genome . " > " . $tmp_motif2 . " 2> error_" . $tmp_motif2;
		$delete->{"error_" . $tmp_motif2} = 1;
		print STDERR $command . "\n";
		`$command`;
		open(my $fh, "<", $tmp_motif2);
		#Save motifs in file
		#Open filehandle for the background file for every motif
#		foreach my $motif (keys %scan_candidates) {
		$filename = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		open(my $fh_out, ">", $filename) or die "Can't open $filename: $!\n";
#			$motifs{$motif} = $i;
#			$fileHandlesBackground[$motifs{$motif}] = $fh;
#			$i++;
#		}
		#Write output of motif scan in all the different files
		while(my $line = <$fh>) {
			chomp $line;
			@split = split('\t', $line);
			@name = split('-', $split[0]);
		#	$fileHandlesBackground[$motifs{$name[0]}]->print($split[1] . "\t" . $split[2] . "\t" . $split[3] . "\n");
			print $fh_out $line . "\n";
		#	$fileHandlesBackground[$motifs{$name[0]}]->print($line . "\n");
		}
		close $fh_out;
#		foreach my $motif (keys %scan_candidates) {
#			close $fileHandlesBackground[$motifs{$motif}];
#		}
		close $fh;

	}
}

sub monitor_scanning_process{
	my $final_number = 0;
	my $all = $_[0];
	my $lines;
	my @n;
	print STDERR "\n\n";
	while($final_number < $all) {
		$lines = `ls consensus_* 2> error_tmp | wc -l`;
		chomp $lines;
		@n = split('\t', $lines);
		$final_number = $n[0];
		print STDERR "\t\t" . $final_number . " of " . $all . " backgrounds generated so far - " . sprintf("%.2f", (($final_number/$all)*100)) . "\r";
	}
}
