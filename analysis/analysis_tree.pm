#!/usr/bin/perl
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
package analysis;
use strict;
use config;
use Statistics::Basic qw(:all);
use Statistics::Distributions;
use Statistics::RankCorrelation;;
#use Statistics::TTest;
use Data::Dumper;
use Storable;
use Statistics::R;

$_ = () for my(@split, %PWM, @tmp_split, %comp);
$comp{'A'} = 'T';
$comp{'T'} = 'A';
$comp{'C'} = 'G';
$comp{'G'} = 'C';
$comp{'-'} = '-';
$comp{'N'} = 'N';

#Write header for mutation distribution plots
sub write_header{
	my @fileHandlesMotif = @{$_[0]};
	my @strains = @{$_[1]};
	my $fc = $_[2];
	my $delta = $_[3];
	for(my $files = 0; $files < @fileHandlesMotif; $files++) {
		$fileHandlesMotif[$files]->print("motif\tpos\t");
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$files]->print("tg " . $strains[$i] . "\t");
		}
		#Check if foldchange is set 
		if($fc == 0) {
			for(my $i = 0; $i < @strains - 1; $i++) {
				for(my $j = $i + 1; $j < @strains; $j++) {
					if($delta == 1) {
						$fileHandlesMotif[$files]->print("delta of tg" . $strains[$i] . " vs " . $strains[$j] . "\t");
					} else {
						$fileHandlesMotif[$files]->print("log2 fc ". $strains[$i] . " vs " . $strains[$j] . "\t");
					}
				}
			}
		}
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$files]->print("motif_score " . $strains[$i] . "\t");
		}
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$files]->print("exist " . $strains[$i] . "\t");
		}
		$fileHandlesMotif[$files]->print("#motifs\n");
	}
}

#This method summarizes the output of HOMER's motif scan for every peak. It generates a hash that contains all motifs for all strains per peak
#The method reports the motif positions from the beginning of the sequence and takes care of shifting through indels within this sequence
#All motif positions are realtive to the reference motif position
sub analyze_motifs{
	$_ = () for my(@split, @header, %block, $tree_tmp, $last_line);
	$_ = 0 for my($pos_beginning, $shift_beginning, $pos_current, $pos_end, $shift_current, $shift_diff, $shift_end, $motif_pos, $chr_num, $motif_start);
	my $motif_file = $_[0];
	my @strains = @{$_[1]};
	my $tree = $_[2];
	my $lookup = $_[3];
	my $last = $_[4];
	my $allele = $_[5];
	my $region = $_[6];
        open FH, "<$motif_file";
        foreach my $line (<FH>) {
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
		if($chr_num !~ /\d+/ && !exists $lookup->{$header[3]}->{$chr_num}) {
			print STDERR "Skip analysis of chromosome " . $chr_num . " in analyze_motifs\n";
			print STDERR $chr_num . "\t" . $header[3] . "\n";
		}
		if(exists $lookup->{$header[3]}->{$chr_num}) {
			$chr_num = $lookup->{$header[3]}->{$chr_num};
		}
		#Check if there is a interval tree for this strain and chromosome
		if(defined $tree->{$header[3]}->{$chr_num}) {
			#Get shifting vector for the beginning of the sequence that was used to scan for the motif
			$tree_tmp = $tree->{$header[3]}->{$chr_num}->fetch($pos_beginning, $pos_beginning + 1);
			if(!exists $tree_tmp->[0]->{'shift'}) {
				$shift_beginning = $last->{$header[3]}->{$header[0]}->{$allele}->{'shift'};
			} else {
				$shift_beginning = $tree_tmp->[0]->{'shift'};
			}
			#Get shifting vector for the position of the motif
			$tree_tmp = $tree->{$header[3]}->{$chr_num}->fetch($motif_start, $motif_start + 1);
			if(!exists $tree_tmp->[0]->{'shift'}) {
				$shift_current = $last->{$header[3]}->{$header[0]}->{$allele}->{'shift'};
			} else {
				$shift_current = $tree_tmp->[0]->{'shift'};
			}
			#Calculate the difference between the shfiting vector at the beginning of the sequence and the shifting vector at the beginning of the motif - if different then there is a InDel in this sequence and shifting is necessary
			$shift_diff = $shift_beginning - $shift_current;
			$motif_pos = $motif_pos + $shift_diff;
		}
		#Save motif, length and orientation for this strain at the position in a hash
		$block{$last_line}{$split[3]}{$motif_pos}{$header[-1]} = $split[-1];
		$block{$last_line}{$split[3]}{$motif_pos}{'length'} = length($split[2]);
		$block{$last_line}{$split[3]}{$motif_pos}{'orientation'} = $split[4];
        }
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
					$shift_vector = 0;
					$chr_num = substr((split('_', $chr_pos))[0], 3);
					if($chr_num !~ /\d+/ && !exists $lookup->{$strains[$i]}->{$chr_num}) {
						print STDERR "Skip analysis of chromosome " . $chr_num . " in merge_block\n";
						next;
					}
					if(exists $lookup->{$strains[$i]}->{$chr_num}) {
						$chr_num = $lookup->{$strains[$i]}->{$chr_num};
					}
					if(!exists $diff{$strains[$i]}) {
						$diff{$strains[$i]} = 0;
					}
					#Motif does not exist at this position in this strain
					if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}) {
						#Check in vicinity (+/- motif_lenght-1) - if there is a motif within motif/2 -> save the current motif under this position, to not count them as separate motifs
						$existance{$strains[$i]} = 1;
						#Check if there is an overlap of a motif nearby in this strain (e.g. in this strain there was a indel that shifted the motif by some basepairs)
						#If overlap is not set, $number_of_bp is 0
						for(my $run_motif = $motif_pos - $num_of_bp; $run_motif < $motif_pos + $num_of_bp; $run_motif++) {
							#There is a motif nearby and the current position we are looking it is not the motif position 
							if(exists $block{$chr_pos}{$motif}{$run_motif} && $run_motif != $motif_pos) {
								#Run through all strains now and shift this motif to the current motif position and then delete the other motif position
								for(my $run_strains = 0; $run_strains < @strains; $run_strains++) {
									if(!exists $block{$chr_pos}{$motif}{$run_motif}{$strains[$run_strains]}) { next; }
									$block{$chr_pos}{$motif}{$motif_pos}{$strains[$run_strains]} += $block{$chr_pos}{$motif}{$run_motif}{$strains[$run_strains]};
									$diff{$strains[$run_strains]} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$run_strains]};
									delete $existance{$strains[$run_strains]};
									delete $block{$chr_pos}{$motif}{$run_motif}{$strains[$run_strains]};
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
					if($block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} eq "" || !exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}) {
						#Define the position of this motif in the strain that is currently investigated
						#Then pull this sequence from the sequence hash and calculate the motif score
						if(defined $tree->{$strains[$i]}->{$chr_num}) {
							$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos + 1);
							$current_shift = $tree_tmp->[$i]->{'shift'};
							$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->fetch($current_pos, $current_pos + 1);
							$prev_shift = $tree_tmp->[$i]->{'shift'};
							$shift_vector = $current_shift - $prev_shift;	
						}
						$block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} = &calculate_motif_score($motif, substr($seq->{$chr_pos . "_" . $strains[$i]}, $motif_pos + $shift_vector, $block{$chr_pos}{$motif}{$motif_pos}{'length'}), $block{$chr_pos}{$motif}{$motif_pos}{'orientation'});
						$diff{$strains[$i]} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]};
					} else {
						$diff{$strains[$i]} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]};
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
	my @fileHandlesMotif = @{$_[1]};
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
        #Now run through all motifs and see if they are in all the other strains
	foreach my $chr_pos (keys %block) {
		@split = split("_", $chr_pos);
		$curr_chr = substr($split[0], 3);
		$curr_pos = $split[1];
		if(exists $recenter_conversion{$curr_chr}{$curr_pos}{'start'}) {
			$curr_pos = $recenter_conversion{$curr_chr}{$curr_pos}{'start'};
		}
		foreach my $motif (keys %{$block{$chr_pos}}) {
			%diff = 0;
			%existance = ();
			$fileHandlesMotif[$index_motifs{$motif}]->print($motif . "\t" . $chr_pos . "\t");
			$fileHandlesMotif[$index_motifs{$motif}]->print($tag_counts{$curr_chr}{$curr_pos});
			foreach my $motif_pos (keys %{$block{$chr_pos}{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					#Filter by motif existance - either motif is not there at all (default) or diff_motif is set
					if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} || ($motif_diff == 0 && $motif_diff_percentage == 0 && $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} == 0)) {
						$existance{$strains[$i]} = 1;
					}
					if($block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} > $max_motif) {
						$max_motif = $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]};
					}
					$diff{$strains[$i]} = $diff{$strains[$i]} + $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]};
				}
				#not binary motif existance, but difference in motif score is used to determine if there is a mutation within the peak
				if($motif_diff != 0) {
					for(my $i = 0; $i < @strains; $i++) {
						if($max_motif - $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} > $motif_diff) {
							$existance{$strains[$i]} = 1;
						}
					}
				}
				#not binary motif existance or absolute difference in motif score is used, but difference percentage between the two motif scores
				if($motif_diff_percentage != 0) {
					for(my $i = 0; $i < @strains; $i++) {
						if($max_motif - $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} > ($max_motif * $motif_diff_percentage)) {
							$existance{$strains[$i]} = 1;
						}
					}
				}
				$max_motif = 0;
			}
			if($fc_exists == 1) {
				$fileHandlesMotif[$index_motifs{$motif}]->print($fc{$curr_chr}{$curr_pos});
			}
			for(my $i = 0; $i < @strains; $i++) {
				$fileHandlesMotif[$index_motifs{$motif}]->print($diff{$strains[$i]} . "\t");
			}
			#output number of mutations in peak for this motif
			for(my $i = 0; $i < @strains; $i++) {
				if(exists $existance{$strains[$i]}) {
					$fileHandlesMotif[$index_motifs{$motif}]->print("1\t");
				} else {
					$fileHandlesMotif[$index_motifs{$motif}]->print("0\t");
				}
			}
			$fileHandlesMotif[$index_motifs{$motif}]->print("" . (keys %{$block{$chr_pos}{$motif}}) . "\n");
		}
	}
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
		$offset = int(($longest_seq - length($seq->{$last_line . "_" . $strains[0]}))/2);
		foreach my $motif (keys %{$block{$last_line}}) {
			foreach my $pos (keys %{$block{$last_line}{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					if(!exists $block{$last_line}{$motif}{$pos}{$strains[$i]}) {
						next;
					}
					#add 1 to counter whenever there is a motif at this position in a peak
					for(my $l = 0; $l < $block{$last_line}{$motif}{$pos}{'length'}; $l++) {
						$dist_plot{$motif}{$strains[$i]}{($pos + $l + $offset)}++;
					}
				}
			}	
		}
	}
	return \%dist_plot;
}

#Method prepares output hash for the background distribution in the distance plots
sub background_dist_plot{
	$_ = () for my (%scan_candidates, %motifs, @fileHandlesBackground, @name, %background, %background_saved, %dist_plot_background);
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
	#Once the motif was scanned in this genome it is saved - so check if the data is already there otherwise add it to the list of motifs to scan for
	print STDERR "Checking for genome wide scan of motifs\n";
	foreach my $motif (keys %index_motifs) {
		$motif_genomewide = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		if(!-e $motif_genomewide) {
			$scan_candidates{$motif} = 1;
		}
	}
	#Scan for motifs that were not already proprocessed
	if(keys %scan_candidates > 0) {
		$tmp_motif = "tmp" . rand(15) . ".txt";
		$delete->{$tmp_motif} = 1;
		open TMP, ">$tmp_motif";
		foreach my $motif (keys %scan_candidates) {
			print TMP ">consensus_" . $motif . "\t$motif\t" . $motif_scan_score->{$motif} . "\n";
			foreach my $pos (sort {$a <=> $b} keys %{$PWM{$motif}}) {
				print TMP $PWM{$motif}{$pos}{'A'} . "\t" . $PWM{$motif}{$pos}{'C'} . "\t" . $PWM{$motif}{$pos}{'G'} . "\t" . $PWM{$motif}{$pos}{'T'} . "\n";
			}
		}
		close TMP;

		print STDERR "Scan motifs genome wide\n";
		$tmp_motif2 = "tmp" . rand(15);
		$delete->{$tmp_motif2} = 1;
		$command = "scanMotifGenomeWide.pl " . $tmp_motif . " " . $genome . " > " . $tmp_motif2;
		`$command`;
		open FH, "<$tmp_motif2";
		#Save motifs in file
		#Open filehandle for the background file for every motif
		foreach my $motif (keys %scan_candidates) {
			$filename = $background_folder . "/" . $genome . "_" . $motif . ".txt";
			open my $fh, ">", $filename or die "Can't open $filename: $!\n";
			$motifs{$motif} = $i;
			$fileHandlesBackground[$motifs{$motif}] = $fh;
			$i++;
		}
		#Write output of motif scan in all the different files
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			@name = split('-', $split[0]);
			$fileHandlesBackground[$motifs{$name[0]}]->print($split[1] . "\t" . $split[2] . "\t" . $split[3] . "\n");
		}
		foreach my $motif (keys %scan_candidates) {
			close $fileHandlesBackground[$motifs{$motif}];
		}
	}
	#Read in all files with background distribution of motifs that are in the list - start with transcription factor that was use for the cip
	$motif_genomewide = $background_folder . "/" . $genome . "_" . $ab . ".txt";
	open FH, "<$motif_genomewide";
	print STDERR "Saving chipped TF " . $motif_genomewide . "\n";
	#Save the positions of the TF used for the chip, which will be the reference for the distance of all the other TF
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if($first == 0) {
			$current_chr = $split[0];
			$first++;
		}
		if($split[0] ne $current_chr) {
			$count = 0;
		}
		$current_chr = $split[0];
		$background_saved{$current_chr}[$count] = $split[1];
		$background{$split[0]}{$split[1]} = $split[2];
		$count++;
	}
	close FH;
	$count = 0;
	#After knowing where the main TF is, calculate the distance distribution of all the other TF
	foreach my $motif (keys %index_motifs) {
		print STDERR "Processing " . $motif . "\n";
		$file = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		open FH, "<$file";
		$count = 0;
		foreach my $line (<FH>) {
			chomp $line;
			#Start with first binding motif of the reference transcription factor and run through second TF till the distance of count is smaller than the distance of count + 1 (we start running away from the motif again - so this is the closest instance) 
			@split = split('\t', $line);
			$current_dist = $split[1] - $background_saved{$split[0]}[$count];
			$next_dist = $split[1] - $background_saved{$split[0]}[$count + 1];
			while(abs($split[1] - $background_saved{$split[0]}[$count]) > abs($split[1] - $background_saved{$split[0]}[$count + 1])) {
				$count++;
			}
			$current_dist = $split[1] - $background_saved{$split[0]}[$count];
			$next_dist = $split[1] - $background_saved{$split[0]}[$count + 1];
			#Determine which motif is closer to the main TF (the one to the left (current_dist) or to the right (next_dist)
			if(abs($current_dist) < abs($next_dist)) {
				$smallest_dist = $current_dist;
				$main_tf_start = $background_saved{$split[0]}[$count];
				$main_tf_end = $background{$split[0]}{$main_tf_start};
			} else {
				$smallest_dist = $next_dist;
				$main_tf_start = $background_saved{$split[0]}[$count + 1];
				$main_tf_end = $background{$split[0]}{$main_tf_start};
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
	}
	return ($delete, \%dist_plot_background);
}

#Method that analyzes the different bases within the motif to find which bases are mutated the most and have the most impact on the TF binding
sub analyze_motif_pos{
	$_ = 0 for my($current_shift, $prev_shift, $shift_vector, $chr_num, $current_pos, $max_motif, $num_of_muts, $max);
	$_ = () for my($tree_tmp, @fc_split, @char_one, @char_two, @header, %matrix_pos_muts);
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
	my $allele = 1;
	foreach my $last_line (keys %block) {
		chomp $last_line;
		@header = split("_", $last_line);
		$chr_num = substr($header[0], 3);
		@fc_split = split('\t', $fc->{$chr_num}->{$header[1]});
		$current_pos = $header[1];
		for(my $i = 0; $i < @strains; $i++) {
			if($chr_num !~ /\d+/ && !exists $lookup->{$strains[$i]}->{$chr_num}) {
				print STDERR "Skip analysis of chromosome " . $chr_num . "in analyze_motifs_pos\n";
			}
			if(exists $lookup->{$strains[$i]}->{$chr_num}) {
				$chr_num = $lookup->{$strains[$i]}->{$chr_num};
			}
		}
		foreach my $motif (keys %{$block{$last_line}}) {
			#Define indels, snps and multiple snps (N == not significant, S == significant) 
			if(!exists $matrix_pos_muts{$motif}{'indel'}) {
				$matrix_pos_muts{$motif}{'indel'}{'N'} = 0;
				$matrix_pos_muts{$motif}{'indel'}{'S'} = 0;
			}
			if(!exists $matrix_pos_muts{$motif}{'multi'}) {
				$matrix_pos_muts{$motif}{'multi'}{'N'} = 0;
				$matrix_pos_muts{$motif}{'multi'}{'S'} = 0;
			}
			#Run through all the motif occurrences in the peak file
			foreach my $motif_pos (keys %{$block{$last_line}{$motif}}) {
				#If not binary decisions but motif score difference define significane, the strain motif with the highest score needs to be found first
				if($motif_diff > 0 || $motif_diff_percentage > 0) {
					$max_motif = 0;
					for(my $i = 0; $i < @strains; $i++) {
						if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $max_motif) {
							$max_motif = $block{$last_line}{$motif}{$motif_pos}{$strains[$i]};
						}
					}
				}
				#Run through all pairwise strain comparisons
				for(my $i = 0; $i < @strains - 1; $i++) {
					for(my $j = $i + 1; $j < @strains; $j++) {
						#Motif is different between the strains
						#motif differences and percentage of motif difference is not used to call different motifs - the score of the motif is 0 in one of the strains
						if(($motif_diff == 0 && $motif_diff_percentage == 0 && ($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} == 0 || $block{$last_line}{$motif}{$motif_pos}{$strains[$j]} == 0)) ||
						#Motif difference is used to call different motifs and the difference between the two motif scores is greater than what is defined as minmum different to be significant
						($motif_diff > 0 && abs($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} - $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) > $motif_diff) ||
						#Percentage of motif score difference is used to call different motifs and the percentual difference is bigger than the defined minimum to call it as significant
						($motif_diff_percentage > 0 && abs($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} - $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) > ($max_motif * $motif_diff_percentage))) {
							#Check if there was an indel between the start of the sequence and the position of the motif in the strain because the original sequences needs to be used in this method
							if(defined $tree->{$strains[$i]}->{$chr_num}) {
								if($current_pos + $motif_pos > $last->{$strains[$i]}->{$chr_num}->{$allele}->{'pos'}) {
									$current_shift = $last->{$strains[$i]}->{$chr_num}->{$allele}->{'shift'};
								} else {
									$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos + 1);
									$current_shift = $tree_tmp->[0]->{'shift'};
								}
								if($current_pos > $last->{$strains[$i]}->{$chr_num}->{$allele}->{'pos'}) {
									$prev_shift = $last->{$strains[$i]}->{$chr_num}->{$allele}->{'shift'};
								} else {
									$tree_tmp = $tree->{$strains[$i]}->{$chr_num}->fetch($current_pos, $current_pos + 1);
									$prev_shift = $tree_tmp->[0]->{'shift'};
								}
								$shift_vector = $current_shift - $prev_shift;	
							}
							#Split the original strain motif sequence into an array, so a base by base comparison is possible
							if($block{$last_line}{$motif}{$motif_pos}{'orientation'} eq "+") {
								@char_one = split("", substr($seq->{$last_line . "_" . $strains[$i]}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'}));
							} else {
								@char_one = split("", &rev_comp(substr($seq->{$last_line . "_" . $strains[$i]}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'})));
							}
							$shift_vector = 0;
							#Check if there was an indel between the start of the sequence and the position of the motif in the second strain
							if(defined $tree->{$strains[$j]}->{$chr_num}) {
								if($current_pos + $motif_pos > $last->{$strains[$j]}->{$chr_num}->{$allele}->{'pos'}) {
									$current_shift = $last->{$strains[$j]}->{$chr_num}->{$allele}->{'shift'};
								} else {
									$tree_tmp = $tree->{$strains[$j]}->{$chr_num}->fetch($current_pos + $motif_pos, $current_pos + $motif_pos + 1);
									$current_shift = $tree_tmp->[0]->{'shift'};
								}
								if($current_pos > $last->{$strains[$j]}->{$chr_num}->{$allele}->{'pos'}) {
									$prev_shift = $last->{$strains[$j]}->{$chr_num}->{$allele}->{'shift'};
								} else {
									$tree_tmp = $tree->{$strains[$j]}->{$chr_num}->fetch($current_pos, $current_pos + 1);
									$prev_shift = $tree_tmp->[0]->{'shift'};
								}
								$shift_vector = $current_shift - $prev_shift;	
							}
							#Split the original sequence of the second strain into an array for base by base comparison
							if($block{$last_line}{$motif}{$motif_pos}{'orientation'} eq "+") {
								@char_two = split("", substr($seq->{$last_line . "_" . $strains[$j]}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'}));
							} else {
								@char_two = split("", &rev_comp(substr($seq->{$last_line . "_" . $strains[$j]}, $motif_pos + $shift_vector, $block{$last_line}{$motif}{$motif_pos}{'length'})));
							}
							#Count number of differences - if more than 2 dismiss this comparison because motifs are not the same and count it as indel
							$num_of_muts = 0;
							for(my $c = 0; $c < @char_one; $c++) {
								if($char_one[$c] ne $char_two[$c]) {
									$num_of_muts++;
								}
							}
							if($num_of_muts > 2) {
								#Number of mutation is greater than 2 -> mutation is counted as indel - check FC in order to define if mutation changed the binding of TF
								if($delta_tag == 0) {
									if($fc_split[$i + $j - 1] > log($fc_significant)/log(2) || $fc_split[$i + $j - 1] < log(1/$fc_significant)/log(2)) {
										$matrix_pos_muts{$motif}{'indel'}{'S'}++;
									} else {
										$matrix_pos_muts{$motif}{'indel'}{'N'}++;
									}
								} else {
									if(abs($fc_split[$i + $j - 1]) > $delta_threshold) {
										$matrix_pos_muts{$motif}{'indel'}{'S'}++;
									} else {
										$matrix_pos_muts{$motif}{'indel'}{'N'}++;
									}
								}
							} elsif($num_of_muts > 1) {
							#Mnumber of mutation is greather than 1 (and smaller than 2) -> count it as multiple mutations within the same motif - check FC in order to define if mutation changed the binding of TF
								if($delta_tag == 0) {
									if($fc_split[$i + $j - 1] > log($fc_significant)/log(2) || $fc_split[$i + $j - 1] < log(1/$fc_significant)/log(2)) {
										$matrix_pos_muts{$motif}{'multi'}{'S'}++;
									} else {
										$matrix_pos_muts{$motif}{'multi'}{'N'}++;
									}
								} else {
									if(abs($fc_split[$i + $j -1]) > $delta_threshold) {
										$matrix_pos_muts{$motif}{'multi'}{'S'}++;
									} else {
										$matrix_pos_muts{$motif}{'multi'}{'N'}++;

									}
								}
							} elsif($num_of_muts == 1) {
								#Number of mutation is exactly 1 - check the kind of subsitution and determine if it was significant or not
								for(my $c = 0; $c < @char_one; $c++) {
									if($char_one[$c] ne $char_two[$c]) {
										if($delta_tag == 0) {
											if($fc_split[$i + $j - 1] > log($fc_significant)/log(2) || $fc_split[$i + $j - 1] < log(1/$fc_significant)/log(2)) {
												#S meaning significant
												if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) {
													$matrix_pos_muts{$motif}{$char_two[$c]}{$c}{'S'}++;
												} else {
													$matrix_pos_muts{$motif}{$char_one[$c]}{$c}{'S'}++;
												}
											} else {
												if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) {
													$matrix_pos_muts{$motif}{$char_two[$c]}{$c}{'N'}++;
												} else {
													$matrix_pos_muts{$motif}{$char_one[$c]}{$c}{'N'}++;
												}
											}
										} else {
											if(abs($fc_split[$i + $j - 1]) > $delta_threshold) {
												#S meaning significant
												if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) {
													$matrix_pos_muts{$motif}{$char_two[$c]}{$c}{'S'}++;
												} else {
													$matrix_pos_muts{$motif}{$char_one[$c]}{$c}{'S'}++;
												}
											} else {
												if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) {
													$matrix_pos_muts{$motif}{$char_two[$c]}{$c}{'N'}++;
												} else {
													$matrix_pos_muts{$motif}{$char_one[$c]}{$c}{'N'}++;
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
	my @base = split('', $local_seq);
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
		@clean2 = split/[\(,\/]+/, $clean1[-1];
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
		@clean1 = split(',', $name);
		if(@clean1 > 1) {
			$name = $clean1[1];
		} else {
			$name = $clean1[0];
		}
		@clean2 = split/[\(,\/]+/, $name;
		$name = $clean2[0];
	}
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
	open FH, "<$_[0]";
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
	my $command = "homer2 find -i " . $seq_file . " -m " . $tf . " -offset 0 > " . $out_file;
	`$command`;	
}

#Get strain sequences for every peak
sub get_seq_for_peaks {
	$_ = 0 for my ($general_offset, $seq_number, $newlines, $longest_seq, $no_shift, $filter_no_mut, $length, $ref_start, $working_start, $working_stop, $byte_offset, $chr_num, $shift_vector, $mut); 
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
	$line_number = $line_number * @strains;
	for(my $i = 0; $i < @strains; $i++) {
		foreach my $chr (keys %peaks) {
			$chr_num = $chr;
			$no_shift = 0;
			#Start with offset for first line
			$general_offset = 5 + length($chr);
			if(exists $lookup->{$strains[$i]}->{$chr_num}) {
				$chr_num = $lookup->{$strains[0]}->{$chr};
			}
			if($chr !~ /\d+/ && !exists $lookup->{$strains[$i]}->{$chr}) {
				print STDERR "Skip peaks for chromosome " . $chr . " in get_seq_for_peaks\n";
				next;
			}
			$current_tree = $tree->{$strains[$i]}->{$chr_num};
			if(!defined $current_tree) {
				$no_shift = 1;
			}
			$filename = $data . "/" . uc($strains[$i]) . "/chr" . $chr . "_allele_" . $allele . ".fa";
			if(!-e $filename) {
				print STDERR "Can't open $filename\n";
				print STDERR "Exclude chromosome $chr from analysis\n";
				delete $peaks{$chr};
				#Also remove it from sequence 
				foreach my $s (%seq) {
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
					$ref_start = $current_tree->fetch($working_start, $working_start + 1);
					if(!exists $ref_start->[0]->{'shift'}) {
						$shift_vector = $last->{$strains[$i]}->{$chr}->{$allele}->{'shift'};
					} else {
						$shift_vector = $ref_start->[0]->{'shift'};
					}
					$working_start = $working_start + $shift_vector;
					$ref_start = $current_tree->fetch($working_stop, $working_stop + 1);
					if(!exists $ref_start->[0]->{'shift'}) {
						$shift_vector = $last->{$strains[$i]}->{$chr}->{$allele}->{'shift'};
					} else {
						$shift_vector = $ref_start->[0]->{'shift'};
					}
					$working_stop = $working_stop + $shift_vector;
				}
				$header = "chr" . $chr . "_" . $start_pos . "_" . $peaks{$chr}{$start_pos} . "_" . $strains[$i];
				#Calculate length and newlines (fastq file saves seq in 50bp lines) 
				$length = $working_stop - $working_start;
				$newlines = int($working_stop/50) - int($working_start/50);
				#Add number of newlines to length (because working on bytes)
				$length = $length + $newlines;
				#Calculate the number of newlines before before working start in sequence
				$newlines = int(($working_start)/50);
				#Calculate offset - general offset is first line
				$byte_offset = $working_start + $general_offset + $newlines;
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
	#Only peaks with mutations are considered in this analysis
	if($mut_only == 1) {
		print STDERR "Filter out peaks without mutations\n";
		foreach my $key (keys %seq) {
			@split_part = split("_", $key);
			$h = $split_part[0] . "_" . $split_part[1] . "_" . $split_part[2];
			if(exists $seen_part{$h}) { next; }
			$seq_first = $seq{$h . "_" . $strains[0]};
			#Compare all strains to first strain - if one sequence is different keep sequences
			for(my $i = 1; $i < @strains; $i++) {
				if($seq_first ne $seq{$h . "_" . $strains[$i]}) {
					$mut = 1;
				} 
			}
			if($mut == 0) {
				$filter_no_mut++;
				#No mutation was found - delete this peak
				for(my $i = 0; $i < @strains; $i++) {
					delete $seq{$h . "_" . $strains[$i]};
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


sub all_vs_all_comparison{
	my $block = $_[0];
	my $recenter_conversion = $_[1];
	my $tag_counts = $_[2];
	my @strains = @{$_[3]};
	my $method = $_[4];
	my $same = 0;
	my $value = 0;
	my $all = keys %{$block};
	my $count = 0;
	$_ = () for my(%correlation, @sum, @tag_sum, $vector_motifs, $vector_tagcounts, @split, $chr, $stddev_m, $stddev_t, $cor, $t, $rank, %without_motif, %with_motif, @no_motif, @motif, $ttest, $max, $max_motif, $min_motif, $max_no_motif, $min_no_motif, @all_motifs, @all_tags, %save_motifs, %save_tags, %save_no_motif, %save_motif, @all_motif, @all_no_motif, %group_with, %group_without, $R);
	my($con_pos, $r, $r2, $distance, $h_tag, $h_sum, $information, $max);
	if($method eq "group" || $method eq "group_all" || $method eq "group_all_scale") {
#		$ttest = new Statistics::TTest;  
#		$ttest->set_significance(90);
		$R = Statistics::R->new();
	}
	foreach my $pos (sort {$a cmp $b} keys %{$block}) {
		print STDERR "Processing: " .  ($count/$all)*100 . "\r";
		$count++;
		@split = split("_", $pos);
		$chr = substr($split[0], 3);
 		$con_pos = $split[1];
		if(exists $recenter_conversion->{$chr}->{$con_pos}) {
			$con_pos = $recenter_conversion->{$chr}->{$con_pos}->{'start'};
		}
		@tag_sum = split('\s+', $tag_counts->{$chr}->{$con_pos});
		foreach my $motif (sort {$a cmp $b} keys %{$block->{$pos}}) {
			@motif = ();
			@no_motif = ();
			%without_motif = ();
			my %tmp_save;
			@sum = ();
			#Sum over motif score for this strain and this peak when there are multiple motifs of the same TF
			foreach my $motif_pos (sort {$a <=> $b} keys %{$block->{$pos}->{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					$sum[$i] += $block->{$pos}->{$motif}->{$motif_pos}->{$strains[$i]};
					$tmp_save{$strains[$i]} = $tag_sum[$i];
					if($block->{$pos}->{$motif}->{$motif_pos}->{$strains[$i]} < 5) {
						$without_motif{$strains[$i]} = 1;
					}
				}
			}
			if($method eq "group_all_scale") {
				$max = 0;
				for(my $i = 0; $i < @strains; $i++) {
					if($max < $tmp_save{$strains[$i]}) { $max = $tmp_save{$strains[$i]}; }
				}
				for(my $i = 0; $i < @strains; $i++) {
					$tmp_save{$strains[$i]} = $tmp_save{$strains[$i]}/$max;
				}
				
			}		

			if($method eq "group_all" || $method eq "group_all_scale") {
				for(my $i = 0; $i < @strains; $i++) {
					if(exists $without_motif{$strains[$i]}) {
						$group_without{$motif}{$tmp_save{$strains[$i]}}++;
					} else {
						$group_with{$motif}{$tmp_save{$strains[$i]}}++;
					}
				}
			}
			$value = $sum[0];
			$same = 0;
			for(my $i = 1; $i < @sum; $i++) {
				if($value != $sum[$i]) {
					$same = 1;
					last;
				}
			}
			if($same == 0) { next; }
			#Convert into a vector - vector of motif scores per strain
			$vector_motifs = vector(@sum);
			#Convert tag counts into a vector - vector of tag counts per strain
			$vector_tagcounts = vector(@tag_sum);
			if($method eq "pearson_all" || $method eq "spearman_all" || $method eq "pearson_all_scale" || $method eq "spearman_all_scale") {
				if(!exists $save_motifs{$motif}) {
					@{$save_motifs{$motif}} = ();
					@{$save_tags{$motif}} = ();
				}
				@all_motifs = @{$save_motifs{$motif}};
				@all_tags = @{$save_tags{$motif}};
				if($method eq "pearson_all_scale") {
					$max = 0;
					for(my $i = 0; $i < @tag_sum; $i++) {
						if($max < $tag_sum[$i]) { $max = $tag_sum[$i]; }
					}
					for(my $i = 0; $i < @tag_sum; $i++) {
						$tag_sum[$i] = $tag_sum[$i]/$max;
					}
				}
				push @all_tags, @tag_sum;
				push @all_motifs, @sum;
				@{$save_tags{$motif}} = @all_tags;
				@{$save_motifs{$motif}} = @all_motifs;
			} elsif($method eq "spearman") {
				$rank = Statistics::RankCorrelation->new(\@sum, \@tag_sum);
				$cor = $rank->spearman;
				$correlation{$motif}{$cor}++;
			} elsif($method eq "pearson") {
				#Calculate standard deviation for motifs and tag counts
				$stddev_m = stddev( $vector_motifs) * 1;
				$stddev_t = stddev( $vector_tagcounts) * 1;
				#Set SD to 0 if very small
				if($stddev_m < 0.00001) { $stddev_m = 0; }
				if($stddev_t < 0.00001) { $stddev_t = 0; }
				#Calcualte pearson correlation between motif vector and tag count vector
				$cor = correlation( $vector_motifs, $vector_tagcounts );
				#If correlation is na - save it
				#If correalation is na either the motif score vector or the tag count vector consisted of identical values so there was no information in it
				if($cor ne "n/a" && $stddev_m > 0 && $stddev_t > 0) {
					$correlation{$motif}{$cor}++;
				}
			} elsif($method eq "group") {
				$max_motif = 0;
				$min_motif = 10000;
				$max_no_motif = 0;
				$min_no_motif = 10000;
				for(my $i = 0; $i < @strains; $i++) {
					if(exists $without_motif{$strains[$i]}) {
						push @no_motif, $tag_sum[$i];
						if($tag_sum[$i] > $max_no_motif) { $max_no_motif = $tag_sum[$i]; }
						if($tag_sum[$i] < $min_no_motif) { $min_no_motif = $tag_sum[$i]; }
					} else {
						push @motif, $tag_sum[$i];
						if($tag_sum[$i] > $max_motif) { $max_motif = $tag_sum[$i]; }
						if($tag_sum[$i] < $min_motif) { $min_motif = $tag_sum[$i]; }
					}
				}
				if(@motif > 0 && @no_motif > 0) {
					if(@motif < 2 || @no_motif < 2) {
						$max = @motif;
						for(my $i = 0; $i < $max; $i++) {
							push @motif, $motif[$i] + rand($motif[$i] * 0.1);
						}
						$max = @no_motif;
						for(my $i = 0; $i < $max; $i++) {
							push @no_motif, $no_motif[$i] + rand($no_motif[$i] * 0.1);
						}
					} else {
						if($max_motif == $min_motif || $max_no_motif == $min_no_motif) {
							for(my $i = 0; $i < @motif; $i++) {
								$motif[$i] = $motif[$i] + rand($motif[$i] * 0.1);
							}
							for(my $i = 0; $i < @no_motif; $i++) {
								$no_motif[$i] = $no_motif[$i] + rand($no_motif[$i] * 0.1);
							}
						}
					}
					$R->set( 'x', \@motif );
					$R->set( 'y', \@no_motif);
					$R->run( q`res = t.test(x,y)` );
					$ttest = $R->get('res');
					$correlation{$motif}{$ttest->[16]}++;
				}
			} elsif($method eq "group_all" || $method eq "group_all_scale" || $method eq "pearson_all" || $method eq "pearson_all_scale" || $method eq "spearman_all" || $method eq "spearman_all_scale") {
			} elsif($method eq "mutual") {
				#Calculate euclidean distance between vectors
				my $distance = &euclidean_distance(\@tag_sum, \@sum);
				my $h_tag = &shannon_entropy(\@tag_sum);
				my $h_sum = &shannon_entropy(\@sum);
				#Calculate mutual information
				my $information = ($h_tag + $h_sum - $distance)/2;
				$correlation{$motif}{$information}++;
			} else {
				print STDERR "Wrong method!\n";
				exit;
			}
		}
	}
	if($method eq "pearson_all" || $method eq "pearson_all_scale") {
		foreach my $motif ( keys %save_tags) {
			$vector_motifs = vector @{$save_motifs{$motif}};
			$vector_tagcounts = vector @{$save_tags{$motif}};
			$cor = correlation( $vector_motifs, $vector_tagcounts );
			$correlation{$motif}{$cor} = 1;
		}	
	}
	if($method eq "spearman_all" || $method eq "spearman_all_scale") {
		foreach my $motif(keys %save_tags) {
			$rank = Statistics::RankCorrelation->new(\@{$save_motifs{$motif}}, \@{$save_tags{$motif}});
			$cor = $rank->spearman;
			$correlation{$motif}{$cor}++;
		}
	}
	if($method eq "group_all" || $method eq "group_all_scale") {
		foreach my $motif (sort {$a cmp $b } keys %group_with) {
			@motif = ();
			@no_motif = ();
			foreach my $v (sort {$a <=> $b} keys %{$group_with{$motif}}) {
				if($v eq "") { next; }
				push @motif, $v foreach (1..$group_with{$motif}{$v});
			}
			foreach my $v (sort {$a <=> $b} keys %{$group_without{$motif}}) {
				if($v eq "") { next; }
				push @no_motif, $v foreach (1..$group_without{$motif}{$v});
			}
			if(@motif < 1 || @no_motif < 5) { next; }
			$R->set( 'x', \@motif);
			$R->set( 'y', \@no_motif);
			$R->run( q`res = t.test(x, y)` );
			$ttest = $R->get('res');
			$correlation{$motif}{$ttest->[16]}++;
		}
	}
	return (\%correlation);
}

sub shannon_entropy{
	my $vector = $_[0];
	my $all = scalar @{$vector};
	my $px = 0;
	my $H = 0;
	for(my $i = 0; $i < @{$vector}; $i++) {
		$px = (($vector->[$i]+1)/$all);
		$H -= $px * (log($px)/log(2));
	}
	return $H;
	
}

sub euclidean_distance{
	my $vector1 = $_[0];
	my $vector2 = $_[1];
	my $distance = 0;
	for(my $i = 0; $i < @{$vector1}; $i++) {
		$distance += ($vector1->[$i] - $vector2->[$i])**2;
	}
	$distance = sqrt($distance);
	return $distance;
}

