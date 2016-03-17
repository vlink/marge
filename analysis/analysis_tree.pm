#!/usr/bin/perl

BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
package analysis;
use strict;
use Getopt::Long;
#require '/Users/verenalink/workspace/strains/general/config.pm';
use config;
use Data::Dumper;
use Statistics::Basic qw(:all);
use Statistics::Distributions;
#require '../db_part/database_interaction.pm';

$_ = "" for my($genome, $file, $tf, $filename, $last_line);
$_ = () for my(@strains, %promoter, %exp, @split, @name, %mutation_pos, %shift, %current_pos, %save_local_shift, %seq, %PWM, @fileHandlesMotif, %index_motifs, %tag_counts, %fc, @header, %block, %existance, %diff, %ranked_order, %mut_one, %mut_two, %delta_score, %matrix_pos_muts, %network_save, %cytoscape_sig, %cytoscape_unsig);
$_ = 0 for my($homo, $allele, $region, $motif_score, $motif_start, $more_motifs, $save_pos, $delta, $tmp_out, $tmp_out_main_motif, $line_number);
my $data = config::read_config()->{'data_folder'};

my %comp = ();
$comp{'A'} = 'T';
$comp{'T'} = 'A';
$comp{'C'} = 'G';
$comp{'G'} = 'C';


print STDERR "Matrix is defined in read_motifs and it is storred globally in this module!\n";
print STDERR "Change that!!!!!\n\n";


sub write_header{
	my @fileHandlesMotif = @{$_[0]};
	my @strains = @{$_[1]};
	my $fc = $_[2];
	for(my $files = 0; $files < @fileHandlesMotif; $files++) {
		$fileHandlesMotif[$files]->print("motif\tpos\t");
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$files]->print("tg " . $strains[$i] . "\t");
		}
		if($fc == 0) {
			for(my $i = 0; $i < @strains - 1; $i++) {
				for(my $j = $i + 1; $j < @strains; $j++) {
					$fileHandlesMotif[$files]->print("fc ". $strains[$i] . " vs " . $strains[$j] . "\t");
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


sub analyze_motifs{
	my @split;
	my @header;
	my $motif_file = $_[0];
	my $seq = $_[1];
	my @strains = @{$_[2]};
	my $tree = $_[3];
	my $lookup = $_[4];
	my $last = $_[5];
	my $allele = $_[6];
	my $region = $_[7];
        open FH, "<$motif_file";
        my $last_line = "";
	my %block = ();
	my $pos_beginning;
	my $shift_beginning;
	my $pos_current;
	my $pos_end;
	my $shift_current;
	my $shift_diff;
	my $shift_end;
	my $tree_tmp;
	my $motif_pos;
	my $chr_num;
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
		if($split[4] eq "-") {
			$motif_start = $motif_start - length($split[2]) + 1;
			$motif_pos = $motif_pos - length($split[2]) + 1;
		}
		$pos_beginning = $header[1] - ($region/2);
		$chr_num = substr($header[0], 3);
		if($chr_num !~ /\d+/ && !exists $lookup->{$header[3]}->{$chr_num}) {
			print STDERR "Skip analysis of chromosome " . $chr_num;
		}
		if(exists $lookup->{$header[3]}->{$chr_num}) {
			$chr_num = $lookup->{$header[3]}->{$chr_num};
		}
		if(defined $tree->{$header[3]}->{$chr_num}) {
			$tree_tmp = $tree->{$header[3]}->{$chr_num}->fetch($pos_beginning, $pos_beginning + 1);
			if(!exists $tree_tmp->[0]->{'shift'}) {
				$shift_beginning = $last->{$header[3]}->{$header[0]}->{$allele}->{'shift'};
			} else {
				$shift_beginning = $tree_tmp->[0]->{'shift'};
			}
			$tree_tmp = $tree->{$header[3]}->{$chr_num}->fetch($motif_start, $motif_start + 1);
			if(!exists $tree_tmp->[0]->{'shift'}) {
				$shift_current = $last->{$header[3]}->{$header[0]}->{$allele}->{'shift'};
			} else {
				$shift_current = $tree_tmp->[0]->{'shift'};
			}
			$shift_diff = $shift_beginning - $shift_current;
			$motif_pos = $motif_pos + $shift_diff;
		}
		$block{$last_line}{$split[3]}{$motif_pos}{$header[-1]} = $split[-1];
		$block{$last_line}{$split[3]}{$motif_pos}{'length'} = length($split[2]);
        }
	return \%block;
}

sub save_all_motifs{
	my $block_ref = $_[0];
	my $last_line = $_[1];
	my %block = %$block_ref;
	foreach my $motif (keys %block) {
		foreach my $pos (keys %{$block{$motif}}) {
			$network_save{$last_line}{$pos}{$motif} = $block{$motif}{$pos};
		}
	}
}

sub merge_block {
	my $block_ref = $_[0];
	my %block = %$block_ref;
	my @strains = @{$_[2]};
	my $overlap = $_[1];
	my $seq = $_[3];
	my $ab_mut = 0;
	my $num_of_bp = 0;
        #Now run through all motifs and see if they are in all the other strains
        my %ignore = ();
	foreach my $chr_pos (keys %block) {
		foreach my $motif (keys %{$block{$chr_pos}}) {
			%existance = ();
			%diff = ();
			foreach my $motif_pos (keys %{$block{$chr_pos}{$motif}}) {
				if($overlap eq "half") {
					$num_of_bp = int($block{$chr_pos}{$motif}{$motif_pos}{'length'}/2)
				} elsif($overlap eq "complete") {
					$num_of_bp = 0;
				} else {
					$num_of_bp = ($overlap*1);
				}
				if(exists $ignore{$motif_pos}) {
					delete $block{$chr_pos}{$motif}{$motif_pos};
					next;
				}
				$save_pos = $motif_pos;
				$motif_score = 0;
				for(my $i = 0; $i < @strains; $i++) {
					if(!exists $diff{$strains[$i]}) {
						$diff{$strains[$i]} = 0;
					}
					if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}) {
						#Check in vicinity (+/- motif_lenght-1) - if there is a motif within motif/2 -> save the current motif under this position, to not count them as separate motifs
						$existance{$strains[$i]} = 1;
						for(my $run_motif = $motif_pos - $num_of_bp; $run_motif < $motif_pos + $num_of_bp; $run_motif++) {
							if(exists $block{$chr_pos}{$motif}{$run_motif} && exists $block{$chr_pos}{$motif}{$run_motif}{$strains[$i]} && $run_motif != $motif_pos) {
								$block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} = $block{$chr_pos}{$motif}{$run_motif}{$strains[$i]};
								$diff{$strains[$i]} += $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]};
								delete $existance{$strains[$i]};
								delete $block{$chr_pos}{$motif}{$run_motif}{$strains[$i]};
								$ignore{$run_motif} = 1;
								if(exists $block{$chr_pos}{$motif}{$run_motif}{'length'}) {
									delete $block{$chr_pos}{$motif}{$run_motif}{'length'};
								}
							}
						}
					}
					if($block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} eq "" || !exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]}) {
						$block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} = &calculate_motif_score($motif, substr($seq->{$chr_pos . "_" . $strains[$i]}, $motif_pos + 1, $block{$chr_pos}{$motif}{$motif_pos}{'length'}));
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

sub output_motifs{
	my $block_ref = $_[0];
	my %block = %$block_ref;
	my @fileHandlesMotif = @{$_[1]};
	my %tag_counts = %{$_[2]};
	my @strains = @{$_[3]};
	my %index_motifs = %{$_[4]};
	my $fc_significant = $_[5];
	my $overlap = $_[6];
	my %fc = %{$_[7]};
	my $fc_exists = 0;
	if((keys %fc) > 1) {
		$fc_exists = 1;
	}
	my %recenter_conversion = %{$_[8]};
	my $curr_chr;
	my $curr_pos;
	my $ab;
	my $ab_mut;
	my %diff;
	my %existance;
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
					if(!exists $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} || $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]} == 0) {
						$existance{$strains[$i]} = 1;
					}
					$diff{$strains[$i]} = $diff{$strains[$i]} + $block{$chr_pos}{$motif}{$motif_pos}{$strains[$i]};
				}
			}
			if($fc_exists == 1) {
				$fileHandlesMotif[$index_motifs{$motif}]->print($fc{$curr_chr}{$curr_pos});
			}
			for(my $i = 0; $i < @strains; $i++) {
				$fileHandlesMotif[$index_motifs{$motif}]->print($diff{$strains[$i]} . "\t");
			}
			for(my $i = 0; $i < @strains; $i++) {
				if(exists $existance{$strains[$i]}) {
					$fileHandlesMotif[$index_motifs{$motif}]->print("1\t");
					if($motif eq $ab) {
						$ab_mut = 1;
					}
				} else {
					$fileHandlesMotif[$index_motifs{$motif}]->print("0\t");
				}
			}
			$fileHandlesMotif[$index_motifs{$motif}]->print("" . (keys %{$block{$chr_pos}{$motif}}) . "\n");
	#		if($mut_pos == 1 && $fc_exists == 1) {
	#			&analyze_motif_pos(\%block, $last_line, $fc{$curr_chr}{$curr_pos}, \@strains, $fc_significant);
	#		}
		}
	}
}


sub distance_plot{
	my %block = %{$_[0]};
	my @strains = @{$_[1]};
	my %fc = %{$_[2]};
	my $fc_significant = $_[3];
	my $effect = $_[4];
	my $longest_seq = $_[5];
	my $seq = $_[6];
	my %dist_plot;
	my $offset;
	foreach my $last_line (keys %block) {
		@split = split("_", $last_line);
		if($effect == 1 && ($fc{substr($split[0], 3)}{$split[1]} < log($fc_significant)/log(2) && $fc{substr($split[0], 3)}{$split[1]} > log(1/$fc_significant)/log(2))) {
			next;
		}
		$offset = int(($longest_seq - length($seq->{$last_line . "_" . $strains[0]}))/2);
		foreach my $motif (keys %{$block{$last_line}}) {
			foreach my $pos (keys %{$block{$last_line}{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					if(!exists $block{$last_line}{$motif}{$pos}{$strains[$i]}) {
						next;
					}
					for(my $l = 0; $l < $block{$last_line}{$motif}{$pos}{'length'}; $l++) {
					#	$dist_plot{$motif}{$strains[$i]}{($pos + $l + $offset - int($longest_seq/2))}++;
						$dist_plot{$motif}{$strains[$i]}{($pos + $l + $offset)}++;
					}
				}
			}	
		}
	}
	return \%dist_plot;
}

sub background_dist_plot{
	my $background_folder = $_[0];	
	my %index_motifs = %{$_[1]};
	my $delete = $_[2];
	my $motif_scan_score = $_[3];
	my $genome = $_[4];
	my $ab = $_[5];
	my $region = $_[6];
	my $longest_seq = $_[7];
	my $motif_genomewide;
	my %scan_candidates;
	my %background_dist;	

	print STDERR "Checking for genome wide scan of motifs\n";
	foreach my $motif (keys %index_motifs) {
		$motif_genomewide = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		if(!-e $motif_genomewide) {
			$scan_candidates{$motif} = 1;
		}
	}


	#Scan for motifs that were not already proprocessed
	if(keys %scan_candidates > 0) {
		my $tmp_motif = "tmp" . rand(15) . ".txt";
		$delete->{$tmp_motif} = 1;
		open TMP, ">$tmp_motif";
		foreach my $motif (keys %scan_candidates) {
			print $motif . "\n";
			print TMP ">consensus_" . $motif . "\t$motif\t" . $motif_scan_score->{$motif} . "\n";
			foreach my $pos (sort {$a <=> $b} keys %{$PWM{$motif}}) {
				print $pos . "\n";
				print TMP $PWM{$motif}{$pos}{'A'} . "\t" . $PWM{$motif}{$pos}{'C'} . "\t" . $PWM{$motif}{$pos}{'G'} . "\t" . $PWM{$motif}{$pos}{'T'} . "\n";
			}
		}
		close TMP;

		print STDERR "Scan motifs genome wide\n";
		my $tmp_motif2 = "tmp" . rand(15);
		$delete->{$tmp_motif2} = 1;
		my $command = "scanMotifGenomeWide.pl " . $tmp_motif . " " . $genome . " > " . $tmp_motif2;
		`$command`;
		open FH, "<$tmp_motif2";

		#Save motifs in file
		my %motifs;
		my $i = 0;
		my @fileHandlesBackground;
		foreach my $motif (keys %scan_candidates) {
			my $filename = $background_folder . "/" . $genome . "_" . $motif . ".txt";
			open my $fh, ">", $filename or die "Can't open $filename: $!\n";
			$motifs{$motif} = $i;
			$fileHandlesBackground[$motifs{$motif}] = $fh;
			$i++;
		}
		my @name;
		foreach my $line (<FH>) {
			print $line;
			chomp $line;
			@split = split('\t', $line);
			@name = split('-', $split[0]);
			$fileHandlesBackground[$motifs{$name[0]}]->print($split[1] . "\t" . $split[2] . "\t" . $split[3] . "\n");
		}
		foreach my $motif (keys %scan_candidates) {
			close $fileHandlesBackground[$motifs{$motif}];
		}
	}

	#Now read in all file - we start with the main transcription factor
	$motif_genomewide = $background_folder . "/" . $genome . "_" . $ab . ".txt";
	open FH, "<$motif_genomewide";
	my %background;
	my %background_saved;
	my @background;
	my $current_chr = "";
	my $first = 0;
	my $count = 0;
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if($first == 0) {
			$current_chr = $split[0];
			$first++;
		}
		if($split[0] ne $current_chr) {
			$background_saved{$current_chr} = \@background;
			@background = ();
			$count = 0;
		}
		$current_chr = $split[0];
		$background[$count] = $split[1];
		$background{$split[0]}{$split[1]} = $split[2];
		$count++;
	}
	close FH;
	$background_saved{$current_chr} = \@background;
	my $last_index = 0;
	$count = 0;
	#Time to start calculation background distribution
	my $current_dist = 0;
	my $next_dist = 0;
	my $smallest_dist = 0;
	my %dist_plot_background;
	my $main_tf_start;
	my $main_tf_end;
	my $length = 0;
	foreach my $motif (keys %index_motifs) {
		my $file = $background_folder . "/" . $genome . "_" . $motif . ".txt";
		open FH, "<$file";
		$count = 0;
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$current_dist = $split[1] - $background_saved{$split[0]}[$count];
			$next_dist = $split[1] - $background_saved{$split[0]}[$count + 1];
			
			while(abs($split[1] - $background_saved{$split[0]}[$count]) > abs($split[1] - $background_saved{$split[0]}[$count + 1])) {
				$count++;
			}
			$current_dist = $split[1] - $background_saved{$split[0]}[$count];
			$next_dist = $split[1] - $background_saved{$split[0]}[$count + 1];
			if(abs($current_dist) < abs($next_dist)) {
				$smallest_dist = $current_dist;
				$main_tf_start = $background[$count];
				$main_tf_end = $background{$split[0]}{$main_tf_start};
			} else {
				$smallest_dist = $next_dist;
				$main_tf_start = $background[$count + 1];
				$main_tf_end = $background{$split[0]}{$main_tf_start};
			}
			if(abs($smallest_dist) < int($longest_seq/2) + $region) {
				$length = ($main_tf_end - $main_tf_start);
				for(my $i = int($length/2) * -1; $i < int($length/2); $i++) {
					$dist_plot_background{$motif}{$smallest_dist + $i}++;
				}
			}
		}
	}
	return ($delete, \%dist_plot_background);
}


sub analyze_motif_pos{
	my %block = %{$_[0]};
	my %fc = %{$_[1]};
	my @strains = @{$_[2]};
	my $fc_significant = $_[3];
	my @fc_split;
	my @char_one;
	my @char_two;
	my $num_of_muts = 0;
	foreach my $last_line (keys %block) {
		@fc_split = split('\t', $fc{$last_line});
		foreach my $motif (keys %{$block{$last_line}}) {
			foreach my $motif_pos (keys %{$block{$last_line}{$motif}}) {
				for(my $i = 0; $i < @strains - 1; $i++) {
					for(my $j = $i + 1; $j < @strains; $j++) {
						if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} != $block{$last_line}{$motif}{$motif_pos}{$strains[$j]}) {
							@char_one = split("", substr($seq{$last_line . "_" . $strains[$i]}, $motif_pos, $block{$last_line}{$motif}{$motif_pos}{'length'}));
							@char_two = split("", substr($seq{$last_line . "_" . $strains[$j]}, $motif_pos, $block{$last_line}{$motif}{$motif_pos}{'length'}));
							#Count number of differences - if more than 2 dismiss this comparison because motifs are not the same
							$num_of_muts = 0;
							for(my $c = 0; $c < @char_one; $c++) {
								if($char_one[$c] ne $char_two[$c]) {
									$num_of_muts++;
								}
							}
							if($num_of_muts > 2) { next; }
							for(my $c = 0; $c < @char_one; $c++) {
								if($char_one[$c] ne $char_two[$c]) {
									if($fc_split[$i + $j - 1] > log($fc_significant)/log(2) || $fc_split[$i + $j - 1] < log(1/$fc_significant)/log(2)) {
										#S meaning significant
										if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $block{$motif}{$motif_pos}{$strains[$j]}) {
											$matrix_pos_muts{$motif}{$char_two[$c]}{$c}{'S'}++;
										} else {
											$matrix_pos_muts{$motif}{$char_one[$c]}{$c}{'S'}++;
										}
									} else {
										if($block{$last_line}{$motif}{$motif_pos}{$strains[$i]} > $block{$motif}{$motif_pos}{$strains[$j]}) {
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
	return \%matrix_pos_muts;
}
sub log_fc {
	my $local_delta = $_[0];
	if($local_delta < 0) {
		$local_delta = $local_delta * (-1);
		$local_delta = log($local_delta+1)/log(2);
		$local_delta = $local_delta * (-1);
	} else {
		$local_delta = log($local_delta+1)/log(2);
	}
	return $local_delta;
}

sub calculate_motif_score{
	my $score_forward = 0;
	my $score_reverse = 0;
	my $motif = $_[0];
	my $local_seq = $_[1];
	my @base = split('', $local_seq);
	for(my $b = 0; $b < @base; $b++) {
		$score_forward += log($PWM{$motif}{$b}{$base[$b]}/0.25); 
	}
	$local_seq = &rev_comp($local_seq);	
	@base = split('', $local_seq);
	for(my $b = 0; $b < @base; $b++) {
		$score_reverse += log($PWM{$motif}{$b}{$base[$b]}/0.25); 
	}
	if($score_reverse < 0) { $score_reverse = 0;}
	if($score_forward < 0) { $score_forward = 0;}
	if($score_reverse > $score_forward) {
		return sprintf("%.6f", $score_reverse);	
	} else {
		return sprintf("%.6f", $score_forward);
	}
}

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

sub get_motif_name{
        my $name = $_[0];
	#First check if it comes from homer motif analysis or database
	if(index($name, "BestGuess") != -1) {
		#Get rid of consensus sequence at the beginning
		my @clean1 = split("BestGuess:", $name);
		#Get rid of best guess
		#Get rif of everything in paraenthesis 
		my @clean2 = split/[\(,\/]+/, $clean1[-1];
		$name = $clean2[0];
		my @underscore = split('_', $clean2[0]);
		for(my $i = @underscore - 1; $i > 0; $i--) {
			if(length($underscore[$i]) > 2 && $underscore[$i] =~ /[a-zA-Z]+/) {
				$name = $underscore[$i];
				return $name;
				last;
			}
		}
		my @paren = split('\(', $name);
		$name = $paren[0];
	} else {
		my @clean1 = split(',', $name);
		if(@clean1 > 1) {
			$name = $clean1[1];
		} else {
			$name = $clean1[0];
		}
		my @clean2 = split/[\(,\/]+/, $name;
		$name = $clean2[0];
	}
	return $name;
}

sub open_filehandles{
	my %motifs = %{$_[0]};
	my %del = %{$_[1]};
	my $output_mut = $_[2];
	my @fileHandlesMotif;
	my $index = 0;
	foreach my $motif (keys %motifs) {
		my $filename = $output_mut . "_" . $motif . ".txt";
		$del{$filename} = 1;
		open my $fh, ">", $filename or die "Can't open $filename: $!\n";
		$fileHandlesMotif[$motifs{$motif}] = $fh;
	}
	return (\@fileHandlesMotif, \%del);
}

sub read_motifs{
	open FH, "<$_[0]";
	my $pos = 0;
	my $motif = "";
	my $index = 0;
	my $name;
	my %index_motifs;
#	my @fileHandlesMotif;
	my %score;
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			$motif = (split('\t', $line))[1];
			$pos = 0;
			#Open file per motif
			$name = &get_motif_name($motif);
#			my $filename = $output_mut . "_" . $name . ".txt";
	#		$del{$filename} = 1;
#			open my $fh, ">", $filename or die "Can't open $filename: $!\n";
			if(!exists $index_motifs{$name}) {
				$index_motifs{$name} = $index;
	#			$fileHandlesMotif[$index] = $fh;
				$index++;
				@split = split('\t', $line);
				$score{$name} = $split[2];
			}
		} else {
			@split = split('\t', $line);
			if($split[0] == 0) { $split[0] += 0.001; }
			if($split[1] == 0) { $split[1] += 0.001; }
			if($split[2] == 0) { $split[2] += 0.001; }
			if($split[3] == 0) { $split[3] += 0.001; }
			$PWM{$name}{$pos}{'A'} = $split[0];
			$PWM{$name}{$pos}{'C'} = $split[1];
			$PWM{$name}{$pos}{'G'} = $split[2];
			$PWM{$name}{$pos}{'T'} = $split[3];
			$matrix_pos_muts{$name}{$pos}{'S'} = 0;
			$matrix_pos_muts{$name}{$pos}{'N'} = 0;
			$pos++;
		}
	}
#	return (\%index_motifs, \%PWM, \@fileHandlesMotif, \%del, \%score);
	return (\%index_motifs, \%PWM, \%score);
}

sub scan_motif_with_homer{
	my $seq_file = $_[0];
	my $out_file = $_[1];
	print STDERR "out file: " . $out_file . "\n";
	my $tf = $_[2];
	my $command = "homer2 find -i " . $seq_file . " -m " . $tf . " -offset 0 > " . $out_file;
	print STDERR $command . "\n";
	`$command`;	
}

sub get_seq_for_peaks {
	open OUT, ">$_[0]";
	print STDERR $_[0] . "\n";
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
	my %current_pos;
	my $equal = 0;
	my $ref_seq = "";
	my $longest_seq = 0;
	my $h;
	my $ref_start;
	my $ref_stop;
	my $stop_pos;
	my $working_start;
	my $working_stop;
	my $current_tree = Set::IntervalTree->new;
	my $header;
	my $byte_offset;
	my $newlines;
	my $chr_num;
	my $shift_vector;
	my $no_shift = 0;
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
				print STDERR "Skip peaks for chromosome " . $chr . "\n";
				next;
			}
			$current_tree = $tree->{$strains[$i]}->{$chr_num};
			if(!defined $current_tree) {
				$no_shift = 1;
			}
			$filename = $data . "/" . uc($strains[$i]) . "/chr" . $chr . "_allele_" . $allele . ".fa";
                        open my $fh, "<", $filename or die "Can't open $filename: $!\n";
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
				$length = $length + $newlines;
				$newlines = int(($working_start)/50);
				#Calculate offset - general offset is first line
				$byte_offset = $working_start + $general_offset + $newlines;
				seek($fileHandles[0], $byte_offset, 0);
				read $fileHandles[0], $seq, $length;
				$seq =~ s/\n//g;
				if(length($seq) > $longest_seq) {
					$longest_seq = length($seq);
				}
				$current_pos{$strains[$i]} = uc($seq);
				$seq{$header} = uc($seq);
				$seq = "";
				$seq_number++;
				print STDERR "" . $seq_number . " of " . $line_number . " gathered\r";
			}
		}
	}
	print STDERR "\n";
	foreach my $key (sort {$a cmp $b} keys %seq) {
		print OUT ">" . $key . "\n";
		print OUT $seq{$key}  . "\n";
	}
	close OUT;
	return(\%seq, $longest_seq);
}

sub all_vs_all_comparison{
	my $block = $_[0];
	my $recenter_conversion = $_[1];
	my $tag_counts = $_[2];
	my @strains = @{$_[3]};
	$_ = () for my(@split, %correlation, @sum, @tag_sum, $vector_motifs, $vector_tagcounts, @split, %correlation, $chr, %pvalue, $p, $stddev_m, $stddev_t, $cor, $t);
	my($con_pos, $r, $r2, $p);
	foreach my $pos (keys %{$block}) {
		@split = split("_", $pos);
		$chr = substr($split[0], 3);
 		$con_pos = $split[1];
		if(exists $recenter_conversion->{$chr}->{$con_pos}) {
			$con_pos = $recenter_conversion->{$chr}->{$con_pos}->{'start'};
		}
		foreach my $motif (keys %{$block->{$pos}}) {
			@sum = ();
			foreach my $motif_pos (keys %{$block->{$pos}->{$motif}}) {
				for(my $i = 0; $i < @strains; $i++) {
					$sum[$i] += $block->{$pos}->{$motif}->{$motif_pos}->{$strains[$i]};
				}
			}
			$vector_motifs = vector(@sum);
			@tag_sum = split('\s+', $tag_counts->{$chr}->{$con_pos});
			$vector_tagcounts = vector(@tag_sum);
			$stddev_m = stddev( $vector_motifs) * 1;
			$stddev_t = stddev( $vector_tagcounts) * 1;
			if($stddev_m < 0.00001) { $stddev_m = 0; }
			if($stddev_t < 0.00001) { $stddev_t = 0; }
			$cor = correlation( $vector_motifs, $vector_tagcounts );
			if($cor ne "n/a" && $stddev_m > 0 && $stddev_t > 0) {
				$correlation{$motif}{$cor}++;
				#Try p-value: 
				if(@strains >= 6) {
					$r = $cor;
					$r2 = $cor * $cor;
					if((1-$r2) < 0 || sqrt((1-$r2)/(@strains - 2)) == 0) {
						$t = 1;
					} else {
						$t = $r/sqrt((1-$r2)/(@strains - 2));
					}
					$p = Statistics::Distributions::tprob(@strains - 1,$t);
					$pvalue{$motif}{$p}++;
				}
			}
		}
	}
	return (\%correlation, \%pvalue);
}
