#!/usr/bin/perl

BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
package analysis;
use strict;
use Getopt::Long;
#require '/Users/verenalink/workspace/strains/general/config.pm';
use config;
use Data::Dumper;
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
	my %seq = %{$_[1]};
	my %save_local_shift = %{$_[2]};
	my $fileHandlesMotif_ref = $_[3];
	my %tag_counts = %{$_[4]};
	my @strains = @{$_[5]};
	my %index_motifs = %{$_[6]};
	my $ab = $_[7];
	my %remove = %{$_[8]};
	my $fc_significant = $_[9];
	my $mut_pos = $_[10];
	my $overlap = $_[11];
	my %fc = ();
	my $ab_mut = 0;
	if(@_ > 12) {
		%fc = %{$_[12]};
	}
        open FH, "<$motif_file";
        my $last_line = "";
	my %block = ();
        foreach my $line (<FH>) {
                chomp $line;
                @split = split('\t', $line);
		$split[3] = &get_motif_name($split[3]);
                @header = split('_', $split[0]);
               #Read in blocks - stop after gathering all information about one position
                if($last_line ne "" && $header[0] . "_" . $header[1] . "_" . $header[2] ne $last_line) {
                        my @a = split("_", $last_line);
			print Dumper %block;
                        $ab_mut = &print_results(substr($a[0], 3), $a[1], \%block, $fileHandlesMotif_ref, $last_line, \%tag_counts, \@strains, \%index_motifs, \%fc, $ab, $fc_significant, $mut_pos, $overlap);
		#	print "\n\n";
		#	print Dumper %fc;
		#	&create_cytoscape_file(\%block, \%fc, $last_line, \@strains);
			&save_all_motifs(\%block, $last_line);
			if($ab_mut == 1) {
				$remove{$last_line} = 1;
			}
			%block=();
                }
                $last_line = $header[0] . "_" . $header[1] . "_" . $header[2];
                my $half = length($seq{$split[0]})/2;
                $half = int($half);
                $motif_start = $header[1] + $half;
                #First get start pos of the motif
                if($split[4] eq "-" && $split[1] < 0) {
                #       $motif_start = $motif_start - ((length($split[2])/2)) + $split[1];      
                        $motif_start = $motif_start - (length($split[2])) + $split[1] + 1;
                }
                if($split[4] eq "-" && $split[1] >= 0) {
                #       $motif_start = $motif_start - ((length($split[2])/2)) + $split[1];
                        $motif_start = $motif_start - (length($split[2])) + $split[1] + 1;
                }
                if($split[4] eq "+") {
                        $motif_start = $motif_start + $split[1];
                }
                $motif_start = $motif_start - $header[1];
                $block{$split[3]}{$motif_start + $save_local_shift{$split[0]}{$motif_start}}{$header[-1]} = $split[-1];
                $block{$split[3]}{$motif_start + $save_local_shift{$split[0]}{$motif_start}}{'length'} = length($split[2]);
        }
        my @a = split("_", $last_line);
	$ab_mut = &print_results(substr($a[0], 3), $a[1], \%block, $fileHandlesMotif_ref, $last_line, \%tag_counts, \@strains, \%index_motifs, \%fc, $ab, $fc_significant, $mut_pos, $overlap);
#	&save_all_motifs(\%block, $last_line);
#	print "network saved\n";
#	&create_cytoscape_file(\%fc, \@strains);
#	exit;
#	print "cytoscape important\n";
#	print Dumper %cytoscape_unsig;
#	print "\n\n";
#	open NET, ">network_significant.sif";
#	my $sum_keys_sig = 0;
#	foreach my $key (keys %cytoscape_sig) {
#		foreach my $key2 (keys %{$cytoscape_sig{$key}}) {
#			$sum_keys_sig += $cytoscape_sig{$key}{$key2};
#		}
#	}
#	foreach my $key (keys %cytoscape_sig) {
#		foreach my $key2 (keys %{$cytoscape_sig{$key}}) {
#		#	print NET $key . " " . (int(($cytoscape_sig{$key}{$key2}/$sum_keys_sig) * 100) + 1). " " . $key2 . "\n";
#			print NET $key . " " . $cytoscape_sig{$key}{$key2} . " " . $key2 . "\n";
#		}
#	}
#	close NET;
#	open NET, ">network_unsignificant.sif";
#	my $sum_keys_unsig = 0;
#	foreach my $key (keys %cytoscape_unsig) {
#		foreach my $key2 (keys %{$cytoscape_unsig{$key}}) {
#			$sum_keys_unsig += $cytoscape_unsig{$key}{$key2};
#		}
#	}
#	foreach my $key (keys %cytoscape_unsig) {
#		foreach my $key2 (keys %{$cytoscape_unsig{$key}}) {
#		#	print NET $key . " " . (int(($cytoscape_unsig{$key}{$key2}/$sum_keys_unsig) * 100) + 1). " " . $key2 . "\n";	
#			print NET $key . " " . $cytoscape_unsig{$key}{$key2} . " " . $key2 . "\n";	
#		}
#	}
#	close NET;
#	&network_analysis(\%fc, $fc_significant, \@strains, \%index_motifs);
	if($ab_mut == 1) {
		$remove{$last_line} = 1;
	}
	return (\%remove, \%matrix_pos_muts);
}

sub create_cytoscape_file{
	my $fc_ref = $_[0];
	my $strain_ref = $_[1];
	my @strains = @{$strain_ref};
	my %fc = %$fc_ref;
	my $sig = 0;
	my @array = ();
	my @pos = ();
	my $first = "";
	my $second = "";
	my $add_first = "";
	my $add_second = "";
	foreach my $peak (sort {$a cmp $b} keys %network_save) {
		@pos = ();
		@array = ();
		foreach my $pos (sort {$a <=> $b} keys %{$network_save{$peak}}) {
			foreach my $motif (sort {$a cmp $b} keys %{$network_save{$peak}{$pos}}) {
				push @array, $motif;
				push @pos, $pos;
			}
		}
		@split = split("_", $peak);
		$sig = 0;
		if($fc{substr($split[0], 3)}{$split[1]} < -1 || $fc{substr($split[0], 3)}{$split[1]} > 1) {
			$sig = 1;
		}
		for(my $i = 0; $i < @array - 1; $i++) {
			for(my $j = $i + 1; $j < @array; $j++) {
				if(!exists $network_save{$peak}{$pos[$i]}{$array[$i]}{$strains[0]} || !exists $network_save{$peak}{$pos[$i]}{$array[$i]}{$strains[1]} || $network_save{$peak}{$pos[$i]}{$array[$i]}{$strains[0]} == 0 || $network_save{$peak}{$pos[$i]}{$array[$i]}{$strains[1]} == 1) { 
					$first = "mut_" . $array[$i];		
				} else {
					$first = $array[$i];
				}
				if(!exists $network_save{$peak}{$pos[$j]}{$array[$j]}{$strains[0]} || !exists $network_save{$peak}{$pos[$j]}{$array[$j]}{$strains[1]} || $network_save{$peak}{$pos[$j]}{$array[$j]}{$strains[0]} == 0 || $network_save{$peak}{$pos[$j]}{$array[$j]}{$strains[1]} == 0) {
					$second = "mut_" . $array[$j];
				} else {
					$second = $array[$j];
				}
				if($array[$i] le $array[$j]) {
					$add_first = $first;
					$add_second = $second;
				} else {
					$add_first = $second;
					$add_second = $first;
				}
				if($sig == 1) {
					$cytoscape_sig{$add_first}{$add_second}++;
				} else {
					$cytoscape_unsig{$add_first}{$add_second}++;
				}
			}
		}
	}
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

sub print_results {
        my $curr_chr = $_[0];
        my $curr_pos = $_[1];
	my $block_ref = $_[2];
	my %block = %$block_ref;
	my @fileHandlesMotif = @{$_[3]};
	my $last_line = $_[4];
	my %tag_counts = %{$_[5]};
	my @strains = @{$_[6]};
	my %index_motifs = %{$_[7]};
	my %fc = %{$_[8]};
	my $fc_exists = 0;
	my $ab = $_[9];
	my $fc_significant = $_[10];
	my $mut_pos = $_[11];
	my $overlap = $_[12];
	my $ab_mut = 0;
	my $num_of_bp = 0;
	if((keys %fc) > 1) {
		$fc_exists = 1;
	}
        #Now run through all motifs and see if they are in all the other strains
        my %ignore = ();
        foreach my $motif (keys %block) {
                %existance = ();
                %diff = ();
                foreach my $motif_pos (keys %{$block{$motif}}) {
			if($overlap eq "half") {
				$num_of_bp = int($block{$motif}{$motif_pos}{'length'}/2)
			} elsif($overlap eq "complete") {
				$num_of_bp = 0;
			} else {
				$num_of_bp = ($overlap*1);
			}
		#	print "overlap: " . $overlap . "\n";
		#	print "num of bp: " . $num_of_bp . "\n";
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
					for(my $run_motif = $motif_pos - $num_of_bp; $run_motif < $motif_pos + $num_of_bp; $run_motif++) {
                                                if(exists $block{$motif}{$run_motif} && exists $block{$motif}{$run_motif}{$strains[$i]} && $run_motif != $motif_pos) {
				#			print "overlap found!\n";
                                                        $block{$motif}{$motif_pos}{$strains[$i]} = $block{$motif}{$run_motif}{$strains[$i]};
							$diff{$strains[$i]} += $block{$motif}{$motif_pos}{$strains[$i]};
                                                        delete $existance{$strains[$i]};
                                                        delete $block{$motif}{$run_motif}{$strains[$i]};
                                                        $ignore{$run_motif} = 1;
                                                        if(exists $block{$motif}{$run_motif}{'length'}) {
                                                                delete $block{$motif}{$run_motif}{'length'};
                                                        }
                                                }
                                        }
                                }
				if($block{$motif}{$motif_pos}{$strains[$i]} eq "" || !exists $block{$motif}{$motif_pos}{$strains[$i]}) {
					$block{$motif}{$motif_pos}{$strains[$i]} = &calculate_motif_score($motif, substr($seq{$last_line . "_" . $strains[$i]}, $motif_pos, $block{$motif}{$motif_pos}{'length'}));
					$diff{$strains[$i]} += $block{$motif}{$motif_pos}{$strains[$i]};
				} else {
					$diff{$strains[$i]} += $block{$motif}{$motif_pos}{$strains[$i]};
				}
			}
                }
                $fileHandlesMotif[$index_motifs{$motif}]->print($motif . "\t" . $last_line . "\t");
                $fileHandlesMotif[$index_motifs{$motif}]->print($tag_counts{$curr_chr}{$curr_pos});
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
                $fileHandlesMotif[$index_motifs{$motif}]->print("" . (keys %{$block{$motif}}) . "\n");
		if($mut_pos == 1 && $fc_exists == 1) {
			&analyze_motif_pos(\%block, $last_line, $fc{$curr_chr}{$curr_pos}, \@strains, $fc_significant);
		}
        }
        %block = ();
	return $ab_mut;
}

sub analyze_motif_pos{
	my %block = %{$_[0]};
	my $last_line = $_[1];
	my $fc = $_[2];
	my @strains = @{$_[3]};
	my $fc_significant = $_[4];
	my @fc_split = split('\t', $fc);
	my @char_one;
	my @char_two;
	my $num_of_muts = 0;
	foreach my $motif (keys %block) {
		foreach my $motif_pos (keys %{$block{$motif}}) {
			for(my $i = 0; $i < @strains - 1; $i++) {
				for(my $j = $i + 1; $j < @strains; $j++) {
					if($block{$motif}{$motif_pos}{$strains[$i]} != $block{$motif}{$motif_pos}{$strains[$j]}) {
						@char_one = split("", substr($seq{$last_line . "_" . $strains[$i]}, $motif_pos, $block{$motif}{$motif_pos}{'length'}));
						@char_two = split("", substr($seq{$last_line . "_" . $strains[$j]}, $motif_pos, $block{$motif}{$motif_pos}{'length'}));
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
									if($block{$motif}{$motif_pos}{$strains[$i]} > $block{$motif}{$motif_pos}{$strains[$j]}) {
										$matrix_pos_muts{$motif}{$char_two[$c]}{$c}{'S'}++;
									} else {
										$matrix_pos_muts{$motif}{$char_one[$c]}{$c}{'S'}++;
									}
								} else {
									if($block{$motif}{$motif_pos}{$strains[$i]} > $block{$motif}{$motif_pos}{$strains[$j]}) {
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

sub read_motifs{
	open FH, "<$_[0]";
	my $output_mut = $_[1];
	my %del = %{$_[2]};
	my $pos = 0;
	my $motif = "";
	my $index = 0;
	my $name;
	my %index_motifs;
	my @fileHandlesMotif;
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			$motif = (split('\t', $line))[1];
			$pos = 0;
			#Open file per motif
			$name = &get_motif_name($motif);
			my $filename = $output_mut . "_" . $name . ".txt";
			$del{$filename} = 1;
			open my $fh, ">", $filename or die "Can't open $filename: $!\n";
			$index_motifs{$name} = $index;
			$fileHandlesMotif[$index] = $fh;
			$index++;
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
	return (\%index_motifs, \%PWM, \@fileHandlesMotif, \%del);
}

sub scan_motif_with_homer{
	my $seq_file = $_[0];
	my $out_file = $_[1];
	my $tf = $_[2];
	my $command = "homer2 find -i " . $seq_file . " -m " . $tf . " > " . $out_file;
	print STDERR $command . "\n";
	`$command`;	
}

sub get_seq_for_peaks {
	open OUT, ">$_[0]";
	my %peaks = %{$_[1]};
	my @strains = @{$_[2]};
	my $data = $_[3];
	my $allele = $_[4];
	my $line_number = $_[5];
	my $mut_only = $_[6];
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
	my $h;
	foreach my $p (keys %peaks) {
		#Start with offset for first line
		$general_offset = 5 + length($p);
		if(length($p) > 3) { next; }
		for(my $i = 0; $i < @strains; $i++) {
			#Read in the offsets for the different strains
			$strains[$i] =~ s/,//g;
			$strains[$i] = uc($strains[$i]);
			$current_pos{$strains[$i]} = 0;
			&create_offset($p, $allele, uc($strains[$i]), $data);
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
				for(my $run = $current_pos{$strains[$i]}; $run < @{$shift{$strains[$i]}}; $run++) {
					#Position of our peak is greater than the last mutation, so we just use the last shifting value
					if($working_start > $mutation_pos{$strains[$i]}[-1]) {
						$strain_spec_start = $working_start + $shift{$strains[$i]}[-1];
						$strain_spec_stop = $working_stop + $shift{$strains[$i]}[-1];
						$shift_start = 1;
						$shift_stop = 1;
					}
					#No shifting of position yet, and mutation pos greater than our current pos - so we use one position earlier and use this shifting vector
					if($shift_start == 0 && $mutation_pos{$strains[$i]}[$run] > $working_start) {
						if($run == 0) {
							$current_pos{$strains[$i]} = 0;
							$strain_spec_start = $working_start;
							$mut_number_start = $run;
						} else {
							$current_pos{$strains[$i]} = $run;
							$strain_spec_start = $working_start + $shift{$strains[$i]}[$run-1];
							$mut_number_start = $run - 1;
						}
						$shift_start = 1;
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
			#	print OUT ">chr" . $p . "_" . $start . "_" . $peaks{$p}{$start} . "_" . $strains[$i] . "\n";
				$length = $strain_spec_stop - $strain_spec_start;
				$newlines = int($strain_spec_stop/50) - int($strain_spec_start/50);
				$length = $length + $newlines;
				$newlines = int(($strain_spec_start)/50);
				$byte_offset = $strain_spec_start + $general_offset + $newlines;
				seek($fileHandles[$i], $byte_offset, 0);
				read $fileHandles[$i], $seq, $length;
				$seq =~ s/\n//g;
				$current_pos{$strains[$i]} = uc($seq);
			#	print OUT uc($seq) . "\n";
				$seq{$header} = uc($seq);
				$seq = "";
				my $local_shift = 0;
				my $insert = 0;
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
							$local_shift = ($shift{$strains[$i]}[$tmp_mut - 1] - $shift{$strains[$i]}[$tmp_mut]);
							$last_shift = $last_shift + $local_shift;
							$tmp_mut++;
							$before_vector = $mutation_pos{$strains[$i]}[$tmp_mut];
							$insert = 0;
						}
						if($last_shift >= 0 && $insert == 0) {
							for(my $k = 0; $k <= $last_shift + 1; $k++) {
								$save_local_shift{$header}{$j - $k} = $last_shift;
							}
							$insert = 1;
						} else {
							if($local_shift < 0) {
								for(my $k = 0; $k > $last_shift; $k--) {
									$save_local_shift{$header}{$j + $index_shift - $k} = $last_shift;
								}
								$save_local_shift{$header}{$j + $index_shift - $last_shift} = $last_shift;
								$index_shift = 0;
								$local_shift = 0;
							} else {
								$save_local_shift{$header}{$j + $index_shift} = $last_shift;
							}
						}
					}
				}
			}
			#$iPrint here
			if($mut_only == 1) {
				$equal = 0;
				$ref_seq = $current_pos{$strains[0]};
				for(my $i = 1; $i < @strains; $i++) {
					if($current_pos{$strains[$i]} ne $ref_seq) {
						$equal = 1;
					}
				}
				if($equal == 0) {
					for(my $i = 0; $i < @strains; $i++) {
						$h = "chr" . $p . "_" . $start . "_" . $peaks{$p}{$start} . "_" . $strains[$i];	
						delete $save_local_shift{$h};
						delete $seq{$h};
					}	
				}
			}
			$seq_number++;
		#	print STDERR "" . int(($seq_number/$line_number) * 100) . "% of sequences are gathered\r";
			print STDERR "" . $seq_number . " of " . $line_number . " gathered\r";
		}
	}
	foreach my $key (sort {$a cmp $b} keys %seq) {
		print OUT ">" . $key . "\n";
		print OUT $seq{$key}  . "\n";
	}
	close OUT;
	return(\%seq, \%save_local_shift);
}

sub create_offset{
	my $chr = $_[0];
	my $a = $_[1];
	my $strain = $_[2];
	my $data = $_[3];
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
