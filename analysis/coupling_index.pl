#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;

$_ = 0 for my($file, $window, $expansion, $help, $index, $last_index, $last_fc, $current_fc, $fc_level, $couple_start, $couple_end, $couple_count, $number_of_peaks, $last_chr, $trend, $last_dir, $fc_expand, $fracOp, $numOp, $addedFrac, $fracThreshold);
$_  = () for my(%peaks, @split, @strains, %coupling, %end);
$fracThreshold = -1;

sub printCMD {
	print STDERR "Usage:\n";
	print STDERR "\t-file: Input file with peaks\n";
	print STDERR "\t-window: window size to start with (default: 2kb)\n";
	print STDERR "\t-expand: maximum expansion (default: unlimited)\n";
	print STDERR "\t-strains: comma seperated list with strains\n";
	print STDERR "\t-fc: fold change level for linking (default: 1.5)\n";
	print STDERR "\t-fc_expand: fold change level for expansion (default: 1.5, option trend overwrites this parameter)\n";
	print STDERR "\t-trend: when expanding peaks just have to go in the right direction, not make the FC threshold anymore (default: off)\n";
	print STDERR "\t-fracOp: fraction of peaks showing fold change in opposing direction (default: 0.1 = 10%)\n";
	print STDERR "\t-fracThreshold: Threshold for difference in foldchange (default: 0.3 - Foldchange between peaks can not be outside of 0.7 - 1.3): If no threshold is desired set it to 0\n";
	print STDERR "\t-h|--help: print help\n";
}

if(@ARGV < 1) {
	&printCMD();
}

GetOptions(	"file=s" => \$file, 
		"window=s" => \$window,
		"expand=s" => \$expansion,
		"strains=s{,}" => \@strains,
		"trend" => \$trend,
		"fc=s" => \$fc_level,
		"fc_expand=s" => \$fc_expand,
		"fracOp=s" => \$fracOp,
		"fracThreshold=s" => \$fracThreshold,
		"h" => \$help, 
		"help" => \$help)
	or die("Error in command line options!\n");

if($help == 1) {
	&printCMD();
}

if($window == 0) { 
	$window = 2000;
}
#Check if window is specified with kb or Mb
if(substr($window, length($window) -2) eq "kb") {
	$window = substr($window, 0, length($window) - 2) * 1000;
}
if(substr($window, length($window) - 2) eq "Mb") {
	$window = substr($window, 0, length($window) - 2) * 1000000;
}

if(substr($fracOp, length($fracOp) - 1) eq "%") {
	$fracOp = substr($fracOp, 0, length($fracOp) - 1) / 100.0;
}
if($fracOp > 0.99999) {
	$fracOp = $fracOp/100.0;
}
#Set default fracThreshold
if($fracThreshold == -1) {
	$fracThreshold = 0.3;
}
if($fc_level == 0) {
	$fc_level = 1.5;
}

if($fc_expand == 0) {
	$fc_expand = $fc_level;
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
}

open FH, "<$file" or die "Could not open $file: $!\n";
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	$index = 0;
	if(substr($line, 0, 1) eq "#" || substr($line, 0, 6) eq "PeakID") {
		next;
	}
	for(my $i = @split - @strains; $i < @split; $i++) {
		$peaks{substr($split[1], 3)}{$split[2]}{$strains[$index]} = $split[$i];	
		$index++;
	}
	$end{substr($split[1], 3)}{$split[2]} = $split[3];
}

#Save the hash as an array so we can extend in both directions
my %array = ();
$index = 0;
my $couple_index;
my $print_couple = 0;
my $couple_count_old = 0;
#Run through it pairwise
for(my $i = 0; $i < @strains - 1; $i++) {
	for(my $j = $i + 1; $j < @strains; $j++) {
		print $strains[$i] . " vs " . $strains[$j] . "\n";
		$number_of_peaks = 1;
		$last_index = 0;
		$couple_count = 0;
		foreach my $chr (sort { $a cmp $b } keys %peaks) {
			%array = ();
			$couple_start = 0;
			$index = 0;
			my $first_fc;
			#Sort peaks by position
			foreach my $pos (sort {$a <=> $b} keys %{$peaks{$chr}}) {
				$array{$strains[$i]}[$index] = $peaks{$chr}{$pos}{$strains[$i]};
				$array{'pos'}[$index] = $pos;
				$array{$strains[$j]}[$index] = $peaks{$chr}{$pos}{$strains[$j]};
				$current_fc = ($peaks{$chr}{$pos}{$strains[$i]} + 1)/($peaks{$chr}{$pos}{$strains[$j]} + 1);
				#Peaks are too far away from each other
				if($last_index < $pos - $window) {
					#Some previous peaks were merged - save
					if($couple_start != 0) {
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'chr'} = $chr;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'start'} = $couple_start;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'end'} = $couple_end;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'peaks'} = $number_of_peaks;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'index'} = $couple_index;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'dir'} = $last_dir;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'fc'} = $last_fc;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'first_fc'} = $first_fc;
						$couple_start = 0;
						$number_of_peaks = 1;
						$last_dir = 0;
						$print_couple = 0;
						$couple_count++;
					}
				} else {
					#Save peaks - they are not too far away from each other to be saved, but they changed their foldchange
					if($print_couple == 1) {
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'chr'} = $chr;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'start'} = $couple_start;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'end'} = $couple_end;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'peaks'} = $number_of_peaks;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'index'} = $couple_index;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'dir'} = $last_dir;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'fc'} = $last_fc;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'first_fc'} = $first_fc;
						$couple_start = 0;
						$number_of_peaks = 1;
						$last_dir = 0;
						$print_couple = 0;
						$couple_count++;

					}
					if($last_fc > $fc_level && $current_fc > $fc_level && $chr eq $last_chr && ($last_dir == 1 || $last_dir == 0)) {
						$last_dir = 1;
					}
					if($last_fc < (1/$fc_level) && $current_fc < (1/$fc_level) && $chr eq $last_chr && ($last_dir == -1 || $last_dir == 0)) {
						$last_dir = -1;
					}
					if(($last_dir == 1 || $last_dir == -1) && $fc_level > $current_fc && $current_fc > (1/$fc_level)) {
						$print_couple = 1;
					}
					if($print_couple == 0 && $last_dir != 0) {
						if($couple_start == 0) {
							$couple_start = $last_index;
							$couple_index = $index - 1;
							$first_fc = $last_fc;
						}
						#Add 150bp because of peak size
						$couple_end = $end{$chr}{$pos};
						$number_of_peaks++;
						
					}
				}
				$last_fc = $current_fc;
				$last_index = $pos;
				$last_chr = $chr;
				$index++;
			}
			#Time to expand coupling index
			for(my $c = $couple_count_old; $c < $couple_count; $c++) {
				my $start_index = $coupling{$strains[$i]}{$strains[$j]}{$c}{'index'};
				my $end_index = $start_index + $coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'};
				my $new_fc;
				#First expand to the left
				my $expansion_start = $coupling{$strains[$i]}{$strains[$j]}{$c}{'start'};
				my $peak_number = $coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'};
				#Foldchange needs to overlap - just right direction?
				#expand to the left
				my $true = 0;
				while($true == 0) {
					$start_index--;
					#Peak is too far away, stop the loop
					if($expansion > 0 && $array{'pos'}[$start_index + 1] - $array{'pos'}[$start_index] > $expansion) {
						$true = 1;
					} else {
						$new_fc = ($array{$strains[$i]}[$start_index] + 1)/($array{$strains[$j]}[$start_index] + 1);
						if($trend == 0 && ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc > $fc_expand) || ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc < (1/$fc_expand))) {
							$peak_number++;
							$addedFrac = 0;
						} elsif($trend == 1 && (($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc > 1) || ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc < 1))) {
							$peak_number++;
							$addedFrac = 0;
						} else {
							#Peak does not follow pattern - check if peaks going in the wrong direction are allowed
							$numOp++;
							if($numOp/$peak_number < $fracOp) {
								#Now check that it does not differ more than defined by fracThreshold
								if($fracThreshold == 0) {
									#No threshold is set
									$addedFrac++;
								} else {
									if($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc - $fracThreshold < (1/$fc_expand)) {
										$addedFrac++;
									} elsif($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc + $fracThreshold > $fc_expand) {
										$addedFrac++;
									} else {
										$true = 1;
									}
								}
							} else {
								$true = 1;
							} 
						}
					}
				}
				#Last added peak was in the wrong direction - need to remove last peak
				while($addedFrac > 0) {
					$peak_number--;
					$start_index++;
					$addedFrac--;
				}

				#expand to the right
				$true = 0;
				while($true == 0) {
					$end_index++;
					if(!defined $array{'pos'}[$end_index] || (defined $array{'pos'}[$end_index] && $expansion > 0 && $array{'pos'}[$end_index] - $array{'pos'}[$end_index - 1] > $expansion)) {
						$true = 1;
						if($addedFrac == 1) {
							$peak_number--;
							$end_index--;
						}
					} else {
						$new_fc = ($array{$strains[$i]}[$end_index] + 1)/($array{$strains[$j]}[$end_index] + 1);
						if($trend == 0 && ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc > $fc_expand) || ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc < (1/$fc_expand))) {
							$peak_number++;
						} elsif($trend == 1 && (($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc > 1) || ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc < 1))) {
							$peak_number++;
						} else {
							#Peak does not follow pattern - check if peaks going in the wrong direction are allowed
							$numOp++;
							if($numOp/$peak_number < $fracOp) {
								#Now check that it does not differ more than defined by fracThreshold
								if($fracThreshold == 0) {
									#No threshold is set
									$addedFrac++;
								} else {
									if($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc - $fracThreshold < (1/$fc_expand)) {
										$addedFrac++;
									} elsif($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc + $fracThreshold > $fc_expand) {
										$addedFrac++;
									} else {
										$true = 1;
									}
								}
							} else {
								$true = 1;
							} 
						}
					}
				}
				while($addedFrac > 0) {
					$peak_number--;
					$start_index++;
					$addedFrac--;
				}

				$numOp = 0;
				#We removed one count from start and added one to end, after we found a peak to expand to - fix that so we don't have to constantly add or remove 1's from here on
				$start_index++;
				$end_index--;
				if($coupling{$strains[$i]}{$strains[$j]}{$c}{'index'} > $start_index) {
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'index'} = $start_index;
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'start'} = $array{'pos'}[$start_index];
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'} = $peak_number;	
				}
				if($coupling{$strains[$i]}{$strains[$j]}{$c}{'index'} + $coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'} < $end_index - 1) {
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'end'} = $array{'pos'}[$end_index];
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'} = $peak_number;
				}
			}
			$couple_count_old = $couple_count;
		}
	}
}



my $output;
foreach my $strain1 (keys %coupling) {
	foreach my $strain2 (keys %{$coupling{$strain1}}) {
		$output = "coupled_peaks_" . $strain1 . "_" . $strain2 . ".txt";
		open OUT, ">$output";
		print STDERR $strain1 . "\t" . $strain2 . "\n";
		foreach my $count (keys %{$coupling{$strain1}{$strain2}}) {
			print OUT "". $strain1 . "_" . $strain2 . "_" . $count . "\tchr" . $coupling{$strain1}{$strain2}{$count}{'chr'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'start'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'end'} . "\t+\t" . $coupling{$strain1}{$strain2}{$count}{'peaks'} . "\n";
		}
		close OUT;
	}
}
