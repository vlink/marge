!/usr/bin/perl -w
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use strict;
use Getopt::Long;
use config;

$_ = 0 for my($file, $window, $expansion, $help, $index, $last_index, $last_fc, $current_fc, $fc_level, $couple_start, $couple_end, $couple_count, $number_of_peaks, $last_chr, $trend, $last_dir, $fc_expand, $fracOp, $numOp, $fracThreshold, $bed, $one_added, $just_added, $no_add, $true, $min, $max, $couple_index, $print_couple, $couple_count_old, $all);
$_  = () for my(%peaks, @split, @strains, %coupling, %end, %array);
$_ = "" for my ($output, $o);
$fracThreshold = -1;
$fracOp = 0.1;

sub printCMD {
	print STDERR "Usage:\n";
	print STDERR "\t-file: Input file with peaks\n";
	print STDERR "\t-output: prefix (default STDOUT )\n";
	print STDERR "\t-window: window size to start with (default: 2kb)\n";
	print STDERR "\t-expand: maximum expansion (default: unlimited)\n";
	print STDERR "\t-strains: comma seperated list with strains\n";
	print STDERR "\t-fc: fold change level for linking (default: 1.5)\n";
	print STDERR "\t-fc_expand: fold change level for expansion (default: 1.5, option trend overwrites this parameter)\n";
	print STDERR "\t-trend: when expanding peaks just have to go in the right direction, not make the FC threshold anymore (default: off)\n";
	print STDERR "\t-fracOp: fraction of peaks showing fold change in opposing direction (default: 0.1 = 10% - to turn it off set it to -1)\n";
	print STDERR "\t-fracThreshold: Threshold for difference in foldchange if peak does not follow the pattern of the other coupled peaks (relevant when fracOp > 0)\n\t\t(default: 0.3 - Foldchange between peaks can not be outside of 0.7 - 1.3): If no threshold is desired set it to 0\n";
	print STDERR "\t-no_add: Does not allow that 1 peak is added, although fracOp is not reached yet (Default: turned on)\n";
	print STDERR "\t-bed: generates a bed graph (if output is defined output.bed.gz, otherwise bedgraph.bed.gz)\n";
	print STDERR "\t-all: shows statistics about how many peaks were merged (print to STDERR)\n";
	print STDERR "\t-h|--help: print help\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

#Check mandatory command line arguments
my %mandatory = ('-file' => 1, '-strains' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

#Prepare commandline varibale for output in plots
my $commandline = $0 . " ". (join " ", @ARGV);
GetOptions(	"file=s" => \$file, 
		"output=s" => \$output,
		"window=s" => \$window,
		"expand=s" => \$expansion,
		"strains=s{,}" => \@strains,
		"trend" => \$trend,
		"fc=s" => \$fc_level,
		"fc_expand=s" => \$fc_expand,
		"fracOp=s" => \$fracOp,
		"fracThreshold=s" => \$fracThreshold,
		"h" => \$help, 
		"bed" => \$bed,
		"no_add" => \$no_add,
		"all" => \$print_all,
		"help" => \$help)
	or die(&printCMD());

if($help == 1) {
	&printCMD();
}
if($window eq "0") { 
	$window = 2000;
}
if($no_add == 1) {
	$just_added = 1;
}
#Check if window is specified with kb or Mb
if(uc(substr($window, length($window) -2)) eq "KB") {
	$window = substr($window, 0, length($window) - 2) * 1000;
}
if(uc(substr($window, length($window) - 2)) eq "MB") {
	$window = substr($window, 0, length($window) - 2) * 1000000;
}

if(uc(substr($expansion, length($expansion) - 2)) eq "KB") {
	$expansion = substr($expansion, 0, length($expansion) - 2) * 1000;
}
if(uc(substr($expansion, length($expansion) - 2)) eq "MB") {
	$expansion = substr($expansion, 0, length($expansion) - 2) * 1000 * 1000;
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
if($trend == 1) {
	$fc_expand = 1;
}
if($fc_expand == 0) {
	$fc_expand = $fc_level;
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
}

#Read in peaks
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
$index = 0;
#Run through it pairwise
for(my $i = 0; $i < @strains - 1; $i++) {
	for(my $j = $i + 1; $j < @strains; $j++) {
		$number_of_peaks = 1;
		$last_index = 0;
		$couple_count = 0;
		foreach my $chr (sort { $a cmp $b } keys %peaks) {
			$couple_start = 0;
			$couple_index = 0;
			$index = 0;
			$last_index = 0;
			#Sort peaks by position
			foreach my $pos (sort {$a <=> $b} keys %{$peaks{$chr}}) {
				$array{$chr}{$strains[$i]}[$index] = $peaks{$chr}{$pos}{$strains[$i]};
				$array{$chr}{'pos'}[$index] = $pos;
				$array{$chr}{$strains[$j]}[$index] = $peaks{$chr}{$pos}{$strains[$j]};
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
						$couple_start = 0;
						$number_of_peaks = 1;
						$last_dir = 0;
						$print_couple = 0;
						$couple_count++;
					}
				} else {
					#Define last direction that was used for shifting
					if($last_fc > $fc_level && $current_fc > $fc_level && $chr eq $last_chr && ($last_dir == 1 || $last_dir == 0)) {
						$last_dir = 1;
					}
					if($last_fc < (1/$fc_level) && $current_fc < (1/$fc_level) && $chr eq $last_chr && ($last_dir == -1 || $last_dir == 0)) {
						$last_dir = -1;
					}
					if(($last_dir == 1 || $last_dir == -1) && $fc_level > $current_fc && $current_fc > (1/$fc_level)) {
						$print_couple = 1;
					}

					#Save peaks - they are not too far away from each other to be saved, but they changed their foldchange
					if($print_couple == 1) {
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'chr'} = $chr;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'start'} = $couple_start;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'end'} = $couple_end;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'peaks'} = $number_of_peaks;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'index'} = $couple_index;
						$coupling{$strains[$i]}{$strains[$j]}{$couple_count}{'dir'} = $last_dir;
						$couple_start = 0;
						$number_of_peaks = 1;
						$last_dir = 0;
						$print_couple = 0;
						$couple_count++;

					}
					if($print_couple == 0 && $last_dir != 0) {
						if($couple_start == 0) {
							$couple_start = $last_index;
							$couple_index = $index - 1;
						}
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
				$just_added = 0;
				$true = 0;
				if($no_add == 0) {
					$one_added = 0;
				}
				while($true == 0) {
					$start_index--;
					if($start_index < 0) { $true = 1; next; }
					#Peak is too far away, stop the loop
					if($expansion > 0 && $array{$chr}{'pos'}[$start_index + 1] - $array{$chr}{'pos'}[$start_index] > $expansion) {
						$true = 1;
					} else {
						#Save current fold change 
						$new_fc = ($array{$chr}{$strains[$i]}[$start_index] + 1)/($array{$chr}{$strains[$j]}[$start_index] + 1);
						#Compare if new fold change follows the direction of the previously coupled peaks
						if(($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc > $fc_expand) || ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc < (1/$fc_expand))) {
							$peak_number++;
							$just_added = 0;
						#Fold change goes in the wrong direction or is not high enough to be coupled
						} else {
							#No peak is allowed that does not follow the pattern
							if($fracOp == -1) {
								$true = 1;
							#Peak does not follow pattern - check if peaks going in the wrong direction are allowed
							} elsif($numOp/$peak_number < $fracOp) {
								#Now check that it does not differ more than defined by fracThreshold
								#Count up just added - we add peaks that go in the wrong direction as long as ratio of peaks to new added peas < fracOp, so when we do not add a correct peak at the end we have to remove all newly added peaks - keep track of how many peaks got added
								if($fracThreshold == 0) {
									#No threshold is set
									$just_added++;
								} else {
									#Add peak if it is within the interval of deviation
									if($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc - $fracThreshold < (1/$fc_expand)) {
										$just_added++;
										$peak_number++;
									} elsif($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc + $fracThreshold > $fc_expand) {
										$peak_number++;
										$just_added++;
									} else {
										$true = 1;
									}
								}
							#One peak is allowed that goes in the wrong direction - add it here
							} elsif($one_added == 0) {
								$peak_number++;
								$one_added = 1;
								$just_added = 1;
								$numOp--;
							} else {
								$numOp--;
								$true = 1;
								#We can not add another peak that does not follow the pattern - also remove all peaks we just added before this one when the do not follow the pattern
								while($just_added > 0) {
									$peak_number--;
									$start_index++;
									$just_added--;
								}
							} 
							$numOp++;
						}
					}
				}
				while($just_added > 0) {
					$peak_number--;
					$start_index++;
					$just_added--;
				}
				#Last added peak was in the wrong direction - need to remove last peak
				#expand to the right
				$true = 0;
				while($true == 0) {
					$end_index++;
					#The last peak was the last peak on this chromosome - no new peak can be added
					if(!defined $array{$chr}{'pos'}[$end_index] || (defined $array{$chr}{'pos'}[$end_index] && $expansion > 0 && $array{$chr}{'pos'}[$end_index] - $array{$chr}{'pos'}[$end_index - 1] > $expansion)) {
						$true = 1;
						if($just_added == 1) {
							$peak_number--;
							$end_index--;
						}
					} else {
						#Check if peak follows the pattern
						$new_fc = ($array{$chr}{$strains[$i]}[$end_index] + 1)/($array{$chr}{$strains[$j]}[$end_index] + 1);
						if(($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc > $fc_expand) || ($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc < (1/$fc_expand))) {
							$peak_number++;
							$just_added = 0;
						} else {
							#No peak is allowed to not follow the pattern
							if($fracOp == -1) {
								$true = 1;
							} elsif($numOp/$peak_number < $fracOp) {
							#Peak does not follow pattern - check if peaks going in the wrong direction are allowed
								#Now check that it does not differ more than defined by fracThreshold
								if($fracThreshold == 0) {
									#No threshold is set
									$just_added++;
								} else {
									if($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == -1 && $new_fc - $fracThreshold < (1/$fc_expand)) {
										$peak_number++;
										$just_added++;
									} elsif($coupling{$strains[$i]}{$strains[$j]}{$c}{'dir'} == 1 && $new_fc + $fracThreshold > $fc_expand) {
										$peak_number++;
										$just_added++;
									} else {
										$true = 1;
									}
								}
							#One peak can be added although it does not follow the pattern - add this peak here
							} elsif($one_added == 0) {
								$numOp--;
								$one_added = 1;
								$just_added = 1;
							} else {
								$numOp--;
								$true = 1;
								#No peak that does not follow the pattern can be added anymore - remove all peaks that were added before this one and do not follow the pattern
								while($just_added > 0) {
									$peak_number--;
									$end_index--;
									$just_added--;
								}
							} 
						}
					}
				}
				while($just_added > 0) {
					$peak_number--;
					$just_added--;
					$start_index++;
				}
				$numOp = 0;
				#We removed one count from start and added one to end, after we found a peak to expand to - fix that so we don't have to constantly add or remove 1's from here on
				$start_index++;
				$end_index--;
				if($coupling{$strains[$i]}{$strains[$j]}{$c}{'index'} > $start_index) {
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'index'} = $start_index;
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'start'} = $array{$chr}{'pos'}[$start_index];
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'} = $peak_number;	
				}
				if($coupling{$strains[$i]}{$strains[$j]}{$c}{'index'} + $coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'} < $end_index - 1) {
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'end'} = $array{$chr}{'pos'}[$end_index];
					$coupling{$strains[$i]}{$strains[$j]}{$c}{'peaks'} = $peak_number;
				}
			}
			$couple_count_old = $couple_count;
		}
	}
}

$min = 1000;

#Print output file
foreach my $strain1 (keys %coupling) {
	foreach my $strain2 (keys %{$coupling{$strain1}}) {
		if ($output ne "") {
			$o = $output . "_" . $strain1 . "_vs_" . $strain2 . ".txt";
			open OUT, ">$o";
			print OUT "ID(cmd = " . $commandline . ")\tchr\tstart\tend\tstrand\tnumber of coupled peaks\n"; 
			foreach my $count (keys %{$coupling{$strain1}{$strain2}}) {
				print OUT "". $strain1 . "_" . $strain2 . "_" . $count . "\tchr" . $coupling{$strain1}{$strain2}{$count}{'chr'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'start'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'end'} . "\t+\t" . $coupling{$strain1}{$strain2}{$count}{'peaks'} . "\n";
				if($coupling{$strain1}{$strain2}{$count}{'peaks'} < $min) { $min = $coupling{$strain1}{$strain2}{$count}{'peaks'}; }
				if($coupling{$strain1}{$strain2}{$count}{'peaks'} > $max) { $max = $coupling{$strain1}{$strain2}{$count}{'peaks'}; }
				$all += $coupling{$strain1}{$strain2}{$count}{'peaks'};
			}
			close OUT;
		} else {
			print "ID(cmd = " . $commandline . ")\tchr\tstart\tend\tstrand\tnumber of coupled peaks\n"; 
			foreach my $count (keys %{$coupling{$strain1}{$strain2}}) {
				print "". $strain1 . "_" . $strain2 . "_" . $count . "\tchr" . $coupling{$strain1}{$strain2}{$count}{'chr'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'start'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'end'} . "\t+\t" . $coupling{$strain1}{$strain2}{$count}{'peaks'} . "\n";
				if($coupling{$strain1}{$strain2}{$count}{'peaks'} < $min) { $min = $coupling{$strain1}{$strain2}{$count}{'peaks'}; }
				if($coupling{$strain1}{$strain2}{$count}{'peaks'} > $max) { $max = $coupling{$strain1}{$strain2}{$count}{'peaks'}; }
				$all += $coupling{$strain1}{$strain2}{$count}{'peaks'};
			}
		}
	}
}

#print statistics about merging
if($print_all == 1) {
	print STDERR $file . "\t#merged peaks: " . $all . "#least peaks merged: " . $min . "\t#most peaks merged: " . $max . "\n";
}

#Create bed graph
if($bed == 1) {
	foreach my $strain1 (keys %coupling) {
		foreach my $strain2 (keys %{$coupling{$strain1}}) {
			if($output ne "") {
				$o = $output . "_" . $strain1 . "_vs_" . $strain2 . ".bed";
			} else {
				$o = "bedgraph_" . $strain1 . "_vs_" . $strain2 . ".bed";
			}
			open OUT, ">$o";
			print OUT "track name=" . $o . " type=bedGraph\n";
			foreach my $count (keys %{$coupling{$strain1}{$strain2}}) {
				print OUT "chr" . $coupling{$strain1}{$strain2}{$count}{'chr'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'start'} . "\t" . $coupling{$strain1}{$strain2}{$count}{'end'} . "\t+\t" . $coupling{$strain1}{$strain2}{$count}{'peaks'} . "\n";
			}
			close OUT;
			`gzip $o`;
			print STDERR "bedgraph for $strain1 versus $strain2 created!\n";
		}
	}
}
