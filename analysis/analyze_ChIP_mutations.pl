#!/usr/bin/perl


use strict;
use Getopt::Long;
#require '/Users/verenalink/workspace/strains/general/config.pm';
require '../general/config.pm';
use Data::Dumper;
#require '../db_part/database_interaction.pm';


$_ = "" for my($genome, $file, $tf, $filename, $last_line);
$_ = () for my(@strains, %peaks, @split, %mutation_pos, %shift, %current_pos, %save_local_shift, %seq, %PWM, @fileHandlesMotif, %index_motifs, %tag_counts, %fc, @header, %block, %analysis_result, %existance, %diff, %ranked_order, %mut_one, %mut_two, %delta_score);
$_ = 0 for my($homo, $allele, $region, $motif_score, $motif_start, $more_motifs, $save_pos, $delta);
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
        print STDERR "\t-strains <strains>: Comma-separated list of strains - Order must overlay with order in annotated peak file\n";
        print STDERR "\t-file <file>: annotated peak file (including tag counts)\n";
	print STDERR "\t-TF <transcription factor motif matrix>: Matrix of the TF that was chipped for\n";
        print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-region: Size of the region used to look for other motifs (Default: 200)\n";
	print STDERR "\t-delta: Uses motif score differences instead of binary existance\n";
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
		"-region=s" => \$region,
		"-delta" => \$delta)
        or die("Error in command line options!\n");
#First step: Get the sequences for the peaks
$allele = 1;
my $ref_save = 0;

#Save motif files
&read_motifs($tf);

#Make sure the reference is in the strains array
#my $reference = config::reference_genome($genome);
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
#	if($strains[$i] eq $reference) {
#		$ref_save = 1;
#	}
}
#if($ref_save == 0) {
#	push @strains, $reference;
#}


print STDERR "Saving peaks\n";
open FH, "<$file";
my $line_number = 0;
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#" || substr($line, 0, 6) eq "PeakID") {
		next;
	}
	@split = split('\t', $line);
	if(@split < 19) {
		print STDERR "This is not an annotated peak file!\n";
		exit;
	}
	if(length($split[1]) > 10) { next; }
	$peaks{substr($split[1], 3)}{$split[2]} = $split[3];
	$line_number++;
	$tag_counts{substr($split[1], 3)}{$split[2]} = "";		
	$fc{substr($split[1], 3)}{$split[2]} = "";
	for(my $i = 19; $i < @split; $i++) {
		$tag_counts{substr($split[1], 3)}{$split[2]} .= $split[$i] . "\t";
	}
	for(my $i = 19; $i < @split - 1; $i++) {
		for(my $j = $i + 1; $j < @split; $j++) {
			$fc{substr($split[1], 3)}{$split[2]} .= log(($split[$i]+1)/($split[$j]+1))/log(2) . "\t";
		} 
	}
#	print "tag counts\n";
#	print $tag_counts{substr($split[1], 3)}{$split[2]} . "\n";
#	print "fold change:\n";
#	print $fc{substr($split[1], 3)}{$split[2]} . "\n";
}
close FH;
#exit;

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
	for(my $i = 0; $i < @strains - 1; $i++) {
		for(my $j = $i + 1; $j < @strains; $j++) {
			$fileHandlesMotif[$files]->print("fc " . $strains[$i] . " vs " . $strains[$j] . "\t");
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
my ($c_chr, $c_pos) = &analyze_motifs();
&print_results($c_chr, $c_pos);

for(my $i = 0; $i < @fileHandlesMotif; $i++) {
	close $fileHandlesMotif[$i];
}

print STDERR "Generating R files!\n";
&generate_R_files();

#Generate R files
sub generate_R_files {
	my @short;
	my $filename;
	my $f = 0;
	my $num_of_muts;

	print STDERR "Add filtering for number of motifs!\n";
	open R, ">output_all_plots.R";
	print R "pdf(\"all_plots.pdf\", width=10, height=5)\n";
	my $cal_pvalue = 0;
	my $count_obs_one = 0;
	my $count_obs_two = 0;
	my $count_obs_both = 0;
	foreach my $motif (keys %index_motifs) {
		@short = split /[\(,\/]/, $motif;
		#Make sure there a mutations in this motif
		$filename = "output_mutation_" . $short[0] . ".txt";
		$num_of_muts = `wc -l $filename`;
		@split = split('\s+', $num_of_muts);
		if($split[0] == 1) {
			print STDERR "No mutations in motif for " . $short[0] . "\n";
			next;
		}	
		open FH, "<$filename" or die "Can't find $filename: $!\n";
		print R "par(oma=c(0,0,0,0))\n";
		$f = 0;
		%ranked_order = ();
		%mut_one = ();
		%mut_two = ();
		my $h;
		my $max = 0;
		my $min = 0;
		my $max_delta = 0;
		my $min_delta = 0;
		foreach my $line (<FH>) {
			if($f == 0) {
				$h = $line;
				$f++;
				next;
			}
			chomp $line;
			@split = split('\t', $line);
			for(my $i = 0; $i < @strains - 1; $i++) {
				for(my $j = $i + 1; $j < @strains; $j++) {
					$ranked_order{$strains[$i] . "_" . $strains[$j]}{$split[1]} = $split[1 + @strains + $i + $j]; 
					if($split[1 + @strains + $i + $j] > $max) {
						$max = $split[1 + @strains + $i + $j];
					}
					if($split[1 + @strains + $i + $j] < $min) {
						$min = $split[1 + @strains + $i + $j];
					}
					if($split[8 + @strains + $i] eq "1") {
						$mut_one{$strains[$i] . "_" . $strains[$j]}{$split[1]} = 1;
					}
					if($split[8 + @strains + $j] eq "1") {
						$mut_two{$strains[$i] . "_" . $strains[$j]}{$split[1]} = 1;
					}
					$delta_score{$strains[$i] . "_" . $strains[$j]}{$split[1]} = $split[8 + $i] - $split[8 + $j];
					if($delta_score{$strains[$i] . "_" . $strains[$j]}{$split[1]} > $max_delta) {
						$max_delta = $delta_score{$strains[$i] . "_" . $strains[$j]}{$split[1]};
					}
					if($delta_score{$strains[$i] . "_" . $strains[$j]}{$split[1]} < $min_delta) {
						$min_delta = $delta_score{$strains[$i] . "_" . $strains[$j]}{$split[1]};
					}
				}
			}
		}
		for(my $i = 0; $i < @strains - 1; $i++) {
			for(my $j = $i + 1; $j < @strains; $j++) {
				my $dist_one = "one <- c(";
				my $dist_two = "two <- c(";
				my $dist_no_mut = "no_mut <- c(";
				my $delta_one = "delta_one <- c(";
				my $delta_two = "delta_two <- c(";
				my $delta_no_mut = "delta_no_mut <- c(";
				$filename = "plot_" . $short[0] . "_" . $strains[$i] . "_vs_" . $strains[$j] . ".R";
			#	open R, ">$filename";
			#	print "output file: " . $filename . "\n";
				my $x_count = 1;
				my $tmp_delta;
				if($delta == 1) {
					$max = $max + &log_fc($max_delta);
					$min = $min + &log_fc($min_delta);
				}
				$cal_pvalue = 0;
				my $add_additional =  ((keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) * 0.25);
				if($add_additional < 20) { $add_additional = 20; }
				my $width_boxplot = $add_additional;
				if($add_additional < 31) { $width_boxplot = 45; }
				my $number_peaks = (keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}});
				print R "plot(c(0," . ((keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) + $add_additional) . "), c(" . $min . ", " . $max . "+1), col=\"white\", xlab=\"rank-ordered peaks\", ylab=\"FC " . $strains[$i] . " vs " . $strains[$j] . "\", main=\"" . $strains[$i] . " vs " . $strains[$j] . "\\n$short[0]\")\n"; 
				#Sort by fold change value between strains
				foreach my $rank (sort {$ranked_order{$strains[$i] . "_" . $strains[$j]}{$a} <=> $ranked_order{$strains[$i] . "_" . $strains[$j]}{$b}} keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) {
					if(exists $mut_one{$strains[$i] . "_" . $strains[$j]}{$rank}) {
						$dist_one .= sprintf("%.2f", $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank}) . ",";
						if($delta == 1) {
							print R "points(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", pch=8, col=\"red\")\n";
						 	print R "segments(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", x1 = " . $x_count . ", y1= " . ($ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} + &log_fc($delta_score{$strains[$i] . "_" . $strains[$j]}{$rank})) . ", col=\"red\")\n";
							$delta_one .= sprintf("%.2f", &log_fc($delta_score{$strains[$i] . "_" . $strains[$j]}{$rank})) . ",";
						} else {
							print R "segments(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . " - 0.5" . ", x1 = " . $x_count . ", y1= " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . "+0.5, col=\"red\")\n";
						}
					}
					if(exists $mut_two{$strains[$i] . "_" . $strains[$j]}{$rank}) {
						$dist_two .= sprintf("%.2f", $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank}) . ",";
						if($delta == 1) {
							print R "points(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", pch=8, col=\"blue\")\n";
							print R "segments(" . $x_count . " + 0.1, " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", x1 = " . $x_count . " + 0.1" . ", y1= " . ($ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} + &log_fc($delta_score{$strains[$i] . "_" . $strains[$j]}{$rank})) . ", col=\"blue\")\n";
							$delta_two .= sprintf("%.2f", &log_fc($ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank})) . ",";

						} else {
							print R "segments(" . $x_count . " + 0.1, " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", x1 = " . $x_count . " + 0.1" . ", y1= " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . "+1, col=\"blue\")\n";
						}
					}
					if(!exists $mut_one{$strains[$i] . "_" . $strains[$j]}{$rank} && !exists $mut_two{$strains[$i] . "_" . $strains[$j]}{$rank}) {
						$dist_no_mut .= sprintf("%.2f", $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank}) . ",";
						if($delta == 1) {
							$delta_no_mut .= sprintf("%.2f", &log_fc($ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank})) . ",";
						}
					}
					$x_count++;	
				}
				print R "legend(\"topleft\", c(\"" . $strains[$i] . "\",\"" . $strains[$j] . "\"), col=c(\"red\", \"blue\"), pch=16)\n";
				if(length($dist_one) > 10) {
					chop $dist_one;
				}
				if(length($dist_two) > 10) {
					chop $dist_two;
				}
				if(length($dist_no_mut) > 10) {
					chop $dist_no_mut;
				}
				$dist_one .= ")";
				$dist_two .= ")";
				$dist_no_mut .= ")";
				$count_obs_one = () = $dist_one =~ /\,/gi;
				$count_obs_two = () = $dist_two =~ /\,/gi;
				$count_obs_both = () = $dist_no_mut =~ /\,/gi;
				if($count_obs_one >= 4 && $count_obs_two >= 4 && $count_obs_both >= 4) {
					$cal_pvalue = 1;
				}
				print R $dist_one . "\n";
				print R $dist_two . "\n";
				print R $dist_no_mut . "\n";
				if($delta == 1) {
					$cal_pvalue = 0;
					if(length($delta_one) > 10) {
						chop $delta_one;
					}
					if(length($delta_two) > 10) {
						chop $delta_two;
					}
					if(length($delta_no_mut) > 10) {
						chop $delta_no_mut;
					}
					$delta_one .= ")";
					$delta_two .= ")";
					$delta_no_mut .= ")";
					$count_obs_one = () = $delta_one =~ /\,/gi;
					$count_obs_two = () = $delta_two =~ /\,/gi;
					$count_obs_both = () = $delta_no_mut =~ /\,/gi;
					if($count_obs_one >= 4 && $count_obs_two >= 4 && $count_obs_both >= 4) {
						$cal_pvalue = 1;
					}
					print R $delta_one . "\n";
					print R $delta_two . "\n";
					print R $delta_no_mut . "\n";
					if($cal_pvalue == 1) {
					#	print R "boxplot(c(delta_no_mut), c(delta_one), c(delta_two), boxwex=" . (int($add_additional/3) - 10) . ", add=TRUE, at=c(" . ($number_peaks + int($add_additional/3)) . "," . ($number_peaks + (int($add_additional/3) * 2)) . "," . ($number_peaks + (int($add_additional/3) * 3)) . "), col=c(\"grey\", \"red\", \"blue\"), outline=FALSE, axes=FALSE)\n";
						print R "boxplot(c(delta_no_mut), c(delta_one), c(delta_two), boxwex=" . (int($width_boxplot/3) - 10) . ", add=TRUE, at=c(" . ($number_peaks + int($add_additional/3)) . "," . ($number_peaks + (int($add_additional/3) * 2)) . "," . ($number_peaks + (int($add_additional/3) * 3)) . "), col=c(\"grey\", \"red\", \"blue\"), outline=FALSE, axes=FALSE)\n";
						print R "p_one <- t.test(delta_no_mut, delta_one)\$p.value\n";
						print R "p_two <- t.test(delta_no_mut, delta_two)\$p.value\n";
						print R "p_both <- t.test(delta_one, delta_two)\$p.value\n";
					}
				} else {
					if($cal_pvalue == 1) {
					#	print R "boxplot(c(no_mut), c(one), c(two), boxwex=" . (int($add_additional/3) - 10) . ", add=TRUE, at=c(" . ($number_peaks + int($add_additional/3)) . "," . ($number_peaks + (int($add_additional/3) * 2)) . "," . ($number_peaks + (int($add_additional/3) * 3)) . "), col=c(\"grey\", \"red\", \"blue\"), outline=FALSE, axes=FALSE)\n";
						print R "boxplot(c(no_mut), c(one), c(two), boxwex=" . (int($width_boxplot/3) - 10) . ", add=TRUE, at=c(" . ($number_peaks + int($add_additional/3)) . "," . ($number_peaks + (int($add_additional/3) * 2)) . "," . ($number_peaks + (int($add_additional/3) * 3)) . "), col=c(\"grey\", \"red\", \"blue\"), outline=FALSE, axes=FALSE)\n";
						print R "p_one <- t.test(no_mut, one)\$p.value\n";
						print R "p_two <- t.test(no_mut, two)\$p.value\n";
						print R "p_both <- t.test(one, two)\$p.value\n";
					}
				}
				if($cal_pvalue == 1) {
					print R "legend(\"bottomright\", c(\"p-values: \", paste(\"" . $strains[$i] . " vs bg: \", round(p_one, digits=6), sep=\"\"), paste(\"" . $strains[$j] . " vs bg: \", round(p_two, digits=6), sep=\"\"), paste(\"" . $strains[$i] . " vs " . $strains[$j] . ": \", round(p_both, digits=6), sep=\"\")), text.col=c(\"black\", \"red\", \"blue\", \"purple\"), cex=0.8, bty=\'n\')\n"; 
				}
			#	close R;
			}
		} 
	}
	print R "dev.off()\n";
	close R;
	
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

sub analyze_motifs{
	open FH, "<$tmp_out_main_motif";
	$last_line = "";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@header = split('_', $split[0]);
		#Read in blocks - stop after gathering all information about one position
		if($last_line ne "" && $header[0] . "_" . $header[1] . "_" . $header[2] ne $last_line) {
			my @a = split("_", $last_line);
			&print_results(substr($a[0], 3), $a[1]);
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
	#	print $line . "\n";
	}
	my @a = split("_", $last_line);
	return (substr($a[0], 3), $a[1]);
}

sub print_results {
	%analysis_result = ();
	my $curr_chr = $_[0];
	my $curr_pos = $_[1];
	#Now run through all motifs and see if they are in all the other strains
#	print "\n\n\n";
#	print Dumper %block;
#	print "\n\n";
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
		$fileHandlesMotif[$index_motifs{$motif}]->print($tag_counts{$curr_chr}{$curr_pos});
	#	for(my $i = 0; $i < @strains; $i++) {
	#		$fileHandlesMotif[$index_motifs{$motif}]->print("tg " . $strains[$i] . "\t");
	#	}
		$fileHandlesMotif[$index_motifs{$motif}]->print($fc{$curr_chr}{$curr_pos});
		for(my $i = 0; $i < @strains; $i++) {
			$fileHandlesMotif[$index_motifs{$motif}]->print($diff{$strains[$i]} . "\t");
		}
		for(my $i = 0; $i < @strains; $i++) {
		#	$fileHandlesMotif[$index_motifs{$motif}]->print($strains[$i]. " ");
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
