#!/usr/bin/perl


use strict;
use Getopt::Long;
#require '/Users/verenalink/workspace/strains/general/config.pm';
require '../general/config.pm';
require 'analysis.pm';
use Data::Dumper;
#require '../db_part/database_interaction.pm';


$_ = "" for my($genome, $file, $tf, $filename, $last_line, $name);
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
my ($index_motif_ref, $PWM_ref, $fileHandlesMotif_ref) = analysis::read_motifs($tf);
%index_motifs = %$index_motif_ref;
%PWM = %$PWM_ref;
@fileHandlesMotif = @$fileHandlesMotif_ref;

#Make sure the reference is in the strains array
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

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
}
close FH;

my $tmp_out = "tmp" . rand(15);
print STDERR "Extracting sequences from strain genomes\n";
print STDERR $tmp_out . "\n";
my ($seq_ref, $save_local_shift_ref) = analysis::get_seq_for_peaks($tmp_out, \%peaks, \@strains, $data, $allele, $line_number);
%seq = %$seq_ref;
%save_local_shift = %$save_local_shift_ref;

my $tmp_out_main_motif = "tmp" . rand(15);
analysis::scan_motif_with_homer($tmp_out, $tmp_out_main_motif, $tf);

analysis::write_header(\@fileHandlesMotif, \@strains, 0);

analysis::analyze_motifs($tmp_out_main_motif, \%seq, \%save_local_shift, \@fileHandlesMotif, \%tag_counts, \@strains, \%index_motifs, \%fc);

for(my $i = 0; $i < @fileHandlesMotif; $i++) {
        close $fileHandlesMotif[$i];
}

print STDERR "Generating R files!\n";
&generate_R_files();


#Generate R files
sub generate_R_files {
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
		#Make sure there a mutations in this motif
		$filename = "output_mutation_" . $motif . ".txt";
		$num_of_muts = `wc -l $filename`;
		@split = split('\s+', $num_of_muts);
		if($split[0] == 1) {
			print "No mutations in motif for " . $motif . "\n";
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
		my $start_pos = 0;
		my $fac = 0;
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
					#Save the foldchange between the strains
					$ranked_order{$strains[$i] . "_" . $strains[$j]}{$split[1]} = $split[1 + @strains + $i + $j]; 
					if($split[1 + @strains + $i + $j] > $max) {
						$max = $split[1 + @strains + $i + $j];
					}
					if($split[1 + @strains + $i + $j] < $min) {
						$min = $split[1 + @strains + $i + $j];
					}
					$fac = &fakrek(@strains - 1); 
					$start_pos = (2 + (@strains * 2) + $fac);
					if($split[$start_pos + $i] eq "1") {
						$mut_one{$strains[$i] . "_" . $strains[$j]}{$split[1]} = 1;
				#		print "mut one\n";
					}
					if($split[$start_pos + $j] eq "1") {
						$mut_two{$strains[$i] . "_" . $strains[$j]}{$split[1]} = 1;
				#		print "mut two\n";
					}
					$start_pos = (2 + @strains + $fac);
				#	print "i: " . $split[$start_pos + $i] . "\tj: " . $split[$start_pos + $j] . "\n";
				#	print "delta score: " . ($split[$start_pos + $i] - $split[$start_pos + $j]) . "\n";
					$delta_score{$strains[$i] . "_" . $strains[$j]}{$split[1]} = $split[$start_pos + $i] - $split[$start_pos + $j];
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
				print R "plot(c(0," . ((keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) + $add_additional) . "), c(" . $min . ", " . $max . "+1), col=\"white\", xlab=\"rank-ordered peaks\", ylab=\"FC " . $strains[$i] . " vs " . $strains[$j] . "\", main=\"" . $strains[$i] . " vs " . $strains[$j] . "\\n$motif\")\n"; 
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

sub fakrek{
        my $number = shift;
        return undef if $number < 0;
        return 1 if $number == 0;
        return ($number * &fakrek ($number -1));
}
