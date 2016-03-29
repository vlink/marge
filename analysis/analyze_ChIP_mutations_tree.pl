#!/usr/bin/perl

use strict;
use Getopt::Long;
use Storable;
use Statistics::Basic qw(:all);
use List::Util 'shuffle';
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
use general;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
#require '/Users/verenalink/workspace/strains/general/config.pm';
use analysis_tree;
use Data::Dumper;
#require '../db_part/database_interaction.pm';
use Set::IntervalTree;

print STDERR "Add pairwise comparison\n";
print STDERR "p-value calculation of pearson correlation taken out for the moment\n";

$_ = "" for my($genome, $file, $tf, $filename, $last_line, $name, $output, $ab, $plots, $overlap, $save, $load, $tmp_out, $tmp_out_no_motif, $tmp_out_far_motif, $data, $seq_no_motif, $seq_far_motif, $seq_recentered, $tmp_center);
$_ = () for my(@strains, %peaks, @split, %mutation_pos, %shift, %current_pos, %save_local_shift, %seq, %seq_far_motif, %seq_no_motif, %PWM, @fileHandlesMotif, %index_motifs, %tag_counts, %fc, @header, %block, %analysis_result, %existance, %diff, %ranked_order, %mut_one, %mut_two, %delta_score, %delete, %remove, %mut_pos_analysis, %dist_plot, %dist_plot_background, %motif_scan_scores, %all_trees, %lookup_strain, %last_strain, %tree, $seq, %peaks_recentered, %seq_recentered, @header_recenter, %recenter_conversion, $correlation, @shuffle_array, $pvalue);
$_ = 0 for my($homo, $allele, $region, $motif_score, $motif_start, $more_motifs, $save_pos, $delta, $keep, $mut_only, $tg, $filter_tg, $fc_significant, $mut_pos, $dist_plot, $effect, $center, $analyze_motif, $analyze_no_motif, $analyze_far_motif, $longest_seq_motif, $longest_seq_no_motif, $longest_seq_far_motif, $shuffle_k, $pvalue_option, $shuffle_between, $shuffle_within);
#$data = "/Users/verenalink/workspace/strains/data/";

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - Order must overlay with order in annotated peak file\n";
        print STDERR "\t-file <file>: annotated peak file (including tag counts)\n";
	print STDERR "\t-AB: Antibody that was used for this ChIP (to exclude mutations in this motif from analysis)\n";
	print STDERR "\t-center: centers peaks on TF specified in -AB\n"; 
	print STDERR "\t-TF <transcription factor motif matrix>: Matrix of the TF that was chipped for\n";
        print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-region: Size of the region used to look for other motifs (Default: 200)\n";
	print STDERR "\t-delta: Uses motif score differences instead of binary existance\n";
	print STDERR "\t-output: Name of the output files\n";
	print STDERR "\n\nAnalysis options for peaks\n";
	print STDERR "\t-motif: analyzes sequences with TF motif\n";
	print STDERR "\t-no_motif: analyzses sequences without TF motif\n";
	print STDERR "\t-far_motif: analyzses sequences with TF motif that is more thatn 25bp away from peak center\n";
	print STDERR "\t-shuffle: <number of repeats for generating bg distribution (default: 10)\n";
	print STDERR "\n\nAdditional options:\n";
	print STDERR "\t-plots: Output name of the plots\n";
	print STDERR "\t-keep: keep temporary files\n";
	print STDERR "\t-save <output name>: saves sequence files\n";
	print STDERR "\t-load <input name>: loads sequence files saved from a previous run\n";
	print STDERR "\t-data_dir <folder to data directory>: default defined in config\n";
	print STDERR "\n\nFiltering options:\n";
	print STDERR "\t-tg <minmal tag count>: Filters out all peaks with less than x tag counts\n";
	print STDERR "\t-mut_only: just keeps peaks where one strains is mutated\n";
	print STDERR "\t-fc_pos: Foldchange threshold to count peaks as strain specific vs not (Default: 2fold)\n";
	print STDERR "\t-overlap: Count motif as not mutated if the overlap n basepairs (complete|half|#bp)\n";
	print STDERR "\t-shuffle_within: Shuffles tag counts per motif to calculate significance all vs all\n";
	print STDERR "\t-shuffle_between: Shuffles tag count vectors and motif vectors to calcualte significane (default)\n";
#	print STDERR "\t-pvalue: Calculates significance all vs all based on pvalue distribution of pearson correlation (N >= 6)\n";
	print STDERR "\n\nPlot options:\n";
	print STDERR "\t-mut_pos: Also analyzes the position of the motif that is mutated\n";
	print STDERR "\t-dist_plot: Plots distance relationships between TF and motif candidates\n";
	print STDERR "\t\t-effect: Just plots distance relationship for peaks that are affected\n";
	print STDERR "Script needs R package seqinr\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}
print STDERR "change output file names for R files\n";
print STDERR "clean up the code, stop giving so much to other methods in the modules, make it more global\n";
print STDERR "CLEAN UP CODE\n";
print STDERR "ADD COMMENTS!!!!\n";

my $param = config::read_config();

#$region = 200;

GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
                "strains=s{,}" => \@strains,
                "homo" => \$homo, 
		"-TF=s" => \$tf, 
		"-region=s" => \$region,
		"-delta" => \$delta, 
#		"-pvalue" => \$pvalue_option,
		"-output=s" => \$output,
		"-plots=s" => \$plots,
		"-AB=s" => \$ab,
		"-keep" => \$keep, 
		"-data_dir=s" => \$data,
		"-tg=s" => \$tg,
		"-shuffle=s" => \$shuffle_k,
		"-shuffle_within" => \$shuffle_within,
		"-shuffle_between" => \$shuffle_between,
		"-mut_only" => \$mut_only, 
		"-mut_pos" => \$mut_pos, 
		"-overlap=s" => \$overlap,
		"-fc_pos=s" => \$fc_significant, 
		"-save=s" => \$save,
		"-load=s" => \$load, 
		"-dist_plot" => \$dist_plot, 
		"-effect" => \$effect,
		"-center" => \$center, 
		"-motif" => \$analyze_motif,
		"-no_motif" => \$analyze_no_motif,
		"-far_motif" => \$analyze_far_motif)
        or die("Error in command line options!\n");
#First step: Get the sequences for the peaks

my $reference_strain = config::reference_genome($genome);
if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}

if($shuffle_k == 0) { $shuffle_k = 10; }
if($shuffle_between == 0 && $shuffle_within == 0) {
	$shuffle_between = 1;
}

if($center == 1) {
	#Add 100bp so we can center on peak and then shorten the sequence, so we don't ahve to pull the seq twice
	$analyze_motif = 1;
} else {
	$analyze_motif = 1;
}

if($ab eq "" && $dist_plot == 1) {
	print STDERR "Distance plots not possible without anchor TF\n";
	print STDERR "Do not generate distance plots!\n";
	$dist_plot = 0;
	$effect = 0;
}

$allele = 1;
my $ref_save = 0;

if($effect == 1) { $dist_plot = 1; }

if($fc_significant == 0) {
	$fc_significant = 2;
}
#Save motif files
my ($index_motif_ref, $PWM_ref, $score_ref) = analysis::read_motifs($tf);
%index_motifs = %$index_motif_ref;
%motif_scan_scores = %$score_ref;

if($ab ne "") {
	if(!exists $index_motifs{$ab}) {
		print STDERR "TF specified in -AB does not exist in motif file\n";
		print STDERR "Please double check your motif file!\n";
		exit;
	}
}
%PWM = %$PWM_ref;

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

print STDERR "Saving peaks\n";
my $filter_out = 0;
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
	#Filter out peaks with too few tag counts
	for(my $i = 19; $i < @split; $i++) {
		if($split[$i] > $tg) {
			$filter_tg = 1;
		} 
	}	
	if($filter_tg == 0) {
		$filter_out++;
		next;
	}
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

print STDERR "" . $filter_out . " peaks filtered because of low tag counts\n\n";

print STDERR "Loading shift vectors\n";
for(my $i = 0; $i < @strains; $i++) {
	my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data, $allele, "ref_to_strain");
	$tree{$strains[$i]} = $tree_ref;
	$lookup_strain{$strains[$i]} = $lookup;
	$last_strain{$strains[$i]} = $last;
}

#Get all sequences from file and save them in hash
&get_files_from_seq();


if($center == 1) {
	&center_peaks();
} else {
	$analyze_no_motif = 0;
	$analyze_far_motif = 0;
	%seq_recentered = %{$seq};
}

my $tmp_out_main_motif = "tmp" . rand(15);

if($plots eq "") {
	$plots = "output_mut";
}

if($output eq "") {
	$output = "output_motif";
}

if($analyze_motif == 1) {
	print STDERR "Run analysis for all sequences with " . $ab . " motif\n";
	if($tmp_center ne "") {
		$tmp_out = $tmp_center . ".recentered_motif";
	}
	&screen_and_plot($output . "_with_motif", $plots . "_with_motif", $tmp_out, $tmp_out_main_motif . "_with_motif", $longest_seq_motif, \%seq_recentered);
	$delete{$tmp_out_main_motif . "_with_motif"} = 1;
}

if($analyze_no_motif == 1) {
	print STDERR "Run analysis for all sequences without " . $ab . " motif\n";
	&screen_and_plot($output . "_without_motif", $plots. "_without_motif", $tmp_center . ".no_motif", $tmp_out_main_motif . "_no_motif", $longest_seq_no_motif, \%seq_no_motif, 1);
	$delete{$tmp_out_main_motif . "_no_motif"} = 1;
}

if($analyze_far_motif == 1) {
	print STDERR "Run analysis for all sequences with " . $ab . " motif that is not in the peak center\n";
	&screen_and_plot($output . "_with_far_motif", $plots . "_with_far_motif", $tmp_center . ".far_motif", $tmp_out_main_motif . "_far_motif", $longest_seq_far_motif, \%seq_far_motif);
	$delete{$tmp_out_main_motif . "_far_motif"} = 1;
}

#Now lets split up the analysis

if($keep == 0) {
	print STDERR "Delete output files\n";
	foreach my $d (keys %delete) {
		`rm $d`
	}
}

sub get_files_from_seq{
        $tmp_out = "tmp" . rand(15);
        $delete{$tmp_out} = 1;
        print STDERR "Extracting sequences from strain genomes\n";
        my ($seq_ref, $l_seq) = analysis::get_seq_for_peaks($tmp_out, \%peaks, \@strains, $data, $allele, $line_number, $mut_only, $region, \%tree, \%lookup_strain, \%last_strain);
	$seq = $seq_ref;
	#Important for distance plots
	$longest_seq_motif = $l_seq;
}

sub screen_and_plot{
	my $output = $_[0];
	my $plots = $_[1];
	my $tmp_out = $_[2];
	my $tmp_out_main_motif = $_[3];
	my $longest_seq = $_[4];
	my $seq = $_[5];
	my $no_motif = 0;
	if(@_ > 6) {
		$no_motif = 1;
	}
	my ($fileHandlesMotif_ref, $delete_ref) = analysis::open_filehandles(\%index_motifs, \%delete, $output);
	@fileHandlesMotif = @$fileHandlesMotif_ref;
	%delete = %$delete_ref;
	analysis::scan_motif_with_homer($tmp_out, $tmp_out_main_motif, $tf);
	$delete{$tmp_out_main_motif} = 1;
	analysis::write_header(\@fileHandlesMotif, \@strains, 0);
	my ($block_ref) = analysis::analyze_motifs($tmp_out_main_motif, $seq, \@strains, \%tree, \%lookup_strain);
	my %block = %$block_ref;
	$block_ref = analysis::merge_block(\%block, $overlap, \@strains, $seq);
	%block = %$block_ref;
	if(@strains > 2) {
		($correlation, $pvalue) = analysis::all_vs_all_comparison(\%block, \%recenter_conversion, \%tag_counts, \@strains);
		if($pvalue_option == 0) {
			for(my $k = 0; $k < $shuffle_k; $k++) {
				#Randomize relationship between motif scores and tag counts
				my @shuffle;
				my $shuffle;
				my %shuffle_tag_counts;
				my $random;
				my @vector;
				if($shuffle_between == 1) {
					my $run_index = 0;
					foreach my $chr (keys %tag_counts) {
						foreach my $pos (keys %{$tag_counts{$chr}}) {
							$shuffle[$run_index] = $tag_counts{$chr}{$pos};
							$run_index++;
						}
					}
				#	print Dumper @shuffle;
					foreach my $chr (keys %tag_counts) {
						foreach my $pos (keys %{$tag_counts{$chr}}) {
							$random = int(rand(@shuffle));
							$shuffle_tag_counts{$chr}{$pos} = $shuffle[$random];
							splice(@shuffle, $random, 1);
						}
					}
				#Randomize relationship within tag count vector
				} else {
					foreach my $chr (keys %tag_counts) {
						foreach my $pos (keys %{$tag_counts{$chr}}) {
							@vector = split('\t', $tag_counts{$chr}{$pos});
							$shuffle = "";
							while(@vector > 0) {
								$random = int(rand(@vector));
								$shuffle .= $vector[$random] . "\t";
								splice(@vector, $random, 1);
							}
							chop $shuffle;
							$shuffle_tag_counts{$chr}{$pos} = $shuffle; 
						}
					}
				}
				my($shuffle_correlation, $pvalue_shuffle) = analysis::all_vs_all_comparison(\%block, \%recenter_conversion, \%shuffle_tag_counts, \@strains);
				$shuffle_array[$k] = $shuffle_correlation;
				print STDERR "round : " . $k . "\n";
				
			}
		}
		&generate_all_vs_all_R_files($plots . ".R", $output, $correlation, \@shuffle_array, $pvalue);
	} else { 
		analysis::output_motifs(\%block, \@fileHandlesMotif, \%tag_counts, \@strains, \%index_motifs, $fc_significant, $overlap, \%fc, \%recenter_conversion);
		for(my $i = 0; $i < @fileHandlesMotif; $i++) {
			close $fileHandlesMotif[$i];
		}
		print STDERR "Generating R files!\n";
		&generate_R_files($plots . ".R", $output);
		if($mut_pos == 1) {
			my ($mut_pos_analysis_ref) = analysis::analyze_motif_pos(\%block, \%fc, \@strains, $fc_significant);
			%mut_pos_analysis = %$mut_pos_analysis_ref; 
			&generate_mut_pos_analysis_file($plots . "_mut_pos_motifs.R");
		}
	}
	if($dist_plot == 1) {
		#Check background distribution
		my ($dist_plot_ref) = analysis::distance_plot(\%block, \@strains, \%fc, $fc_significant, $effect, $longest_seq, $seq);
		%dist_plot = %{$dist_plot_ref};
		my $bg_folder = config::read_config()->{'motif_background_path'};
		my ($delete_ref, $dist_plot_background_ref) = analysis::background_dist_plot($bg_folder, \%index_motifs, \%delete, \%motif_scan_scores, $genome, $ab, $region, $longest_seq);
		%delete = %{$delete_ref};
		%dist_plot_background = %{$dist_plot_background_ref};
		&generate_dist_plot($plots . "_distance.R", $no_motif, $longest_seq);	
	}

	if($delta == 1) {
		&generate_delta_files($plots . "_delta.R");
	}

	if($ab ne "") {
		if(@strains > 2) { 
			exit;
		} else {
			open FH, "<", $output . "_" . $ab . ".txt";
			foreach my $line (<FH>) {
				chomp $line;
				@split = split('\t', $line);
				if($split[-3] > 0 || $split[-2] > 0) {
					$remove{$split[1]} = 1;
				}
			}
			@fileHandlesMotif = ();
			my ($fileHandlesMotif_ref, $delete_ref) = analysis::open_filehandles(\%index_motifs, \%delete, $output . "_removed");
			%delete = %$delete_ref;
			@fileHandlesMotif = @$fileHandlesMotif_ref;
			#Remove all positions from the files with mutations in chipped motif
			foreach my $key (keys %index_motifs) {
				open FH, "<", $output . "_" . $key . ".txt";
				my $filename  = $output . "_removed_" . $key . ".txt";
				$delete{$filename} = 1;
				open my $fh, ">", $filename or die "Can't open $filename: $!\n";
				$fileHandlesMotif[$index_motifs{$key}] = $fh;
				foreach my $line (<FH>) {
					chomp $line;
					@split = split("\t", $line);
					if(!exists $remove{$split[1]}) {
						$fileHandlesMotif[$index_motifs{$key}]->print($line . "\n");
					}
				}
			}
		#	&generate_R_files($plots . "_removed.R", 0, $output . "_removed");
			&generate_R_files($plots . "_removed.R", $output . "_removed");
			print STDERR "add motif mutation plots for removed data set\n";
			foreach my $pos (keys %block) {
				if(exists $remove{$pos}) {
					delete $remove{$pos};
				}
			}
			if($mut_pos == 1) {
				my ($mut_pos_analysis_ref) = analysis::analyze_motif_pos(\%block, \%fc, \@strains, $fc_significant);
				%mut_pos_analysis = %$mut_pos_analysis_ref;
				&generate_mut_pos_analysis_file($plots . "_mut_pos_motifs_removed.R");
			}
		}
	} else {
		print STDERR "No antibody specified, analysis ends here\n";
	}
}


#Generate distance plots
sub generate_dist_plot{
	my $output = $_[0];
	my $no_motif = $_[1];
	my $longest_seq = $_[2];
	print STDERR "open R file: " . $output . "\n";
	open R, ">$output";
	my $x = "x <- c(";
	if($longest_seq % 2 == 1) {
		$longest_seq++;
	}

	print STDERR "longest seq: " . $longest_seq . "\n";
	print STDERR "region: " . $region . "\n";
	for(my $i = 0; $i < $longest_seq + $region; $i++) {
		$x .= ($i - int($longest_seq/2)) . ",";
	}
	if(length($x) > 8) {
		chop $x;
	}
	$x .= ")";
	print R "pdf(\"" . substr($output, 0, length($output) - 2) . ".pdf\", width=10, height=5)\n";
	print R $x . "\n";
	#Generate colors
	print R "color <- colorRampPalette(c(\"blue\", \"red\"))(" . (@strains) . ")\n";
	my $y = "";
	my $first = 0;
	my $legend_col;
	my $legend_class;
	my $max_ab = 0;
	my $max = 0;
	my $strain_number = 1;
	my $max_main = "max_main <- max(";
	#First generate vectors for main TF
	foreach my $strain (keys %{$dist_plot{$ab}}) {
		$y = $ab . "_" . $strain . " <- c(";
		$max_main .= $ab . "_" . $strain . ",";
		for(my $i = 0; $i < $longest_seq; $i++) {
	#	for(my $i = int(($longest_seq/2) * (-1)); $i < ($longest_seq/2); $i++) {
			if(!exists $dist_plot{$ab}{$strain}{$i}) {
				$y .= "0,";
			} else {
				$y .= $dist_plot{$ab}{$strain}{$i} . ",";
				if($dist_plot{$ab}{$strain}{$i} > $max_ab) {
					$max_ab = $dist_plot{$ab}{$strain}{$i};
				}	
			}	
		}
		if(length($y) > length($a) + length($strain) + 8) {
			chop $y;
		}
		$y .= ")";
		print R $y . "\n";
	} 
	print R $max_main . "1)\n";
	#Add background
	$y = $ab . "_background <- c(";
	for(my $i = (int($longest_seq + $region)/2)* -1; $i < int($longest_seq + $region)/2; $i++) {
		if(!exists $dist_plot_background{$ab}{$i}) {
			$y .= "0,";
		} else {
			$y .= $dist_plot_background{$ab}{$i} . ",";
		}
	}
	chop $y;
	$y .= ")";
	print R $y . "\n";

	my $second_factor;
	my $motif_print;
	my %second_factor;
	my $max_second = "max_second <- max(";
	foreach my $motif (keys %dist_plot) {
		$max_second = "max_second <- max(";
		if($motif eq $ab) { next; }
		$motif_print = $motif;
		$motif_print =~ s/\-//g;
		$motif_print =~ s/://g;
		$motif_print =~ s/\+//g;
		$max = 0;
		$first = 0;
		$strain_number = 1;
		$y = $motif_print . "_background <- c(";
		for(my $i = int(($longest_seq + $region)/2) * -1; $i < int(($longest_seq + $region)/2); $i++) {
			if(!exists $dist_plot_background{$motif}{$i}) {
				$y .= "0,";
			} else {
				$y .= $dist_plot_background{$motif}{$i} . ",";
			}
		}
		chop $y;
		print R $y . ")\n";
		
		foreach my $strain (keys %{$dist_plot{$motif}}) {
			$y = $motif_print . "_" . $strain . " <- c(";
			$max_second .= $motif_print . "_" . $strain . ",";
		#	for(my $i = int(($longest_seq/2) * (-1)); $i < ($longest_seq/2); $i++) {
			for(my $i = 0; $i < $longest_seq; $i++) {
				if(!exists $dist_plot{$motif}{$strain}{$i}) {
					$y .= "0,";
				} else {
					$y .= $dist_plot{$motif}{$strain}{$i} . ",";
					if($dist_plot{$motif}{$strain}{$i} > $max) {
						$max = $dist_plot{$motif}{$strain}{$i};
					}	
				}	
			}
			chop $y;
			$y .= ")";
			print R $y . "\n";
			if($first == 0) {
				if($no_motif == 0) {
					print R "plot(x, " . $ab . "_" . $strain . ", type=\"l\", col=\"forestgreen\", main=\"Distance plot for " . $ab . " and " . $motif_print . "\", ylim=c(0, " . $max_ab . "), ylab=\"Motif frequence\")\n";
				}
				$legend_class = "legend_class <- c(\"" . $ab . "_" . $strain. "\",";
				$legend_col = "legend_col <- c(\"forestgreen\",";
			#	if($no_motif == 0) {
					$second_factor{$strain} = "plot(x, " . $motif_print . "_" . $strain . ", col=\"red\", xaxt=\"n\",yaxt=\"n\",xlab=\"\",ylab=\"\", type=\"l\")\n";
			#	} else {
			#		$second_factor{$strain} = "plot(x, " . $motif_print . "_" . $strain . ", col=color[" . $strain_number . "], yaxt=\"n\",xlab=\"\",ylab=\"Motif frequency\", type=\"l\", main=\"Distance plot for motif " . $motif_print . " for peaks without " . $ab . "\")\n";
			#	}
				$legend_class .= "\"" . $motif_print . "_" . $strain . "\",";
				$legend_col .= "color[" . $strain_number . "],";
				$first++;
				$strain_number++;
			} else {
				if($no_motif == 0) {
					print R "lines(x, " . $ab . "_" . $strain . ", col=\"green\")\n";
			#		print R "par(new=TRUE)\n";
				}
				$legend_class .= "\"" . $ab . "_" . $strain . "\",";
				$legend_col .= "\"green\",";
				print R $second_factor . "\n";
			#	print R "axis(4)\n";
				$second_factor{$strain} =  "lines(x, " . $motif_print . "_" . $strain . ", col=color[" . $strain_number . "])\n";
				$legend_class .= "\"" . $motif_print . "_" . $strain . "\",";
				$legend_col .= "color[" . $strain_number . "],";
				$strain_number++;
			}
		}
		chop $max_second;
		if($no_motif == 0) {
			print R "lines(x, " . "((" . $ab . "_background/max(" . $ab . "_background)) * max_main), col=\"black\", lty=3)\n";
		} else {
			print R "plot(x, " . "((" . $ab . "_background/max(" . $ab . "_background)) * max_main), col=\"black\", lty=3, type=\"l\", main=\"Distance plot for motif " . $motif_print . "for peaks without " . $ab . "\", ylim=c(0, " . $max_ab . "), ylab=\"Motif frequence\")\n";


		}
		print R $max_second . ")\n";
		print R "par(new=TRUE)\n";
		foreach my $strain (keys %{$dist_plot{$motif}}) {
			print R $second_factor{$strain};
		}
		print R "lines(x, " . "((" . $motif_print . "_background/max(" . $motif_print . "_background))*max_second), col=\"darkgrey\", lty=3)\n";
		print R "axis(4)\n";
	#	chop $legend_class;
	#	chop $legend_col;	
		print R $legend_class . "\"" . $ab . " background\", \"" . $motif_print . " background\")\n";
		print R $legend_col . "\"black\", \"darkgrey\")\n";
		print R  "legend(\"topright\", legend_class, col=legend_col, pch=20, bty=\'n\')\n";
	}
	print R "dev.off()\n";
	close R;
	`Rscript $output`;	
}

#Gnereate motif position mutation plots
sub generate_mut_pos_analysis_file{
	my $output = $_[0];
	my $divide = 0;
	open R, ">$output";
	my $run = 0;
	#Step one: generate logo sequence for motif
	print R "library(seqLogo)\n";
	print R "library(grid)\n";
	print R "library(gridBase)\n";
	print R "mySeqLogo = seqLogo::seqLogo\n";
	print R "bad = (sapply( body(mySeqLogo), \"==\", \"grid.newpage()\") | sapply( body(mySeqLogo), \"==\", \"par(ask=FALSE)\"))\n";
	print R "body(mySeqLogo)[bad] = NULL\n";
	print R "pdf(\"" . substr($output, 0, length($output) - 2) . ".pdf\", width=10, height=5)\n";
	my %col;
	$col{'A'} = "green";
	$col{'C'} = "blue";
	$col{'G'} = "orange";
	$col{'T'} = "red";
	my $first = 0;
	$_ = "" for my($mut_freq_sig, $mut_freq_unsig, $y_sig, $y_unsig);
	foreach my $motif (keys %mut_pos_analysis) {
		my $A = "A <- c(";
		my $C = "C <- c(";
		my $G = "G <- c(";
		my $T = "T <- c(";
		#Check if PWM is smaller 1
		foreach my $pos (sort {$a <=> $b} keys %{$PWM{$motif}}) {
			if($PWM{$motif}{$pos}{'A'} > 1 || $PWM{$motif}{$pos}{'C'} > 1 || $PWM{$motif}{$pos}{'G'} > 1 || $PWM{$motif}{$pos}{'T'} > 1) {
				$divide = ($PWM{$motif}{$pos}{'A'} + $PWM{$motif}{$pos}{'C'} + $PWM{$motif}{$pos}{'G'} + $PWM{$motif}{$pos}{'T'});
				$A .= "" . ($PWM{$motif}{$pos}{'A'}/$divide) . ",";
				$C .= "" . ($PWM{$motif}{$pos}{'C'}/$divide) . ",";
				$G .= "" . ($PWM{$motif}{$pos}{'G'}/$divide) . ",";
				$T .= "" . ($PWM{$motif}{$pos}{'T'}/$divide) . ",";
			} else {
				$A .= $PWM{$motif}{$pos}{'A'} . ",";
				$C .= $PWM{$motif}{$pos}{'C'} . ",";
				$G .= $PWM{$motif}{$pos}{'G'} . ",";
				$T .= $PWM{$motif}{$pos}{'T'} . ",";
			}
		}
		chop $A;
		chop $C;
		chop $G;
		chop $T;
		$A .= ")";
		$C .= ")";
		$G .= ")";
		$T .= ")";
		$A =~ s/,+\)/\)/g;
		$C =~ s/,+\)/\)/g;
		$G =~ s/,+\)/\)/g;
		$T =~ s/,+\)/\)/g;
		print R $A . "\n";
		print R $C . "\n";
		print R $G . "\n";
		print R $T . "\n";
		print R "pwm = data.frame(A, C, G, T)\n";
		print R "pwm = t(pwm)\n";
		print R "par(mar=c(2.5,2.5,1,1), oma=c(2.5,2.5,1,1))\n";
		$first = 0;
		my $max = 0;
		my $count = 0;
		my $step = (1/(keys %{$PWM{$motif}}));
		foreach my $base (keys %col) {
			$max = 0;
			$mut_freq_sig = "mut_freq_sig_" . $base . " <- c(";
			$mut_freq_unsig = "mut_freq_unsig_" . $base . " <- c(";
			$y_sig = "y_sig <- c(";
			$y_unsig = "y_unsig <- c(";
			foreach my $pos (sort {$a <=> $b} keys %{$PWM{$motif}}) {
				if(!exists $mut_pos_analysis{$motif}{$base}{$pos}{'S'}) {
					$mut_freq_sig .= "0,";
				} else {
					$mut_freq_sig .= $mut_pos_analysis{$motif}{$base}{$pos}{'S'} . ",";
					if($mut_pos_analysis{$motif}{$base}{$pos}{'S'} > $max) {
						$max = $mut_pos_analysis{$motif}{$base}{$pos}{'S'};
					}
				}
				if(!exists $mut_pos_analysis{$motif}{$base}{$pos}{'N'}) {
					$mut_freq_unsig .= "0,";
				} else {
					$mut_freq_unsig .= $mut_pos_analysis{$motif}{$base}{$pos}{'N'} . ",";
					if($mut_pos_analysis{$motif}{$base}{$pos}{'N'} > $max) {
						$max = $mut_pos_analysis{$motif}{$base}{$pos}{'N'};
					}
				}
				if($pos > 1) {
					$y_sig .= $pos . " + ($step  * ($pos - 1)) + (0.05 * $count),";
					$y_unsig .= $pos . " + ($step * ($pos - 1)) + 0.2 + (0.05 * $count),";
				} else {
					$y_sig .= $pos . " + (0.05 * $count),";
					$y_unsig .= $pos . " + 0.2 + (0.05 * $count),";

				}
			}
			$count++;
			if(length($mut_freq_sig) > 20) {
				chop $mut_freq_sig;
			} else {
				$mut_freq_sig .= "0";
			}
			if(length($mut_freq_unsig) > 22) {
				chop $mut_freq_unsig;
			} else {
				$mut_freq_unsig .= "0";
			}
			chop $y_sig;
			chop $y_unsig;
			$mut_freq_sig .= ")";
			$mut_freq_unsig .= ")";
			$y_sig .= ")";
			$y_unsig .= ")";
			print R $mut_freq_sig . "\n";
			print R $mut_freq_unsig . "\n";
			print R $y_sig . "\n";
			print R $y_unsig . "\n";
			if($first == 0) { 
				print R "plot(y_sig, mut_freq_sig_" . $base . ", xlab=NA, ylim=c(-0.5, " . $max . " + 1), axes=FALSE, main=\"" . $motif . "\", col=\"" . $col{$base} . "\", pch=20, xlim=c(0, ". (keys %{$PWM{$motif}}) . "))\n";
				$first++;
			} else {
				print R "points(y_sig, mut_freq_sig_" . $base . ", col=\"" . $col{$base} . "\", pch=20)\n";
			}
			print R "points(y_unsig, mut_freq_unsig_" . $base . ", col=\"" . $col{$base} . "\", pch=8)\n";
			print R "axis(2)\n";
		}
		print R "legend(\"topleft\", c(\"sig\", \"unsig\", \"A\", \"C\", \"G\", \"T\"), col=c(\"black\", \"black\", \"" . $col{'A'} . "\", \"" . $col{'C'} . "\", \"" . $col{'G'} . "\", \"" . $col{'T'} . "\"), pch=c(20, 8, 16, 16, 16, 16))\n";
		print R "opar <- par(las=1)\n";
		print R "par(opar)\n";	
		print R "mtext(\"Frequencies\", 2, 3)\n";
		print R "vp1 <- viewport(x=0.03, y=0, width=1, height=0.3, just=c(\"left\", \"bottom\"))\n";
		print R "pushViewport(vp1)\n";
		print R "par(new=TRUE, mar=c(2.5,1,1,1.5), oma=c(2.5,1,1,1.5))\n";
		print R "mySeqLogo(pwm, xaxis=FALSE, yaxis=FALSE, ic.scale=FALSE)\n"; 
		print R "popViewport()\n";
		print R "grid.newpage()\n";
		$run++;
	}
	print R "dev.off()\n";
	close R;
	`Rscript $output`;
}

#Generate R files
sub generate_R_files {
#	my $output = $_[2];
	my $output = $_[1];
	my $filename;
	my $f = 0;
	my $num_of_muts;
	my $output_R_file = $_[0];
#	my $removed = $_[1];
	print STDERR "Add filtering for number of motifs!\n";
	open R, ">", $output_R_file;
	print R "pdf(\"" . substr($output_R_file, 0, length($output_R_file) - 2) . ".pdf\", width=10, height=5)\n";
	my $cal_pvalue = 0;
	my $count_obs_one = 0;
	my $count_obs_two = 0;
	my $count_obs_both = 0;
	foreach my $motif (sort {$index_motifs{$a} cmp $index_motifs{$b}} keys %index_motifs) {
		$filename = $output . "_" . $motif . ".txt";
		$num_of_muts = `wc -l $filename`;
		@split = split('\s+', $num_of_muts);
		if($split[0] == 1) {
			print STDERR "No occurrences of " . $motif . " found close to " . $ab . "\n";
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
					}
					if($split[$start_pos + $j] eq "1") {
						$mut_two{$strains[$i] . "_" . $strains[$j]}{$split[1]} = 1;
					}
					$start_pos = (2 + @strains + $fac);
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
					$max = $max + $max_delta;
					$min = $min + $min_delta;
				}
				$cal_pvalue = 0;
				my $add_additional =  ((keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) * 0.25);
				if($add_additional < 20) { $add_additional = 20; }
				my $width_boxplot = $add_additional;
				if($add_additional < 31) { $width_boxplot = 45; }
				my $number_peaks = (keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}});
				print R "plot(c(0," . ((keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) + $add_additional) . "), c(" . $min . ", " . $max . "+1), col=\"white\", xlab=\"rank-ordered peaks\", ylab=\"FC " . $strains[$i] . " vs " . $strains[$j] . "\", main=\"" . $strains[$i] . " vs " . $strains[$j] . "\\n$motif\")\n"; 
				print R "abline(c(0,0), c(0,0))\n";
				#Sort by fold change value between strains
				foreach my $rank (sort {$ranked_order{$strains[$i] . "_" . $strains[$j]}{$a} <=> $ranked_order{$strains[$i] . "_" . $strains[$j]}{$b}} keys %{$ranked_order{$strains[$i] . "_" . $strains[$j]}}) {
					if(exists $mut_one{$strains[$i] . "_" . $strains[$j]}{$rank}) {
						$dist_one .= sprintf("%.2f", $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank}) . ",";
						if($delta == 1) {
							print R "points(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", pch=8, col=\"red\")\n";
						 	print R "segments(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", x1 = " . $x_count . ", y1= " . ($ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} + $delta_score{$strains[$i] . "_" . $strains[$j]}{$rank}) . ", col=\"red\")\n";
						} else {
							print R "segments(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . " - 0.5" . ", x1 = " . $x_count . ", y1= " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . "+0.5, col=\"red\")\n";
						}
					}
					if(exists $mut_two{$strains[$i] . "_" . $strains[$j]}{$rank}) {
						$dist_two .= sprintf("%.2f", $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank}) . ",";
						if($delta == 1) {
							print R "points(" . $x_count . ", " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", pch=8, col=\"blue\")\n";
							print R "segments(" . $x_count . " + 0.1, " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", x1 = " . $x_count . " + 0.1" . ", y1= " . ($ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} + $delta_score{$strains[$i] . "_" . $strains[$j]}{$rank}) . ", col=\"blue\")\n";
						} else {
							print R "segments(" . $x_count . " + 0.1, " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . ", x1 = " . $x_count . " + 0.1" . ", y1= " . $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank} . "+1, col=\"blue\")\n";
						}
					}
					if(!exists $mut_one{$strains[$i] . "_" . $strains[$j]}{$rank} && !exists $mut_two{$strains[$i] . "_" . $strains[$j]}{$rank}) {
						$dist_no_mut .= sprintf("%.2f", $ranked_order{$strains[$i] . "_" . $strains[$j]}{$rank}) . ",";
					}
					$x_count++;	
				}
				if(length($dist_one) > 10) {
					chop $dist_one;
				}
				if(length($dist_two) > 10) {
					chop $dist_two;
				}
				if(length($dist_no_mut) > 15) {
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
				print R "legend(\"topleft\", c(paste(\"" . $strains[$i] . ": \", length(one), \" muts\", sep=\"\"), paste(\"" . $strains[$j] . ": \", length(two), \" muts\", sep=\"\")), col=c(\"red\", \"blue\"), pch=16)\n";
				if($cal_pvalue == 1) {
					print R "boxplot(c(no_mut), c(one), c(two), boxwex=" . (int($width_boxplot/3) - 10) . ", add=TRUE, at=c(" . ($number_peaks + int($add_additional/3)) . "," . ($number_peaks + (int($add_additional/3) * 2)) . "," . ($number_peaks + (int($add_additional/3) * 3)) . "), col=c(\"grey\", \"red\", \"blue\"), outline=FALSE, axes=FALSE)\n";
					print R "p_one <- t.test(no_mut, one)\$p.value\n";
					print R "p_two <- t.test(no_mut, two)\$p.value\n";
					print R "p_both <- t.test(one, two)\$p.value\n";
					print R "legend(\"bottomright\", c(\"p-values: \", paste(\"" . $strains[$i] . " vs bg: \", round(p_one, digits=6), sep=\"\"), paste(\"" . $strains[$j] . " vs bg: \", round(p_two, digits=6), sep=\"\"), paste(\"" . $strains[$i] . " vs " . $strains[$j] . ": \", round(p_both, digits=6), sep=\"\")), text.col=c(\"black\", \"red\", \"blue\", \"purple\"), cex=0.8, bty=\'n\')\n"; 
				}
			}
		} 
	}
	print R "dev.off()\n";
	close R;
	`Rscript $output_R_file`;	
}

sub generate_delta_files{
	my $output_R = $_[0];
	open R, ">$output_R";
	$_ = "" for my($x, $y, $f);
	print R "pdf(\"" . substr($output_R, 0, length($output_R) - 2) . ".pdf\", width=10, height=5)\n";
	foreach my $motif (sort {$index_motifs{$a} cmp $index_motifs{$b}} keys %index_motifs) {
		#Make sure there a mutations in this motif
		$x = "x <- c(";
		$y = "y <- c(";
		$filename = $output . "_" . $motif . ".txt";
		open FH, "<$filename" or die "Can't find $filename: $!\n";
		$f = 0;
		foreach my $line (<FH>) {
			if($f == 0) { $f++; next;}
			chomp $line;
			@split = split('\t', $line);
			$x .= $split[4] . ",";
			$y .= ($split[5] - $split[6]) . ",";
		}
		if(length($x) > 8) {
			chop $x;
		}
		$x .= ")";
		if(length($y) > 8) {
			chop $y;
		}
		$y .= ")";
		print R $x . "\n";
		print R $y . "\n";
		print R "plot(x, y, main=\"" . $motif . "\ motif diff vs fc\", xlim=c(min(x, y, 0), max(x, y, 1)), ylim=c(min(x, y, 0), max(x, y, 1)), xlab=\"Foldchange " . $strains[0] . " vs " . $strains[1] . "\", ylab=\"motif difference " . $strains[0] . " vs " . $strains[1] . "\", pch=20)\n";
		print R "abline(c(0,0), c(1,1))\n";
			
	}
	print R "dev.off()\n";
	close R;
	`Rscript $output`;
}

sub center_peaks{
	print STDERR "Center the peaks!\n";
	#100 was added to grab the initial sequences
	if($ab eq "") {
		print STDERR "Not possible to center peak - no transcription factor is specified\n";
		print STDERR "Skipping...\n";
	} else {
		open TMP_MOTIF, ">tmp_scan_motif_$ab.txt";
		$delete{"tmp_scan_motif_$ab.txt"} = 1;
		print TMP_MOTIF ">consensus\t$ab\t$motif_scan_scores{$ab}\n";
		foreach my $pos (sort {$a <=> $b} keys %{$PWM{$ab}}) {
			print TMP_MOTIF $PWM{$ab}{$pos}{'A'} . "\t" . $PWM{$ab}{$pos}{'C'} . "\t" . $PWM{$ab}{$pos}{'G'} . "\t" . $PWM{$ab}{$pos}{'T'} . "\n";
		}
		close TMP_MOTIF;
	}
	my $command = "homer2 find -i " . $tmp_out . " -m tmp_scan_motif_" . $ab . ".txt -offset 0 > output_tmp.txt";
	`$command`;
	open FH, "<output_tmp.txt";
	my ($block_ref) = analysis::analyze_motifs("output_tmp.txt", \%save_local_shift, \@strains, \%tree, \%lookup_strain, \%last_strain, $allele, $region);
	%block = %$block_ref;
	my ($block_ref) = analysis::merge_block(\%block, $overlap, \@strains, \%seq);
	%block = %$block_ref;
	my $peak_center;
	my $dist_to_center;
	my $closest_motif;
	my $start;
	my $tmp_header;
	my $length;
	#Now recenter peaks
	$tmp_center = "tmp" . rand(15);
	open OUT, ">$tmp_center.far_motif";
	$delete{$tmp_center . ".far_motif"} = 1;
	foreach my $position (keys %block) {
		$peak_center = 10000;
		foreach my $motif (keys $block{$position}) {
			foreach my $pos (keys $block{$position}{$motif}) {
				for(my $i = 0; $i < @strains; $i++) {
					$dist_to_center = int(length($seq->{$position . "_" . $strains[$i]})/2 - ($pos + ($block{$position}{$motif}{$pos}{'length'}/2)));
					if(abs($dist_to_center) < abs($peak_center)) {
						$peak_center = $dist_to_center;
						$closest_motif = $pos;
					}
				}
			}
		}		
		#Split them up in no motif at all, far motif, or motif in middle and center
		if(abs($peak_center) < 25) {
			@header_recenter = split("_", $position);
			$peaks_recentered{substr($header_recenter[0], 3)}{$header_recenter[1] - $peak_center} = ($header_recenter[2] - $peak_center);
			$recenter_conversion{substr($header_recenter[0], 3)}{$header_recenter[1] - $peak_center}{'start'} = $header_recenter[1];
			$recenter_conversion{substr($header_recenter[0], 3)}{$header_recenter[1] - $peak_center}{'end'} = $header_recenter[2];
			for(my $i = 0; $i < @strains; $i++) {
				delete $seq->{$position . "_" . $strains[$i]};
			}
		} else {
			for(my $i = 0; $i < @strains; $i++) {
				print OUT ">" . $position . "_" . $strains[$i] . "\n";
				print OUT $seq->{$position . "_" . $strains[$i]} . "\n";
				$seq_far_motif{$position . "_" . $strains[$i]} = $seq->{$position . "_" . $strains[$i]};
				if(length($seq_far_motif{$position . "_" . $strains[$i]}) > $longest_seq_far_motif) {
					$longest_seq_far_motif = length($seq_far_motif{$position . "_" . $strains[$i]});
				}
				delete $seq->{$position . "_" . $strains[$i]};
			}
		}
	}
	close OUT;

	open OUT, ">$tmp_center.no_motif";
	$delete{$tmp_center . ".no_motif"} = 1;
	foreach my $position (keys %{$seq}) {
		print OUT ">" . $position . "\n";
		print OUT $seq->{$position} . "\n";
		$seq_no_motif{$position} = $seq->{$position};
		if(length($seq_no_motif{$position}) > $longest_seq_no_motif) {
			$longest_seq_no_motif = length($seq_no_motif{$position});
		}
		delete $seq->{$position};
	}
	close OUT;
	#Time to get recentered peaks
	my ($seq_ref, $l_seq) = analysis::get_seq_for_peaks($tmp_out . ".recenter", \%peaks_recentered, \@strains, $data, $allele, $line_number, $mut_only, $region, \%tree, \%lookup_strain, \%last_strain);
	$delete{$tmp_out . ".recenter"} = 1;
	%seq_recentered = %{$seq_ref};
	$longest_seq_motif = $l_seq;
	#Write sequences to files 
	$delete{$tmp_center . ".recentered_motif"} = 1;
	open OUT, ">$tmp_center.recentered_motif";
	foreach my $header (keys %seq_recentered) {
		print OUT ">" . $header . "\n";
		print OUT $seq_recentered{$header} . "\n";
	}
	close OUT;
}

sub write_seqs{
	my $file = $_[0];
	my %hash = %{$_[1]};
	my @strains = @{$_[2]};
	open OUT, ">$file";
	my %seen = ();
	my @split;
	my $tmp_header;
	foreach my $seq (sort {$a <=> $b} keys %hash) {
		@split = split("_", $seq);
		$tmp_header = $split[0] . "_" . $split[1] . "_" . $split[2];
		if(exists $seen{$tmp_header}) {
			next;
		}	
		for(my $i = 0; $i < @strains; $i++) {	
			print OUT ">" . $tmp_header . "_" . $strains[$i] . "\n";
			print OUT $hash{$tmp_header . "_" . $strains[$i]} . "\n";
		}
		$seen{$tmp_header} = 1;
	}
        close OUT;
}

sub generate_all_vs_all_R_files{
	my $plots = $_[0];
	print "plots: " . $plots . "\n";
	my $output = $_[1] . "_all_vs_all.pdf";
	my $correlation = $_[2];
	my $correlation_shuffle = $_[3];
	my $pvalue = $_[4];
	$_ = "" for my ($x, $y);
	open OUT, ">", $plots;
	print OUT "pdf(\"" . $output . "\", width=5, height=5)\n";
	print OUT "breaks <- seq(-1, 1, by=0.1)\n";
	print OUT "par(oma=c(0,0,0,0))\n";
	foreach my $motif (keys %{$correlation}) {
		my $ks = "ks <- c(";
		my $x = "x <- c(";
		my $y = "y <- c(";
		foreach my $pos (sort {$a <=> $b} keys %{$correlation->{$motif}}) {
			$x .= $pos . ",";
			$y .= $correlation->{$motif}->{$pos} . ",";
			for(my $i = 0; $i < $correlation->{$motif}->{$pos}; $i++) {
				$ks .= $pos . ",";
			}
		}
		if(length($x) > 10) {
			chop $x;
			chop $ks;
		}
		if(length($y) > 10) {
			chop $y;
		}
		print OUT $x . ")\n";
		print OUT $ks . ")\n";
		print OUT $y . ")\n";
		if($pvalue_option == 0) {
			for(my $k = 0; $k < $shuffle_k; $k++) {
				my $x_shuffle = "x_shuffle_" . $k . " <- c(";
				my $y_shuffle = "y_shuffle_" . $k . " <- c(";
				my $ks_shuffle = "ks_shuffle_" . $k . " <- c(";
				foreach my $pos (sort {$a <=> $b} keys %{$correlation_shuffle->[$k]->{$motif}}) {
					$x_shuffle .= $pos . ",";
					$y_shuffle .= $correlation_shuffle->[$k]->{$motif}->{$pos} . ",";
					for(my $i = 0; $i < $correlation_shuffle->[$k]->{$motif}->{$pos}; $i++) {
						$ks_shuffle .= $pos . ",";
					}
				}
				if(length($x_shuffle) > 10) {
					chop $x_shuffle;
					chop $ks_shuffle;
				}	
				if(length($y_shuffle) > 10) {
					chop $y_shuffle;
				}
				print OUT $x_shuffle . ")\n";
				print OUT $ks_shuffle . ")\n";
				print OUT $y_shuffle . ")\n";
				print OUT "ks_test_" . $k . " <- ks.test(ks, ks_shuffle_" . $k . ")\n";
			}	
			print OUT "p_value_avg <- mean(ks_test_0\$p.value";
			for(my $k = 1; $k < $shuffle_k; $k++) {
				print OUT ", ks_test_" . $k . "\$p.value";
			}
			print OUT ")\n";
			print OUT "print(paste(\"Motif: \", \"" . $motif . "\", \"   p value: \", p_value_avg, sep=\"\"))\n";
			print OUT "plot(x, y, main=paste(\"" . $motif . "\\nP-value: \", p_value_avg, sep=\"\"), type=\"l\", xlab=\"Pearson Correlation\", ylab=\"Frequency\", ylim=c(0, max(y, y_shuffle_0)))\n";
		} else {
			my $pvalues_list = "pvalues <- c(";
			my $pvalues_uni = "pvalues_uni <- c(";
			foreach my $values (sort {$a <=> $b} keys %{$pvalue->{$motif}}) {
				for(my $i = 0; $i < $pvalue->{$motif}->{$values}; $i++) {
					$pvalues_list .= $values . ",";
					$pvalues_uni .= rand(1) . ",";
				}
			}
			if(length($pvalues_list) > 15) {
				chop $pvalues_list;
				chop $pvalues_uni;
			}
			print OUT $pvalues_list . ")\n";
			print OUT $pvalues_uni . ")\n";
			print OUT "pvalues_uni_sort <- sort(pvalues_uni)\n";
			print OUT "ks_test <- ks.test(pvalues, pvalues_uni_sort)\n";
			print OUT "p_value_avg <- ks_test\$p.value\n";
			print OUT "print(paste(\"Motif: \", \"" . $motif . "\", \"\\tp value: \", p_value_avg, sep=\"\"))\n";
			print OUT "plot(x, y, main=paste(\"" . $motif . "\\nP-value: \", p_value_avg, sep=\"\"), type=\"l\", xlab=\"Pearson Correlation\", ylab=\"Frequency\")\n";
		}
	#	print OUT "ks_test <- ks.test(ks, ks_shuffle)\n";
		if($pvalue_option == 0) {
			for(my $k = 0; $k < $shuffle_k; $k++) {
				print OUT "lines(x_shuffle_" . $k . ", y_shuffle_" . $k . ", col=\"lightgrey\", lty=3)\n";
			}
		}
		print OUT "lines(x, y, lwd=2, col=\"black\")\n";
	}
	print OUT "dev.off()\n";
	close OUT;
	`Rscript $plots 2> /dev/null`;

}
sub fakrek{
        my $number = shift;
        return undef if $number < 0;
        return 1 if $number == 0;
        return ($number * &fakrek ($number -1));
}
