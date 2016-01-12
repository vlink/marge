#!/usr/bin/perl


use strict;
use Getopt::Long;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
#require '/Users/verenalink/workspace/strains/general/config.pm';
use analysis;
use Data::Dumper;
#require '../db_part/database_interaction.pm';


$_ = "" for my($genome, $file, $tf, $filename, $last_line, $name, $output, $ab, $seqfile, $plots, $overlap);
$_ = () for my(@strains, %peaks, @split, %mutation_pos, %shift, %current_pos, %save_local_shift, %seq, %PWM, @fileHandlesMotif, %index_motifs, %tag_counts, %fc, @header, %block, %analysis_result, %existance, %diff, %ranked_order, %mut_one, %mut_two, %delta_score, %delete, %remove, %mut_pos_analysis);
$_ = 0 for my($homo, $allele, $region, $motif_score, $motif_start, $more_motifs, $save_pos, $delta, $keep, $mut_only, $tg, $filter_tg, $fc_significant, $mut_pos);
my $data = config::read_config()->{'data_folder'};
#$data = "/Users/verenalink/workspace/strains/data/";

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - Order must overlay with order in annotated peak file\n";
        print STDERR "\t-file <file>: annotated peak file (including tag counts)\n";
	print STDERR "\t-seqfile <file>: File with sequences for peaks\n";
	print STDERR "\t-AB: Antibody that was used for this ChIP (to exclude mutations in this motif from analysis)\n";
	print STDERR "\t-TF <transcription factor motif matrix>: Matrix of the TF that was chipped for\n";
        print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-region: Size of the region used to look for other motifs (Default: 200)\n";
	print STDERR "\t-delta: Uses motif score differences instead of binary existance\n";
	print STDERR "\t-output: Name of the output files\n";
	print STDERR "\t-plots: Output name of the plots\n";
	print STDERR "\t-keep: keep temporary files\n";
	print STDERR "\n\nFiltering options:\n";
	print STDERR "\t-tg <minmal tag count>: Filters out all peaks with less than x tag counts\n";
	print STDERR "\t-mut_only: just keeps peaks where one strains is mutated\n";
	print STDERR "\t-mut_pos: Also analyzes the position of the motif that is mutated\n";
	print STDERR "\t-fc_pos: Foldchange threshold to count peaks as strain specific vs not (Default: 2fold)\n";
	print STDERR "\t-overlap: Count motif as not mutated if the overlap n basepairs (complete|half|#bp)\n";
	print STDERR "Script needs R package seqinr\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $param = config::read_config();

$region = 200;

$output = "output_motif_" . rand(5);
GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
		"seqfile=s" => \$seqfile,
                "strains=s{,}" => \@strains,
                "homo" => \$homo, 
		"-TF=s" => \$tf, 
		"-region=s" => \$region,
		"-delta" => \$delta, 
		"-output=s" => \$output,
		"-plots=s" => \$plots,
		"-AB=s" => \$ab,
		"-keep" => \$keep, 
		"-tg=s" => \$tg,
		"-mut_only" => \$mut_only, 
		"-mut_pos" => \$mut_pos, 
		"-overlap=s" => \$overlap,
		"-fc_pos=s" => \$fc_significant)
        or die("Error in command line options!\n");
#First step: Get the sequences for the peaks
$allele = 1;
my $ref_save = 0;

if($fc_significant == 0) {
	$fc_significant = 2;
}
#Save motif files
my ($index_motif_ref, $PWM_ref, $fileHandlesMotif_ref, $del_ref) = analysis::read_motifs($tf, $output, \%delete);
%index_motifs = %$index_motif_ref;
%PWM = %$PWM_ref;
@fileHandlesMotif = @$fileHandlesMotif_ref;
%delete = %$del_ref;

#Make sure the reference is in the strains array
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
my $tmp_out = "tmp" . rand(15);
$delete{$tmp_out} = 1;

print STDERR "Extracting sequences from strain genomes\n";
print STDERR $tmp_out . "\n";
my ($seq_ref, $save_local_shift_ref) = analysis::get_seq_for_peaks($tmp_out, \%peaks, \@strains, $data, $allele, $line_number, $mut_only);
%seq = %$seq_ref;
%save_local_shift = %$save_local_shift_ref;

my $tmp_out_main_motif = "tmp" . rand(15);
$delete{$tmp_out_main_motif} = 1;
analysis::scan_motif_with_homer($tmp_out, $tmp_out_main_motif, $tf);

analysis::write_header(\@fileHandlesMotif, \@strains, 0);

my ($remove_ref, $mut_pos_analysis_ref) = analysis::analyze_motifs($tmp_out_main_motif, \%seq, \%save_local_shift, \@fileHandlesMotif, \%tag_counts, \@strains, \%index_motifs, $ab, \%remove, $fc_significant, $mut_pos, $overlap, \%fc);
%remove = %$remove_ref;
%mut_pos_analysis = %$mut_pos_analysis_ref;

for(my $i = 0; $i < @fileHandlesMotif; $i++) {
        close $fileHandlesMotif[$i];
}

if($plots eq "") {
	$plots = "output_all_motifs";
}
print STDERR "Generating R files!\n";
&generate_R_files($plots . ".R", 1);
&generate_mut_pos_analysis_file($plots . "_mut_pos_motifs.R");
if($delta == 1) {
	&generate_delta_files($plots . "_delta.R");
}

print STDERR "change output file names for R files\n";
print STDERR "clean up the code, stop giving so much to other methods in the modules, make it more global\n";
print STDERR "CLEAN UP CODE\n";
print STDERR "ADD COMMENTS!!!!\n";


if($ab ne "") {
	@fileHandlesMotif = ();

	#Remove all positions from the files with mutations in chipped motif
	foreach my $key (keys %index_motifs) {
		open FH, "<", $output . "_" . $key . ".txt";
		my $filename  = $output . "_" . $key . "_removed.txt";
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
	&generate_R_files($plots . "_removed.R", 0);
	print STDERR "add motif mutation plots for removed data set\n";
	&generate_mut_pos_analysis_file($plots . "_mut_pos_motifs_removed.R");
} else {
	print STDERR "No antibody specified, analysis ends here\n";
}

if($keep == 0) {
	print STDERR "Delete output files\n";
	foreach my $d (keys %delete) {
		`rm $d`;
	}
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
	my $filename;
	my $f = 0;
	my $num_of_muts;
	my $output_R_file = $_[0];
	my $removed = $_[1];
	print STDERR "Add filtering for number of motifs!\n";
	open R, ">", $output_R_file;
	print R "pdf(\"" . substr($output_R_file, 0, length($output_R_file) - 2) . ".pdf\", width=10, height=5)\n";
	my $cal_pvalue = 0;
	my $count_obs_one = 0;
	my $count_obs_two = 0;
	my $count_obs_both = 0;
	foreach my $motif (sort {$index_motifs{$a} cmp $index_motifs{$b}} keys %index_motifs) {
		#Make sure there a mutations in this motif
		if($removed == 0) {
			$filename = $output . "_" . $motif . "_removed.txt";
		} else {
			$filename = $output . "_" . $motif . ".txt";
		}
		$num_of_muts = `wc -l $filename`;
		@split = split('\s+', $num_of_muts);
		if($split[0] == 1 && $removed == 1) {
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

sub fakrek{
        my $number = shift;
        return undef if $number < 0;
        return 1 if $number == 0;
        return ($number * &fakrek ($number -1));
}
