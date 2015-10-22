#n!/usr/bin/perl 
#author: Verena Link

#Software to analyze motif differences of genes that are differently expressed between different mouse strains

use strict;
use Getopt::Long;
require '../general/config.pm';
use Data::Dumper;
require '../db_part/database_interaction.pm';

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
	print STDERR "\t-method: rna or genomic: Analyzes TSS or genomic loci\n";
	print STDERR "\t-strains <strains>: Comma-separated list of strains - if not defined only the reference is used - if all is specified all strains are used\n";
	print STDERR "\t-file <file>: Genes for which the promoter should be annotated\n";
	print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-d <tag dirs>: List with tag directories (RNA-Seq)\n";
	print STDERR "\t-ann <file>: File with TSS annotation!\n";
	print STDERR "\t-no-first: Skips annotation of first expressed exon\n";
	print STDERR "\t-no-diff: Skips annotation of first differently expressed exon\n";
	print STDERR "\t-fc: foldchange that defines differently expressed exons (default: 4.0)\n";
	print STDERR "\t-tss <down>,<up>: Extends promoter from TSS or first expressed base (default: -500,500)\n";
	print STDERR "\t-quant: Instead of looking for difference in motif occurrence it checks difference in motif scores (more sensitive method)\n";
	print STDERR "\t-motif-diff: Minimal difference between motif scores to consider (default: there is a differences in general)\n";
	print STDERR "TODO percentage of perfect motif score!\n";
	print STDERR "\t\tExample: -tss -500,500: Extends TSS 500 in both directions\n";
	print STDERR "\t\tExample: -tss 500: Extends TSS 500 in both directions\n";
	print STDERR "\nFilter methods:\n";
	print STDERR "\t-filter <count>: Filters out genes with less than <count> tag counts in all exons together\n";
	print STDERR "\nOutput:\n";
	print STDERR "\t-html: Outputs a html file with the sequence and the motifs\n";
	print STDERR "\t-bed: Outputs a bedgraph with candidates motif per strain\n";
	print STDERR "\n\nBackground (Default: all TSS):\n";
	print STDERR "\t-bgfile: <File with RefSeq Identifiers to use as background\n";
	print STDERR "\t-bgtagdir: <reference and strain RNA-Seq tag dir to calculate equally expressed genes\n";
	print STDERR "\n\nRun methods (Default: find pattern and calculate statistics for it)\n";
	print STDERR "\t-no-stats (Does not calculate statistics for results)\n";
	print STDERR "\t-no-pattern (Des not try to find a pattern in motif mutations that overlaps with gene expression pattern\n";
	print STDERR "\n\nMotif operations:\n";
	print STDERR "\t-motif-file <file>: Path to file with PWM for motifs (default: homers all.motif)\n";
	print STDERR "\t-zscore <number>: The z-score used to define an event as significant (default: 1.5)\n";
	exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $noise = 0.00000001;

#define variables
my %diff_motifs = ();
my %motif_score = ();
my %motif_pos = ();
my $promoter_start = 0;
my $promoter_end = 0;
my $current_chr;
my %save_motifs = ();
my %save_motifs_aligned = ();
my %save_alignment = ();
my $f = 0;
my %pos;
my $check = 0;
my $done_tasks = 0;
my %target = ();
my @b;
my $done_tss = 0;
my $m;
my $zscore = 0;
my $cand_tricky = "";
my $pattern_up = "";
my $pattern_down = "";
my %strain_is_different_up = ();
my %strain_is_different_down = ();


my $tmp_bgfile;
my %background = ();
my $output;
my @detail_motifs;
my $m;
my $count = 0;
my $number_of_tss = 0;
my $sum = 0;
my %mean;
my %stddev;
my $sum;
my $sum_square;
my $number;
my %perfect_motif_score = ();

$_ = "" for my ($genome, $table_name, $ann_file, $name, $file, $seqs, $command, $bgfile, $bgtagdir, $motif_file, $method);
$_ = 0 for my ($homo, $num, $skip_first, $skip_diff, $tss_up, $tss_down, $filter, $html, $html_output, $tasks, $stats, $pattern, $occurrence, $zscore_significant_value, $quant, $get_detailed_information, $bed_output, $motif_diff, $pot);
my ($dbh, $sth, @strains, @dirs, %tss, %first, %diff, @split, %to_del, @part, %genes, @tss, %filter, %del, %file_genes, %seqs, $cand_up, $cand_down, @candidate_strains_up, @candidate_strains_down, @candidate_strains_equal, @fileHandles, %strains_down_for_stats, %strains_up_for_stats);
my $homo = 1;
my %cand_motif = ();
my %cands;
my %save_pos = ();
my %seq_motif = ();
my %orientation = ();


my %rev_comp = ();
$rev_comp{'A'} = 'T';
$rev_comp{'T'} = 'A';
$rev_comp{'G'} = 'C';
$rev_comp{'C'} = 'G';


my %motif_matrix = ();

my %mandatory = ('-genome' => 1, '-strains' => 1, "-method" => '1');
my %convert = map { $_ => 1 } @ARGV;
my $fc = 4;

my $param = config::read_config();

config::check_parameters(\%mandatory, \%convert);

GetOptions(     "genome=s" => \$genome,
		"method=s" => \$method,
		"file=s" => \$file,
                "strains=s{,}" => \@strains,
                "homo" => \$homo,
                "d=s{,}" => \@dirs, 
		"ann=s{,}" => \$ann_file, 
		"-no-first" => \$skip_first,
		"-no-diff" => \$skip_diff,
		"-fc=s" => \$fc, 
		"-tss=s{,}" => \@tss, 
		"-filter=s" => \$filter, 
		"-html" => \$html_output,
		"-bed" => \$bed_output, 
		"-bgfile=s" => \$bgfile,
		"-bgtagdir" => \$bgtagdir,
		"-no-stats" => \$stats, 
		"-no-pattern" => \$pattern, 
		"-motif_file=s" => \$motif_file, 
		"-zscore=s" => \$zscore_significant_value, 
		"-quant" => \$quant, 
		"-motif-diff" => \$motif_diff)
	or die("Error in command line options!\n");

my $ref_in_array = 0;
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	if($strains[$i] eq "reference") {
		$ref_in_array = 1;
	}
}

if($ref_in_array == 0) {
	unshift @strains, "reference";
}


if($html_output == 1 || $quant == 1 || $bed_output == 1) {
	$get_detailed_information = 1;
}

if($motif_file eq "") {
	$motif_file = "all.motifs";
	print STDERR "PATH IS NOT CORRECT - THAT NEEDS TO BE FIXED!\n";
}
database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");
if($zscore_significant_value == 0) {
	$zscore_significant_value = 1.5;
}
@tss = split(',', $tss[0]);
if(@tss == 0) {
	$tss_down = 500;
	$tss_up = 500;
} elsif(@tss == 1) {
	$tss_down = $tss[0];
	$tss_up = $tss[0];
} else {
	if($tss[0] < 0) {
		$tss_down = ($tss[0]) * -1;
	} else {
		$tss_down = $tss[0];
	}
	$tss_up = $tss[1];
}


my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}", $con->{'user'},$con->{'pw'}, {'RaiseError' => 1});

if($method eq "rna") {
	#Report TSS
	$table_name = $genome . "_annotations";
	if($ann_file eq "" && &check_if_table_exists($table_name) == 1) {
		&tss_from_db()
	} elsif($ann_file ne "") {
		&tss_from_file();
	} else {
		&tss_from_homer();
	}
	print STDERR "Officially annotated TSS is stored!\n";
	print STDERR "Go to first expressed exon!\n";

	#Report first expressed exon
	if(@dirs == 0) {
		print STDERR "No tag directories are specified!\n";
		print STDERR "Skip first expressed exon and differently expressed exon!\n";
		$skip_first = 1;
		$skip_diff = 1;	
	} else {
		if(@dirs != @strains) {
			print STDERR "Number of tag directories and strains are not the same!\n";
			print STDERR "Can not define what tag directory belongs to which strains!\n";
			print STDERR "Skip first expressed exon and differently expressed exon annotation!\n";
			$skip_first = 1;
			$skip_diff = 1;
		}
	}

} else {
	print STDERR "Add enhancer routine!\n";
	$skip_first = 1;
	$skip_diff = 1;
	exit;
}


#Open filehanldes for output bed file
my $filename;
if($bed_output == 1) {
	for(my $i = 0; $i < @strains; $i++) {
		$filename = $strains[$i];
		open my $fh, ">", "diff_motifs_$filename.bed" or die "Can't open $filename: $!";
		push @fileHandles, $fh;
		$fileHandles[$i]->print("track name=\"candidate-mutations-$filename\"\n");
	}
}

if($skip_first == 0) {
	&first_expressed_exon();
}
#Report first differently expressed exon
if($skip_diff == 0) {
	&diff_expressed_exon();
}

#Filter out all genes that have less than $filter tag counts
if($filter > 0) {
	print STDERR "Collect tag counts per gene\n";
	print STDERR "TODO\n";
	exit;
	my $command = "analyzeRepeats.pl rna " . $genome . " -count exons -d";
	for(my $i = 0; $i < @dirs; $i++) {
		$command .= " " . $dirs[$i];
	}
	$command .= " > gene_counts.tmp";
	$to_del{'gepne_counts.tmp'} = 1;

	`$command`;

	open FH, "<gene_counts.tmp";
	my $f = 0;
	my $del = 0;
	foreach my $line (<FH>) {
		chomp $line;
		if($f == 0) { $f++; next; }
		$del = 0;
		for(my $i = 8; $i < @split; $i++) {
			if($split[$i] > $filter) { $del = 1; } 
		}	
		if($del == 0) { $filter{$split[0]} = 1; }
	}		
	close FH;
}

#Read in file with genes that should be considered - do not run over all genes
%file_genes = ();
if($file ne "" && $method eq "rna") {
	print STDERR "Read in file with genes\n";
	open FH, "<$file";
	foreach my $line (<FH>) {
		@split = split('\t', $line);
		$file_genes{$split[0]} = $line;
		$tasks++;
	}
	close FH;
} 
else {
	print STDERR "Again check here!\n";
	exit;
	%file_genes = %tss;
	$tasks = keys %tss;
}


#Calculate the background
#Check if alternative background is used
my $background_path = $param->{'homer_path'} ;

#Calculate the z-score of every motif that is mutated between two strains
if($stats == 0) { 
	print STDERR "We also need to check that with the methods!\n";
	&read_in_background();
	print STDERR "Calculate mean and stddev\n";
	&cal_mean_and_stddev();
}


open PROMOTER, ">promoter_analysis.txt";
&add_header();
if($quant == 1) {
	&read_in_motif_files();
}
my $chr;
my $tasks = keys %file_genes;
my $done_tasks = 0;
foreach my $line (keys %tss) {
	if(exists $filter{$line}) { next; }
	@split = split('\t', $tss{$line}->{'tss'});
	#Gene in file - get motif differences
	if(exists $file_genes{$line}) {
		$done_tss++;
		chomp $file_genes{$line};
		print PROMOTER $file_genes{$line};
		$chr = (split('\t', $file_genes{$line}))[1];
		&scan_seq_for_motifs($split[0], ($split[1] - $tss_down), ($split[1] + $tss_up));
		&save_motifs();
		&align_motifs();
		if($quant == 0) {
			%diff_motifs = ();
			&count_motifs();
		} else {
			%motif_score = ();
			&calculate_motif_differences();
		}
		#Now find pattern
		$cand_up = "";
		$cand_down = "";
		$cand_tricky = "";
		&analyze_patterns($file_genes{$line});
		if($pot == 1) {
			print PROMOTER "UP: " . $cand_up . "\tDOWN: " . $cand_down . "\tTRICKY: " . $cand_tricky;
		} else {
			print PROMOTER "no candidate";
		}
		if($stats == 0) {
			if($quant == 1) {

			#	print STDERR "TODO STATS FOR MOTIF DIFFERENCES\n";
			} else {
				&calculate_stats($cand_up, \@candidate_strains_up);
				&calculate_stats($cand_down, \@candidate_strains_down);
			}
		}

		if($bed_output == 1) {
			&print_bedGraph_output(($split[1] - $tss_down));
		}

		if($html_output == 1) {
			&print_html_output("tss_" . $line);
		}

		#Add logic for the other promoters!
	}
	if(exists $file_genes{$line} || exists $first{$line} || exists $diff{$line}) {
		print PROMOTER "\n";
		$done_tasks++;
	}
	print STDERR "Status\tdone: $done_tasks\t vs \ttaks: $tasks Completed\r";

}
print STDERR "TSS is done!\n";
my %seen = ();
foreach my $line (keys %first) {
	print "show: " . $line . "\n";
	if(exists $file_genes{$line}) {
		%seen = ();
		foreach my $strain (keys %{$first{$line}}) {
			print $first{$line}{$strain}->{'exon_start'} . "\t" . $first{$line}{$strain}->{'exon_chr'} . "\t" . $first{$line}{$strain}->{'exon_strand'} . "\n";;
			if(exists $seen{$first{$line}{$strain}->{'exon_start'}}) {
				next;
			}
			$seen{$first{$line}{$strain}->{'exon_start'}} = 1;
			if($first{$line}{$strain}->{'exon_strand'} eq "+") {
				&scan_seq_for_motifs($first{$line}{$strain}->{'exon_chr'}, $first{$line}{$strain}->{'exon_start'} - (2 * $tss_up), $first{$line}{$strain}->{'exon_start'});
			} else {
				&scan_seq_for_motifs($first{$line}{$strain}->{'exon_chr'}, $first{$line}{$strain}->{'exon_start'}, $first{$line}{$strain}->{'exon_start'} + (2 * $tss_up));
			}
			&save_motifs();
			print Dumper %save_motifs . "\n";
			exit;
#			&align_motifs();
#			if($quant == 0) {
#				%diff_motifs = ();
#				&count_motifs();
#			} else {
#				%motif_score = ();
#				&calculate_motif_differences();
#			}
#			#Now find pattern
#			$cand_up = "";
#			$cand_down = "";
#			$cand_tricky = "";
#			&analyze_patterns($file_genes{$line});
#			if($pot == 1) {
#				print PROMOTER "UP: " . $cand_up . "\tDOWN: " . $cand_down . "\tTRICKY: " . $cand_tricky;
#			} else {
#				print PROMOTER "no candidate";
#			}
#			if($stats == 0) {
#				if($quant == 1) {
#
#				#	print STDERR "TODO STATS FOR MOTIF DIFFERENCES\n";
#				} else {
#					&calculate_stats($cand_up, \@candidate_strains_up);
#					&calculate_stats($cand_down, \@candidate_strains_down);
#				}
#			}
#
#			if($bed_output == 1) {
#				&print_bedGraph_output(($split[1] - $tss_down));
#			}
#
#			if($html_output == 1) {
#				&print_html_output("tss_" . $line);
#			}
#
#			#Add logic for the other promoters!
#		}
#		if(exists $file_genes{$line} || exists $first{$line} || exists $diff{$line}) {
#			print PROMOTER "\n";
#			$done_tasks++;
#		}
#		print STDERR "Status\tdone: $done_tasks\t vs \ttaks: $tasks Completed\r";
		}	
	}
}
exit;
print "\n\n";
print STDERR "Now add rest of TSS and enhancers!!!!\n";

if($bed_output == 1) {
	for(my $i = 0; $i < @strains; $i++) {
		close $fileHandles[$strains[$i]];
	}

}

#Connects to the database and saves the annotated TSS for this genome
sub tss_from_db{
	$sth = $dbh->prepare("SELECT * FROM "  . $genome . "_annotations where annotation = \'P\'");
	$sth->execute();
	while (my $gene = $sth->fetchrow_hashref()) { 
	#	$tss{$gene->{'gene'}} = $gene->{'chr'} . "\t" . $gene->{'start'} . "\t" . $gene->{'stop'};
		$tss{substr($gene->{'gene'}, 0, length($gene->{'gene'}) - 1)}->{'tss'} = $gene->{'chr'} . "\t" . $gene->{'start'} . "\t" . $gene->{'stop'};
		$tss{substr($gene->{'gene'}, 0, length($gene->{'gene'}) - 1)}->{'strand'} = $gene->{'strand'};
	}
}

#Reads in the genomic coordinates for the transcription start site from a file defined by the user
#The format of this file is determined in the config
sub tss_from_file{
	my ($split, $ann, $nf) = config::promoter_format();	
	open FH, "<$ann_file";

	my @str_name = split('\t', $split);
	my $final_name;
	print STDERR "processing $ann_file\n";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if($split[$nf->{'annotation'}] eq $ann) {
			if(@str_name == 1) {
				$final_name = substr($split[$nf->{'name'}], $str_name[0]);			
			} else {
				$final_name = substr($split[$nf->{'name'}], $str_name[0], length($split[$nf->{'name'}]) - $str_name[1] - $str_name[0]);
			}
			if(substr($split[$nf->{'chr'}], 0, 3) eq "chr") {
				$split[$nf->{'chr'}] = substr($split[$nf->{'chr'}], 3);
			}
			$tss{$final_name}->{'tss'} = $split[$nf->{'start'}] . "\t" . $split[$nf->{'stop'}] . "\t" . $split[$nf->{'chr'}];
			$tss{$final_name}->{'strand'} = $split[$nf->{'strand'}];
		}
	}
}

#checks if the used genome has a genomic folder in homer and if this folder contains a annotation file
#If available it uses this annotation file to find all TSS for this genome
#Format has to follow the standard homer annotation file format
sub tss_from_homer{
	print STDERR "Annotation for " . $genome . " is not saved in the database!\n";
	print STDERR "\nCheck homer for genome annotation!\n";
	
	$ann_file = $param->{'homer_path'} . "/data/genomes/" . $genome . "/"  . $genome . ".full.annotation";
	if(!(-e $ann_file)) {
		print STDERR "No genome annotation found in homer!\n";
		next;
	}
	print STDERR "\nGenome annotation found!\n";
	open FH, "<$ann_file";
	print "processing $ann_file\n";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if($split[5] eq "P") {
			$tss{substr($split[0], 14, length($split[0]) - 15)}->{'tss'} = substr($split[1], 3) . "\t" . $split[2] . "\t" . $split[3];			
			$tss{substr($split[0], 14, length($split[0]) - 15)}->{'strand'} = $split[4];
		}
	}	
}

#Needs tag directories to function
#Runs analyzeRepeats (a homer program) to get the tag counts per exon (-noCondensingParts)
#Saves the first expressed exon in a hash called first
sub first_expressed_exon{
	$command = "analyzeRepeats.pl rna " . $genome . " -count exons -d";
	for(my $i = 0; $i < @dirs; $i++) {
		$command .= " " . $dirs[$i];
	}
	$command .= " -noCondensingParts > genes.tmp";
	$to_del{'genes.tmp'} = 1;
	print STDERR "COMMENT COMMAND IN AGAIN!!!\n";
#	`$command`; 

	open FH, "<genes.tmp";
	my $f = 0;

	foreach my $line (<FH>) {
		chomp $line;
		if($f == 0) {
			$f++;
			next;
		}	
		@split  = split('\t', $line);
		@part = split('--', $split[0]);	
		for(my $i = 8; $i < @split; $i++) {
			if($split[$i] > 0) {
				#No entry for this gene yet or gene on positive strand and this part is before the saved part or gene on negative strand and this part is after the saved part
				if(!exists $first{$part[0]}->{$strains[$i-8]} || ($first{$part[0]}->{$strains[$i-8]}->{'part'} > substr($part[1], 4) && $split[4] eq "+") || ($first{$part[0]}->{$strains[$i-8]}->{'part'} < substr($part[1], 4) && $split[4] eq "-")) {
					$first{$part[0]}->{$strains[$i-8]}->{'part'} = substr($part[1], 4);
					$first{$part[0]}->{$strains[$i-8]}->{'exon_start'} = $split[2];
					$first{$part[0]}->{$strains[$i-8]}->{'exon_chr'} = substr($split[1], 3);
					$first{$part[0]}->{$strains[$i-8]}->{'exon_strand'} = $split[4];
				}
			}
		}	
	}
	print STDERR "First expressed exon is saved!\n";
	close FH;
}

#Needs tag directories to function
#Runs analyzeRepeats (a homer program) to get the tag counts per exon (-noCondensingParts)
#Saves the first exon that has a greater foldchange between to strains than defined in $fc
sub diff_expressed_exon{
	if(!exists $to_del{'genes.tmp'}) {
		$command = "analyzeRepeats.pl rna " . $genome . " -count exons -d";
		for(my $i = 0; $i < @dirs; $i++) {
			$command .= " " . $dirs[$i];
		}
		$command .= " -noCondensingParts > genes.tmp";
		$to_del{'genes.tmp'} = 1;

		`$command`; 
	}
	print STDERR "Get first differently expressed exon!\n";
	open FH, "<genes.tmp";
	my $f = 0;

	foreach my $line (<FH>) {
		if($f == 0) {
			$f++;
			next;
		}
		chomp $line;
		@split = split('\t', $line);
		@part = split('--', $split[0]);
		for(my $i = 8; $i < @split - 1; $i++) {
			for(my $j = $i + 1; $j < @split; $j++) {
				if(($split[$i] + 1)/($split[$j] + 1) > $fc || ($split[$j] + 1)/($split[$i] + 1) > $fc) {
					if(!exists $diff{$part[0]}->{$strains[$i-8]}->{$strains[$j-8]} || $diff{$part[0]}->{$strains[$i-8]}->{$strains[$j-8]}->{'part'} > substr($part[1], 4)) {
						$diff{$part[0]}->{$strains[$i-8]}->{$strains[$j-8]}->{'part'} = substr($part[1], 4);
						$diff{$part[0]}->{$strains[$i-8]}->{$strains[$j-8]}->{'exon_start'} = $split[2];
						$diff{$part[0]}->{$strains[$i-8]}->{$strains[$j-8]}->{'exon_chr'} = substr($split[1], 3);
						$diff{$part[0]}->{$strains[$i-8]}->{$strains[$j-8]}->{'exon_strand'} = $split[4];
					}
				}
			}
		}
	}
	print STDERR "Differently expressed exon is saved!\n";
	close FH;
}
my @k;
my %ref;
my %gapless;
my $tmp_seq;
my $length;
my $first = 0;
my $print;
my $equal = 0;
my $val;
my @pos;
my @chr;
my $middle;
my $start;
my @array;
my @name;
my $num = 0;
my %output = ();
my %save = ();
my @split;
my @second;

sub calculate_motif_differences {
	my $ref_score;
	my $score;
	my $motif_seq;	
	my $run;
	my $l;
	my @save_tmp;
	my $save_in_array = 0;
	my $current_orientation = "";
	foreach my $motif (keys %motif_pos) {
		foreach my $pos (keys %{$motif_pos{$motif}}) {
			for(my $i = 0; $i < @strains; $i++) {
				if($current_orientation eq "") {
					$current_orientation = $save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'orientation'};
				}
				if(!exists $save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'score'}) {
					$save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'score'} = "not found";
				}

			}
			$ref_score = 0;
			if(!exists $perfect_motif_score{$motif}) {
				$perfect_motif_score{$motif} = &calculate_perfect_motif_score($motif, $motif_pos{$motif}{$pos});
			}
			for(my $i = 0; $i < @strains; $i++) {
				$run = $pos;
				$l = $pos;
				$motif_seq = "";
				#Generate motif sequence from aligned sequence without gaps
				while($l < $pos + $motif_pos{$motif}{$pos}) {
					if(substr($seqs{$strains[$i]}, $run, 1) ne "-") {
						$motif_seq .= substr($seqs{$strains[$i]}, $run, 1);
						$l++;
					}
					$run++;
				}
				if(!exists $save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'orientation'}) {
					$save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'orientation'} = $current_orientation;
				}
				$score = &calculate_motif_score($motif, $motif_seq, $save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'orientation'});
				if($ref_score == 0) { $ref_score = $score; }
				$save_tmp[$i] = $score;
				if($ref_score != $score) { 
					$save_in_array = 1;
				}
			}
			if($save_in_array == 1) {
				for(my $i = 0; $i < @save_tmp; $i++) {
					$motif_score{$motif}{$pos}{$strains[$i]} = $save_tmp[$i];
				}
			}
			$save_in_array = 0;
			$current_orientation = "";
		}
	}
}

sub count_motifs {
	my $ref_count;
	foreach my $motif (keys %save_motifs_aligned) {
		$ref_count = keys %{$save_motifs_aligned{$motif}{'reference'}};
		foreach my $strain (keys %{$save_motifs_aligned{$motif}}) {
			if($strain eq "reference") { next; }
			if($ref_count != (keys %{$save_motifs_aligned{$motif}{$strain}})) {
				$diff_motifs{$motif} = 1;
			}
		}
	}
}


sub align_motifs {
	my $half;
	my $motif_start;
	my $motif_seq_from_file;
	my $aligned_seq;
	my @seq_with_gaps;
	my @motif_array;
	my $l;
	my $run;
	%save_motifs_aligned = ();
	if($html_output == 1) {
		%save_alignment = ();
	}
	if($quant == 1) {
		%motif_pos = ();
	}
	my $strain;
	foreach my $motif (keys %save_motifs) {
		for(my $i = 0; $i < @strains; $i++) {
			$strain = $strains[$i];
			$half = int(length($gapless{$strain})/2);
			$aligned_seq = "";
			$start = 0;
			$run = 0;
			$l = 0;
			@seq_with_gaps = split('', $seqs{$strain});
			foreach my $pos (sort {$a <=> $b} keys %{$save_motifs{$motif}{$strain}}) {
				$motif_start = $half + $pos;
				$motif_seq_from_file = $save_motifs{$motif}{$strain}{$pos}{'seq'};
				if($save_motifs{$motif}{$strain}{$pos}{'orientation'} eq "-") {
					$motif_start = $motif_start - length($motif_seq_from_file) + 1;
				}
				#Motif also appears on the negative strand - do we really want to skip it?
				if($motif_start < $l) { next; }
				@motif_array = split('', substr($gapless{$strain}, $motif_start, length($motif_seq_from_file)));
				while($l < $motif_start) {
					if($seq_with_gaps[$run] eq "-") {
						$aligned_seq .= "-";
					} else {
						$aligned_seq .= " ";
						$l++;
					}
					$run++;
				}
				$save_motifs_aligned{$motif}{$strain}{$run}{'orientation'} = $save_motifs{$motif}{$strain}{$pos}{'orientation'};
				$save_motifs_aligned{$motif}{$strain}{$run}{'seq'} = $save_motifs{$motif}{$strain}{$pos}{'seq'};
				$save_motifs_aligned{$motif}{$strain}{$run}{'score'} = $save_motifs{$motif}{$strain}{$pos}{'score'};
				while($l < $motif_start + length($motif_seq_from_file)) {
					if($seq_with_gaps[$run] eq "-") {
						$aligned_seq .= "-";
					} else {
						$aligned_seq .= $seq_with_gaps[$run];
						$l++;
					}
					$run++;
				}
				if($quant == 1) {
					$motif_pos{$motif}{$run - length($motif_seq_from_file)} = length($motif_seq_from_file);
				}
			}
			if($html_output == 1) {
				for(my $i = length($aligned_seq); $i < length($seqs{'reference'}); $i++) {
					if($seq_with_gaps[$i] eq "-") {
						$aligned_seq .= "-";
					} else {
						$aligned_seq .= " ";
					}
				}
				$save_alignment{$motif}{$strain} = $aligned_seq;
			}
		}
	}
}

sub extract_seqs_bg {
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");
	my $ref = database_interaction::get_multiple_alignment($_[1], $_[2], $_[0], "+", 0, 0, 0, 0);
	for(my $i = 0; $i < @{$ref}; $i++) {
		print BG_SEQ ">" . $strains[$i] . "_" . $_[0] . "_" . $_[1] . "_" . $_[2] . "\n";
		$tmp_seq = $ref->[$i];
		$tmp_seq =~ s/-//g;
		print BG_SEQ $tmp_seq . "\n";
	}
}

#Method to count the differences between motifs for each strain 
sub calculate_bg {
	my %sum_up_motifs;
	open FH, "<tmp_bgseqs_with_motifs.txt";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@second = split('_', $split[0]);
		if(!exists $sum_up_motifs{$second[2]}{$split[3]}{$second[0]}) {
			$sum_up_motifs{$second[2]}{$split[3]}{$second[0]} = 0;
		}
		$sum_up_motifs{$second[2]}{$split[3]}{$second[0]}++;	
	}
	close FH;
	foreach my $pos (keys %sum_up_motifs) {
		foreach my $motif (keys %{$sum_up_motifs{$pos}}) {
			foreach my $strain (keys %{$sum_up_motifs{$pos}{$motif}}) {
				if($strain eq "reference") { next; }
				$background{$strain}{$motif}{$sum_up_motifs{$pos}{$motif}{'reference'}-$sum_up_motifs{$pos}{$motif}{$strain}}++;
			}
		}
	}
}

sub scan_seq_for_motifs {
	#Extract sequences from database and also get MSA
	open OUT, ">seqs.txt";
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");
	my $ref = database_interaction::get_multiple_alignment($_[1], $_[2], $_[0], "+", 0, 0, 0, 0);
	$length = 0;
	#Get the longest sequences from MSA
	#Save the MSA and also the sequences without any gaps
	for(my $i = 0; $i < @{$ref}; $i++) {
		$seqs{$strains[$i]} = $ref->[$i];
		$tmp_seq = $ref->[$i];
		$tmp_seq =~ s/-//g;
		print OUT ">" . $strains[$i] . "\n";
		print OUT $tmp_seq . "\n";
#		if($html_output == 1) {
		$gapless{$strains[$i]} = $tmp_seq;
#		}
	}
	close OUT;

	#Step 2: Scan all gapless sequences for all motifs in the motif file
	my $command = "homer2 find -i seqs.txt -m all.motifs > seqs_with_motifs.txt 2> output_homer.tmp";
	`$command`;
}

sub save_motifs {
	open FH, "<seqs_with_motifs.txt";
	%save_motifs = ();
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		#save motifs{motif}{strain}{position}
		#also save in this hash: orientation, sequence and score (maybe useful?)
		$save_motifs{$split[3]}{$split[0]}{$split[1]}{'orientation'} = $split[4];
		$save_motifs{$split[3]}{$split[0]}{$split[1]}{'score'} = $split[5];
		$save_motifs{$split[3]}{$split[0]}{$split[1]}{'seq'} = $split[2];
	}
}

my %mot_pos = ();
#Method to extract sequences from database and scan them for motif 
#It also saves the MSA for all strains, so we can align the motifs later
#Saves motifs for every strain and every position and compares the number of motifs for each strain - saves motifs with different occurrences in a hash calld %cands

sub calculate_motif_score{
	my $motif_score = 0;
	my $local_motif = $_[0];
	my $strain_seq = $_[1];
	my $orientation = $_[2];
	my $greater = 0;
	#That needs to move someqhere else
	if($orientation eq "-") {
		$strain_seq = &rev_comp($strain_seq);
	}
	@split = split('', $strain_seq);
	$motif_score = 0;
	for(my $i = 0; $i < @split; $i++) {
		$motif_score += log(($motif_matrix{$local_motif}{length($strain_seq)}{$i}{$split[$i]}+$noise)/0.25);
	}
	$motif_score = sprintf("%.6f", $motif_score);
	return $motif_score;
}

sub calculate_perfect_motif_score {
	my $perfect_score = 0;
	my $local_motif = $_[0];
	my $max = 0;
	my $length = $_[1];
#	%{$motif_matrix{$save_header}{$save_length}} = %matrix;
#	print $local_motif . "\n";
	foreach my $pos (sort {$a <=> $b } keys %{$motif_matrix{$local_motif}{$length}}) {
		$max = 0;
		foreach my $base (keys %{$motif_matrix{$local_motif}{$length}{$pos}}) {
	#		print $pos . "\t" . $base . "\t" . $motif_matrix{$local_motif}{$length}{$pos}{$base} . "\n";
			if($max < log(($motif_matrix{$local_motif}{$length}{$pos}{$base}+$noise)/0.25)) {
				$max = log(($motif_matrix{$local_motif}{$length}{$pos}{$base}+$noise)/0.25);
			}
		}
	#	print $max . "\n";
		$perfect_score += $max;
	#	print "perfect: " . $perfect_score . "\n";
	}
        $perfect_score = sprintf("%.6f", $perfect_score);
#	print $perfect_score . "\n";
	return $perfect_score;
}
#Function to create the reverse complementary motif seqeunce
sub rev_comp{
	my $seq = $_[0];
	my $rev_seq = reverse($seq);
	@split = split('', $rev_seq);
	my $rev_comp_seq = "";
	for(my $i = 0; $i < @split; $i++) {
		$rev_comp_seq .= $rev_comp{$split[$i]};
	}
	return $rev_comp_seq;
}


#Reads in the motif files and saves it in a hash matrix so it is possible to later calculate the motif scores
sub read_in_motif_files{
	my %matrix = ();
	my $motif_counter = 0;
	my $save_header = "";
	my $save_length = 0;
	$motif_counter = 0;
	open FH, "<$motif_file";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if(substr($line, 0, 1) eq ">") {
			if(keys %matrix == 0) { 
				next; 
			} else {
				%{$motif_matrix{$save_header}{$save_length}} = %matrix;
				%matrix = ();
				$motif_counter = 0;	
			}
			$save_header = $split[1];
			$save_length = length($split[0]) - 1;
		} else {
			$matrix{$motif_counter}{'A'} = $split[0];	
			$matrix{$motif_counter}{'C'} = $split[1];	
			$matrix{$motif_counter}{'G'} = $split[2];	
			$matrix{$motif_counter}{'T'} = $split[3];	
			$motif_counter++;
		}
	}
	
	%{$motif_matrix{$save_header}{$save_length}} = %matrix;
}

#Function to output a html file with the alignment of the sequences for the different strains and the motif alignments
sub print_html_output{
	my $name = $_[0];
	open OUT, ">seq_$name.html";
	my %ref = %seqs;
	my $ref_strain = $ref{'reference'};
	
	while(length($ref{$strains[0]}) > 0) {
		for(my $i = 0; $i < @strains; $i++) {
			while(length($ref{$strains[$i]}) > 0 && length($ref{$strains[$i]}) < 100) {
				$ref{$strains[$i]} .= " ";
			}
			print OUT "<pre>";
			for(my $r = 0; $r < 10; $r++) {
				print OUT substr($ref{$strains[$i]}, $r * 10, 10) . "|";
			}
			print OUT "\t" . $strains[$i] . "\n";
			$ref{$strains[$i]} = substr($ref{$strains[$i]}, 100);
		}
		$ref_strain = substr($ref_strain, 100);
		foreach my $motif (keys %save_alignment) {
			for(my $i = 0; $i < @strains; $i++) {
				while(length($save_alignment{$motif}{$strains[$i]}) > 0 && length($save_alignment{$motif}{$strains[$i]}) < 100) {
					$save_alignment{$motif}{$strains[$i]} .= " ";
				}
				print OUT "<pre>";
				for(my $r = 0; $r < 10; $r++) {
					print OUT substr($save_alignment{$motif}{$strains[$i]}, $r * 10, 10) . "|";
				}
				$save_alignment{$motif}{$strains[$i]} = substr($save_alignment{$motif}{$strains[$i]}, 100);
				print OUT "\t" . $motif . " - " . $strains[$i] . "\n";
			}
		}
	}
	close OUT;
}

sub print_bedGraph_output{
	my $tss_start = $_[0];
	my %combined = (%strains_down_for_stats, %strains_up_for_stats);
	my $strain;
	my $length = 0;
	my $gaps = 0;
	foreach my $motif (keys %combined) {
		my @short_name = split('/', $motif);
		if($quant == 0) {
			for(my $i = 0; $i < @strains; $i++) {
				$strain = $strains[$i];
				foreach my $pos (keys %{$save_motifs_aligned{$motif}{$strain}}) {
					$gaps = () = substr($seqs{'reference'}, 0, $pos) =~ /-/g;
					$fileHandles[$i]->print($chr . "\t" . ($tss_start + $pos - $gaps) . "\t" . ($tss_start + $pos - $gaps + length($save_motifs_aligned{$motif}{$strain}{$pos}{'seq'})) . "\t" . $short_name[0] . "($save_motifs_aligned{$motif}{$strain}{$pos}{'seq'}) \t1\t+\n");
				}
			}
		} else {
			foreach my $pos (keys %{$combined{$motif}}) {
				$gaps = () = substr($seqs{'reference'}, 0, $pos) =~ /-/g;
				for(my $i = 0; $i < @strains; $i++) {
					if(exists $save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'seq'}) { $length = length($save_motifs_aligned{$motif}{$strains[$i]}{$pos}{'seq'}); }
				}
				for(my $i = 0; $i < @strains; $i++) {
					$fileHandles[$i]->print($chr . "\t" . ($tss_start + $pos - $gaps) . "\t" . ($tss_start + $pos - $gaps + $length) . "\t" . $short_name[0] . ":" . $motif_score{$motif}{$pos}{$strains[$i]} . "\t1\t+\n");
				}
			}
		}
	}	
}




#Add header to promoter analysis output file
sub add_header{
	print PROMOTER "ID\tchr\tstart\tstop\tstrand\tlength\tcopy\tannotation";
	for(my $i = 0; $i < @strains; $i++) {
		print PROMOTER "\t$strains[$i]";
	}
	for(my $i = 1; $i < @strains; $i++) {
		print PROMOTER "\ttss mut " . $strains[$i];
	}
	for(my $i = 0; $i < @strains; $i++) {
		print PROMOTER "\tlost in $strains[$i]";
	}
	print PROMOTER "\n";
}

#Method to calculate mean and standard deviation for the motifs in the background
sub cal_mean_and_stddev {
	foreach my $strain (keys %background) {
		foreach my $motif (keys %{$background{$strain}}) {
			$sum = 0;
			$number = 0;
			$sum_square = 0;
			foreach my $pos (keys %{$background{$strain}{$motif}}) {
				$sum += ($pos * $background{$strain}{$motif}{$pos});
				$sum_square += ($pos*$pos) * $background{$strain}{$motif}{$pos};
				$number += $background{$strain}{$motif}{$pos};
			}
			$mean{$strain}{$motif} = ($sum/$number);
			$stddev{$strain}{$motif} = sqrt(($sum_square/$number) - ($sum/$number)**2);
		}
	}
}

#Method to read in the background distribution of the motifs or create a new background
my $all_bg_exists = 0;
sub read_in_background {
	#Background file exists 
	#File with genes used for background exists - read in genes and calculate the bg
	if($bgfile ne "") {
		print STDERR "Reading in file with genes for background!\n";
		my $number_of_lines = &get_lines_in_file($bgfile); 
		open FH, "<$bgfile";
		foreach my $line (<FH>) {
			chomp $line;
			$count++;
			print STDERR "Reading in file: " . $count . " of " . $number_of_lines . " done\r";
			if(substr($line, 0, 1) ne "N") { next; } 
			@split = split('\t', $line);
			#Get sequences for every gene
			$output = &extract_and_analyze_motifs(substr($split[1],3), ($split[2] - $tss_down), ($split[2] + $tss_up), 0);
			foreach my $strain (keys %{$output}) {
				if($strain eq "reference") { next; }
				@split = split('\),', $output->{$strain});
				#Split at the last ( to get number and mutations and motif name
				my @a = split(/([^\(]+)$/, $split[0]);
				chop $a[0];
				$background{$strain}{$a[0]}{$a[1]}++;
			}
		}
		print STDERR "\n";
	} elsif($bgtagdir ne "") {
		print STDERR "TODO\n";
		#TODO
		exit;
	} else {
		#Use all genes as background
		#Check if there is a saved background
		for(my $i = 0; $i < @strains; $i++) {
			if($strains[$i] eq "reference") { next; }
			$tmp_bgfile = $param->{'motif_background_path'} . $genome . "_" . $strains[$i] . "_complete_bg";
			if(!-e $tmp_bgfile) {
				$all_bg_exists = 1;	
			}
		}
		if($all_bg_exists == 0) {
			for(my $i = 0; $i < @strains; $i++) {
				if($strains[$i] eq "reference") { next; }
					$tmp_bgfile = $param->{'motif_background_path'} . $genome . "_" . $strains[$i] . "_complete_bg";
					print STDERR "Precalculated background exists!\n";
					print STDERR "Reading in...\n";
					#Read in bg file
					open FH, "<$tmp_bgfile";
					foreach my $line (<FH>) {
						chomp $line;
						@split = split('\t', $line);
						my @motif = split('/', $split[0]);
						$motif[0] =~ s/ //g;
						$background{$strains[$i]}{$motif[0]}{$split[1]} = $split[2];
					}
				}
		} else {
			print STDERR "Background does not exist so far!\n";
			print STDERR "Genrating the background might take several hours.\n";
			print STDERR "Are you sure you want to proceed?\n";
			print STDERR "Hit Ctrl + C for interruption\n";
			print STDERR "Waiting for 10 seconds\n";
			for(my $i = 0; $i < 10; $i++) {
				print STDERR " . ";
				sleep(1);
			}
			print STDERR "\n";
			open BG_SEQ, ">tmp_bgseqs.txt";
			$del{"tmp_bgseqs.txt"} = 1;
			#Generate bg
			$number_of_tss = keys %tss;
			foreach my $line (keys %tss) {
				$count++;
				@split = split('\t', $tss{$line}->{'tss'});
				&extract_seqs_bg($split[0], ($split[1] - $tss_down), ($split[1] + $tss_up));
				print STDERR "Status: " . $count . " vs " . $number_of_tss . " Completed \r";
			}
			print STDERR "\n";
			print STDERR "Scanning sequences for motifs\n";
			my $command = "homer2 find -i tmp_bgseqs.txt -m all.motifs > tmp_bgseqs_with_motifs.txt";
			`$command`;
			$del{"tmp_bgseqs_with_motifs.txt"} = 1;
			print STDERR "Calculate background\n";
			&calculate_bg();
			foreach my $strain (keys %background) {
				$tmp_bgfile = $param->{'motif_background_path'} . $genome . "_" . $strain . "_complete_bg";
				open BG, ">$tmp_bgfile";
				foreach my $m (keys %{$background{$strain}}) {
					$sum = 0;
					@split = split('/', $m);
					$m = $split[0];
					foreach my $motif_position (keys %{$background{$strain}{$m}}) {
						print BG $m . "\t" . $motif_position . "\t" . $background{$strain}{$m}{$motif_position} . "\n";
						$sum += $background{$strain}{$m}{$motif_position};		
					}
					print BG $m . "\t" . 0 . "\t" . ($number_of_tss - $sum) . "\n";
				}
				close BG;
			}
		}
	}
}

#Function to overlap gene epxression patterns with motif pattern
#It uses the first datapoint in the file as reference and calcuates foldchange between the experiemnts (TODO add better method for that here)
#When 2fold different it adds the strain to the candidate pool for up or downregulated genes
#Then it analyzes gain or loss of motifs that follow the up or downregulated pattern of genes
sub analyze_patterns {
	$_ = 0 for my($cand_in_up_are_equal, $cand_in_down_are_equal, $first_cand_in_up, $first_cand_in_down, $occ_in_motif_for_up, $occ_in_motif_for_down, $rest_is_diff_up, $rest_is_diff_down);
	@split = split('\t', $_[0]);
#	print $split[0] . "\n";
	$pot = 0;
	if(@strains != (@split - 8)) {
		print STDERR "Number of strains does not match number of entries in file!\n";
		exit;
	} else {
	
		%cand_motif = ();
		$pattern_up = "";
		$pattern_down = "";
		%strains_down_for_stats = ();
		%strains_up_for_stats = ();
		#Set every gene expression to 16, so we do not include artefacts
		#TODO think about that
		if($split[8] < 16) { $split[8] = 16; }
		@candidate_strains_up = ();
		@candidate_strains_down = ();
		@candidate_strains_equal = ();
		for(my $i = 9; $i < @split; $i++) {
			if($split[$i] < 16) {
				$split[$i] = 16;
			}
			#Gene is upregulated
			#TODO add statistis to that? Use Homer output file?
			if($split[$i]/$split[8] > 2) {
				push(@candidate_strains_up, $strains[$i-8]);
			}
			#Gene is downregulated
			#TODO add statistis to that? Use Homer output file?
			if($split[$i]/$split[8] < 0.5) {
				push(@candidate_strains_down, $strains[$i-8]);
			}
			if($split[$i]/$split[8] < 2 && $split[$i]/$split[8] > 0.5) {
				push(@candidate_strains_equal, $strains[$i-8]);
			}
		}
		push(@candidate_strains_equal, "reference");
		print PROMOTER "\tUP: ";
		for(my $i = 0; $i < @candidate_strains_up; $i++) {
			print PROMOTER $candidate_strains_up[$i] . ",";
			$strain_is_different_up{$candidate_strains_up[$i]} = 1;
		}
		print PROMOTER "\tDOWN: ";
		for(my $i = 0; $i < @candidate_strains_down; $i++) {
			print PROMOTER $candidate_strains_down[$i] . ",";
			$strain_is_different_down{$candidate_strains_down[$i]} = 1;
		}
		my $follows_up = 0;
		my $follows_down = 0;
		my $follows_equal = 0;
		my $p_up = 0;
		my $p_down = 0;
		my $ref;
		my $cand_score;
		if($quant == 0) {
			foreach my $motif (keys %diff_motifs) {
				$follows_up = $follows_down = $follows_equal  = $p_up = $p_down = 0;
			#	print $motif . "\n";
			#	print Dumper %{$save_motifs_aligned{$motif}};
			#	print "\n\n\n";
				$ref = keys %{$save_motifs_aligned{$motif}{'reference'}};
				for(my $i = 0; $i < @candidate_strains_down; $i++) {
					$cand_score = keys %{$save_motifs_aligned{$motif}{$candidate_strains_down[$i]}};
					if($cand_score == $ref) {
						$follows_down = 1;
						last;
					}
					if($i == 0) {
						if($cand_score < $ref) {
							$p_down = -1;
						} elsif($cand_score > $ref) {
							$p_down = 1;
						} else {
							$follows_down = 1;
							last;
						}
					} else {
						if(($cand_score < $ref && $p_down == 1) || ($cand_score > $ref && $p_down == -1)) {
							$follows_down = 1;
							last;
						}
					}
				}
				for(my $i = 0; $i < @candidate_strains_up; $i++) {
					$cand_score = keys %{$save_motifs_aligned{$motif}{$candidate_strains_up[$i]}};
					if($cand_score == $ref) {
						$follows_up = 1;
						last;
					}
					if($i == 0) {
						if($cand_score < $ref) {
							$p_up = -1;
						} elsif($cand_score > $ref) {
							$p_up = 1;
						} else {
							$follows_up = 1;
							last;
						}
					} else {
						if(($cand_score < $ref && $p_up == 1) || ($cand_score > $ref && $p_up == -1)) {
							$follows_up = 1;
							last;
						}
					}
				}
				for(my $i = 0; $i < @candidate_strains_equal; $i++) {
					$cand_score = keys %{$save_motifs_aligned{$motif}{$candidate_strains_equal[$i]}};
					if($cand_score != $ref) {
						$follows_equal = 1;
						last;
					}
				}
				#Now check if we have up and down candidates - do they change their motif in the same direction?
				if($p_up != 0 && $p_down != 0) {
					if(($p_up == -1 && $p_down == -1) || ($p_up == 1 && $p_down == 1)) {
						if(keys %{$save_motifs_aligned{$motif}{$candidate_strains_down[0]}} != keys %{$save_motifs_aligned{$motif}{$candidate_strains_up[0]}}) {
							my @short_name = split('/', $motif);
							$cand_tricky .= $short_name[0] . "(up: " .  ((keys %{$save_motifs_aligned{$motif}{$candidate_strains_up[0]}}) - $ref) . " vs down: " . ((keys %{$save_motifs_aligned{$motif}{$candidate_strains_down[0]}}) - $ref) . "),";
							if(($stats == 0 || $bed_output == 1) && $p_up != 0) {
								$strains_up_for_stats{$motif} = 1;
							}
							if(($stats == 0 || $bed_output == 1) && $p_down != 0) {
								$strains_down_for_stats{$motif} = 1;
							}
						} 
						$follows_up = 1;
						$follows_down = 1;
					}
				}
				if($follows_down == 0 && $follows_up == 0 && $follows_equal == 0) {
					my @short_name = split('/', $motif);
					if($p_up != 0) {
				 		$cand_up .= $short_name[0] . "(" .  ((keys %{$save_motifs_aligned{$motif}{$candidate_strains_up[0]}}) - $ref) . "),";
					}
				#	print "\tdown: "; 
					if($p_down != 0) {
						 $cand_down .= $short_name[0] . "(" . ((keys %{$save_motifs_aligned{$motif}{$candidate_strains_down[0]}}) - $ref) . "),";
					} 
					$pot = 1;
					if(($stats == 0 || $bed_output == 1) && $p_up != 0) {
						$strains_up_for_stats{$motif} = 1;
					}
					if(($stats == 0 || $bed_output == 1) && $p_down != 0) {
						$strains_down_for_stats{$motif} = 1;
					}
				}
			}
		} else {
			foreach my $motif (keys %motif_score) {
				foreach my $pos (keys %{$motif_score{$motif}}) {
					$follows_up = $follows_down = $follows_equal  = $p_up = $p_down = 0;
					$ref = $motif_score{$motif}{$pos}{'reference'};
					for(my $i = 0; $i < @candidate_strains_down; $i++) {
						$cand_score = $motif_score{$motif}{$pos}{$candidate_strains_down[$i]};
						#Check if score different from ref - if not stop
						if(abs($cand_score-$ref) < $motif_diff) {
							$follows_down = 1;
							last;
						}
						#Check the direction of the change
						if($i == 0) {
							if($cand_score < $ref) {
								$p_down = -1;
							} elsif($cand_score > $ref) {
								$p_down = 1;
							} else {
								$follows_down = 1;
								last;
							}
						} else {
							#Check if all candidates have the same direction - if not stop
							if(($cand_score < $ref && $p_down == 1) || ($cand_score > $ref && $p_down == -1)) {
								$follows_down = 1;
								last;
							} 
						}
					}
					for(my $i = 0; $i < @candidate_strains_up; $i++) {
						$cand_score = $motif_score{$motif}{$pos}{$candidate_strains_up[$i]}; 
						if(abs($cand_score-$ref) < $motif_diff) {
							$follows_up = 1;
							last;
						}
						if($i == 0) {
							if($cand_score < $ref) {
								$p_up = -1;
							} elsif($cand_score > $ref) {
								$p_up = 1;
							} else {
								$follows_up = 1;
								last;
							}
						} else {
							if(($cand_score < $ref && $p_up == 1) || ($cand_score > $ref && $p_up == -1)) {
								$follows_up = 1;
							}
						}
					}
					for(my $i = 0; $i < @candidate_strains_equal; $i++) {
						$cand_score = $motif_score{$motif}{$pos}{$candidate_strains_equal[$i]};
						if($cand_score != $ref) {
							$follows_equal = 1;
							last;
						}
					}
					#Now check if we have up and down candidates - do they change their motif in the same direction?
					if($p_up != 0 && $p_down != 0) {
						if(($p_up == -1 && $p_down == -1) || ($p_up == 1 && $p_down == 1)) {
							if($motif_score{$motif}{$pos}{$candidate_strains_down[0]} != $motif_score{$motif}{$pos}{$candidate_strains_up[0]}) {
								my @short_name = split('/', $motif);
								$cand_tricky .= $short_name[0] . "(up: " .  $motif_score{$motif}{$pos}{$candidate_strains_down[0]} . " vs down: " . $motif_score{$motif}{$pos}{$candidate_strains_up[0]} . " vs ref: " . $ref . " - perfect: " . $perfect_motif_score{$motif} . "),";
								if(($stats == 0 || $bed_output == 1) && $p_up != 0) {
						#			print STDERR "TODO stats for motif score\n";
									$strains_up_for_stats{$motif}{$pos} = 1;
								}
								if(($stats == 0 || $bed_output == 1) && $p_down != 0) {
									$strains_down_for_stats{$motif}{$pos} = 1;
						#			print STDERR "TODO stats for motif score\n";
								}
							} 
							$follows_up = 1;
							$follows_down = 1;
						}
					}
					if($follows_down == 0 && $follows_up == 0 && $follows_equal == 0) {
						my @short_name = split('/', $motif);
						if($p_up != 0) {
				 			$cand_up .= $short_name[0] . "(" .  $motif_score{$motif}{$pos}{$candidate_strains_up[0]} . " vs " . $ref . " - perfect: " . $perfect_motif_score{$motif} . " - homer: " . $save_motifs_aligned{$motif}{'reference'}{$pos}{'score'} . "),";
						}
					#	print "\tdown: "; 
						if($p_down != 0) {
							 $cand_down .= $short_name[0] . "(" . $motif_score{$motif}{$pos}{$candidate_strains_down[0]} . " vs " . $ref . " - perfect: " . $perfect_motif_score{$motif} . " - homer: " . $save_motifs_aligned{$motif}{'reference'}{$pos}{'score'} . "),";
						}
						$pot = 1; 
						if(($stats == 0 || $bed_output == 1) && $p_up != 0) {
							$strains_up_for_stats{$motif}{$pos} = 1;
						}
						if(($stats == 0 || $bed_output == 1) && $p_down != 0) {
							$strains_down_for_stats{$motif}{$pos} = 1;
						}
					}
				}
			}
		}
	}
}

#Calculates the z-score for one motif and uses this method to decide if the result is significant
#TODO Probably not the best statistical method to use here - Talk to someone!!!
sub calculate_stats {
	my $diff = $_[0];
	my @strains_diff = @{$_[1]};
	my @a;
	my %seen = ();
	if(@_ > 0) {
		if($diff eq "") {
			return 0;
		}
		$diff = $_[0];
		@strains_diff = @{$_[1]};
	} else {
		@strains_diff = @strains;
	}
	for(my $i = 0; $i < @strains_diff; $i++) {
		if($diff ne "") {
			@a = split('\),', $diff);
		} else {
			@a = split('\),', $output->{$strains_diff[$i]});
		}
		for(my $j = 0; $j < @a; $j++) {
			$a[$j] =~ s/ //g;
			@b = split(/([^\(]+)$/, $a[$j]);
			chop $b[0];
			$m = $b[0];
			#Mean for this motif was not calculated - this motif was not part of the background
			if(!exists $mean{$strains_diff[$i]}{$m}) { next; }
			if($stddev{$strains_diff[$i]}{$m} == 0) {
				$zscore = $b[-1];
			} else {
				$zscore =  (($b[-1]- $mean{$strains_diff[$i]}{$m})/($stddev{$strains_diff[$i]}{$m}));
			}
			
			if($zscore > $zscore_significant_value || $zscore < ($zscore_significant_value * -1)) {
				if(!exists $seen{$m}) {
					print PROMOTER $m . ",";
					$seen{$m}  = 1;
				}
				if(!exists $target{$m}{$strains_diff[$i]}{$b[-1]}) {
					$target{$m}{$strains_diff[$i]}{$b[-1]} = 1;
				} else {
					$target{$m}{$strains_diff[$i]}{$b[-1]}++;
				}
			}
		}
	}
	return 0;
}

#Count the number of lines in a file
sub get_lines_in_file {
	my $number_of_lines = `wc -l $_[0]`;
	return $number_of_lines;
} 

#check if table exists in the database
sub check_if_table_exists{
        $sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
        $sth->execute();
        $num = $sth->fetchrow_hashref();
        $sth->finish();
        return $num->{'c'};
}
