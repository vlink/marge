#!/usr/bin/perl 
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
	print STDERR "\t-strains <strains>: Comma-separated list of strains - if not defined only the reference is used - if all is specified all strains are used\n";
	print STDERR "\t-file <file>: Genes for which the promoter should be annotated\n";
	print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-d <tag dirs>: List with tag directories (RNA-Seq)\n";
	print STDERR "\t-ann <file>: File with TSS annotation!\n";
	print STDERR "\t-no-first: Skips annotation of first expressed exon\n";
	print STDERR "\t-no-diff: Skips annotation of first differently expressed exon\n";
	print STDERR "\t-fc: foldchange that defines differently expressed exons (default: 4.0)\n";
	print STDERR "\t-tss <down>,<up>: Extends promoter from TSS or first expressed base (default: -500,500)\n";
	print STDERR "\t\tExample: -tss -500,500: Extends TSS 500 in both directions\n";
	print STDERR "\t\tExample: -tss 500: Extends TSS 500 in both directions\n";
	print STDERR "\nFilter methods:\n";
	print STDERR "\t-filter <count>: Filters out genes with less than <count> tag counts in all exons together\n";
	print STDERR "\nOutput:\n";
	print STDERR "\t-html: Outputs a html file with the sequence and the motifs\n";
	print STDERR "\n\nBackground (Default: all TSS):\n";
	print STDERR "\t-bgfile: <File with RefSeq Identifiers to use as background\n";
	print STDERR "\t-bgtagdir: <reference and strain RNA-Seq tag dir to calculate equally expressed genes\n";
	print STDERR "\n\nRun methods (Default: find pattern and calculate statistics for it)\n";
	print STDERR "\t-no-stats (Does not calculate statistics for results)\n";
	print STDERR "\t-no-pattern (Des not try to find a pattern in motif mutations that overlaps with gene expression pattern\n";
	exit;
}

if(@ARGV < 1) {
        &printCMD();
}

$_ = "" for my ($genome, $table_name, $ann_file, $name, $file, $seqs, $command, $bgfile, $bgtagdir);
$_ = 0 for my ($homo, $num, $skip_first, $skip_diff, $tss_up, $tss_down, $filter, $html, $html_output, $tasks, $stats, $pattern);
my ($dbh, $sth, @strains, @dirs, %tss, %first, %diff, @split, %to_del, @part, %genes, @tss, %filter, %del, %file_genes, %seqs);
my $homo = 1;


my %mandatory = ('-genome' => 1, '-strains' => 1);
my %convert = map { $_ => 1 } @ARGV;
my $fc = 4;

my $param = config::read_config();

config::check_parameters(\%mandatory, \%convert);

GetOptions(     "genome=s" => \$genome,
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
		"-bgfile=s" => \$bgfile,
		"-bgtagdir" => \$bgtagdir,
		"-no-stats" => \$stats, 
		"-no-pattern" => \$pattern)
	or die("Error in command line options!\n");

database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");

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
if($file ne "") {
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
	%file_genes = %tss;
	$tasks = keys %tss;
}

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

#Calculate the background
#Check if alternative background is used
my $background_path = $param->{'homer_path'} ;

#Calculate the z-score of every motif that is mutated between two strains
if($stats == 0) { 
	&read_in_background();
	print STDERR "Calculate mean and stddev\n";
	&cal_mean_and_stddev();
}


my $f = 0;
open PROMOTER, ">promoter_analysis.txt";
my %pos;
&add_header();
my $check = 0;
my $done_tasks = 0;
my %target = ();
my @b;
my $done_tss = 0;
my $m;
my $zscore = 0;
my $cand_down = "";
my $cand_up = "";
my $cand_tricky = "";
my $pattern_up = "";
my $pattern_down = "";
my @candidate_strains_up;
my @candidate_strains_down;
my %strain_is_different_up = ();
my %strain_is_different_down = ();

foreach my $line (keys %tss) {
	if(exists $filter{$line}) { next; }
	$done_tss++;
	@split = split('\t', $tss{$line}->{'tss'});

	#Gene in file - get motif differences
	if(exists $file_genes{$line}) {
	#	print $line . "\n";
		$output = &extract_and_analyze_motifs($split[0], ($split[1] - $tss_down), ($split[1] + $tss_up));
		chomp $file_genes{$line};
		print PROMOTER $file_genes{$line};
		$cand_up = "";
		$cand_down = "";
		$cand_tricky = "";
		$check = 1;
		chomp $file_genes{$line};
		for(my $j = 0; $j < @strains; $j++) {
			if($strains[$j] eq "reference") { next; }
			$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$j] . "_allele_1 WHERE chr = \'" . $split[0] . "\' AND pos >= " . ($split[1] - $tss_down) . " AND pos <= " . ($split[1] + $tss_up));
			$sth->execute();
			print PROMOTER "\t";
			while(my $res = $sth->fetchrow_hashref()) {
				print PROMOTER $res->{'pos'} . ":" . $res->{'reference'} . "->" . $res->{'strain'} . ";";
			}
		}
		#Try to find a pattern in motif mutations that corresponds to the pattern in gene expression between the different strains
		if($pattern == 0) {
			&analyze_patterns($file_genes{$line});
		}
		#Calculate the z-score of every motif
		if($stats == 0) {
			&calculate_stats();		
		}
		if($html_output == 1) {
			&print_html_output();
		}
	}
	#Alternative protmoer needs same logic
	if(exists $first{$line}) {
		for(my $i = 0; $i < @strains; $i++) {
			if(exists $first{$line}->{$strains[$i]}) {
				print PROMOTER "chr" . $first{$line}->{$strains[$i]}->{'exon_chr'} . ":" . ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down) . "-" . ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up) . "\t";
				$pos{$first{$line}->{$strains[$i]}->{'exon_chr'} . "," . ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down) . "," . ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up)} = 1;
				for(my $j = 0; $j < @strains; $j++) {
					if($strains[$j] eq "reference") { next; }
					$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$j] . "_allele_1 WHERE chr = \'" . $first{$line}->{$strains[$i]}->{'exon_chr'} . "\' AND pos >= " . ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down) . " AND pos <= " . ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up));
					$sth->execute();
					print PROMOTER "\t";
					while(my $res = $sth->fetchrow_hashref()) {
						print PROMOTER $res->{'pos'} . ":" . $res->{'reference'} . "->" . $res->{'strain'} . ";";
					}
				}
				&extract_and_analyze_motifs($first{$line}->{$strains[$i]}->{'exon_chr'}, ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down), ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up));
				if($html_output == 1) {
					&print_html_output();
				}
			}
		}
	}
	if(exists $diff{$line}) {
		for(my $i = 0; $i < @strains - 1; $i++) {
			if($strains[$i] eq "reference") { next; }
			for(my $j = 1; $j < @strains; $j++) {
				if($strains[$j] eq "reference") { next; }
				if(exists $diff{$line}->{$strains[$i]}->{$strains[$j]}) {
					print PROMOTER "chr" . $diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'} . ":" . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down) . "-" . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up) . "\t";
					$pos{$diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'} . "," . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down) . "," . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up)} = 1;
					for(my $k = 0; $k < @strains; $k++) {
						if($strains[$k] eq "reference") { next; }
						$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$k] . "_allele_1 WHERE chr = \'" . $diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'} . "\' AND pos >= " . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down) . " AND pos <= " . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up));
						$sth->execute();
						while(my $res = $sth->fetchrow_hashref()) {
							print PROMOTER $res->{'pos'} . ":" . $res->{'reference'} . "->" . $res->{'strain'} . ";";
						}
						print OUT "\t";
					} 
				}
				&extract_and_analyze_motifs($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'},  ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down), ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up));
				if($html_output == 1) {
					&print_html_output();
				}
			}
		}
	}
	if(exists $file_genes{$line} || exists $first{$line} || exists $diff{$line}) {
		print PROMOTER "\n";
		$done_tasks++;
	}
	print STDERR "Status\tdone: $done_tasks\tvs\ttaks: $tasks Completed\r";
}

sub tss_from_db{
	$sth = $dbh->prepare("SELECT * FROM "  . $genome . "_annotations where annotation = \'P\'");
	$sth->execute();
	while (my $gene = $sth->fetchrow_hashref()) { 
	#	$tss{$gene->{'gene'}} = $gene->{'chr'} . "\t" . $gene->{'start'} . "\t" . $gene->{'stop'};
		$tss{substr($gene->{'gene'}, 0, length($gene->{'gene'}) - 1)}->{'tss'} = $gene->{'chr'} . "\t" . $gene->{'start'} . "\t" . $gene->{'stop'};
		$tss{substr($gene->{'gene'}, 0, length($gene->{'gene'}) - 1)}->{'strand'} = $gene->{'strand'};
	}
}

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

sub first_expressed_exon{
	$command = "analyzeRepeats.pl rna " . $genome . " -count exons -d";
	for(my $i = 0; $i < @dirs; $i++) {
		$command .= " " . $dirs[$i];
	}
	$command .= " -noCondensingParts > genes.tmp";
	$to_del{'genes.tmp'} = 1;

	`$command`; 

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

sub diff_expressed_exon{
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
my %cands;
my @pos;
my @chr;
my $middle;
my $start;
my @array;
my @name;
my $num = 0;
my %output = ();
my %save = ();
my %save_pos = ();
my @split;
my @second;
my %seq = ();
my %orientation = ();

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

sub extract_and_analyze_motifs {
#	print STDERR "Extract sequences for every strain for every region from database\n";
	open OUT, ">seqs.txt";
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");
	my $ref = database_interaction::get_multiple_alignment($_[1], $_[2], $_[0], "+", 0, 0, 0, 0);
	$length = 0;
	for(my $i = 0; $i < @{$ref}; $i++) {
		if($length < length($ref->[$i])) {
			$length = length($ref->[$i]);
		}
	}
	for(my $i = 0; $i < @{$ref}; $i++) {
			$seqs{$strains[$i]} = $ref->[$i];
			$tmp_seq = $ref->[$i];
			$tmp_seq =~ s/-//g;
			print OUT ">" . $strains[$i] . "_" . $_[0] . "_" . $_[1] . "_" . $_[2] . "\n";
			print OUT $tmp_seq . "\n";
			if($html_output == 1) {
				$gapless{$strains[$i]} = $tmp_seq;
			}
	}
	close OUT;

	#Next step search for motifs
	my $command = "homer2 find -i seqs.txt -m all.motifs > seqs_with_motifs.txt 2> output_homer.tmp";
	`$command`;
	#Analyze motifs (work_on_find.pl just smarter)
	open FH, "<seqs_with_motifs.txt";
	%save = ();
	%save_pos = ();
	%orientation = ();
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@second = split("_", $split[0]);
		$seq{$split[3] . "_" . $split[1] . "_" . $second[0]} = $split[2];
		if(exists $save{$second[1]}{$split[3]}{$second[0]}) {
			$save{$second[1]}{$split[3]}{$second[0]}++;
			if($html_output == 1) {
				$save_pos{$second[1]}{$split[3]}{$second[0]} .= "," . $split[1];
				$orientation{$second[1]}{$split[3]}{$second[0]}{$split[1]} = $split[4];
			}
		} else {
			$save{$second[1]}{$split[3]}{$second[0]} = 1;
			if($html_output == 1) {
				$save_pos{$second[1]}{$split[3]}{$second[0]} = $split[1];
				$orientation{$second[1]}{$split[3]}{$second[0]}{$split[1]} = $split[4];
			}
		}
	}
	close FH;

	for(my $i = 0; $i < @strains; $i++) {
		$output{$strains[$i]} = "";
	}
	foreach my $pos (keys %save) {
		$first = 0;
		foreach my $motif (keys $save{$pos}) {
			$equal = 0;
			for(my $i = 0; $i < @strains; $i++) {
				my $k = $strains[$i];
				if(!exists $save{$pos}{$motif}{$k}) { $save{$pos}{$motif}{$k} = 0; }
			
				if($save{$pos}{$motif}{$k} > 10 && $first == 0) { last; }
				#Check if occurrence of motifs is different in the different strains
				if($first == 0) {
					$val = $save{$pos}{$motif}{$k};
				} else {
					if($val != $save{$pos}{$motif}{$k}) {
						$equal = 1;
						my @short_name = split('/', $motif);
						$output{$k} .= $short_name[0] . "(" . ($val - $save{$pos}{$motif}{$k}) . "), ";
					} else {
						my @short_name = split('/', $motif);
						$output{$k} .= $short_name[0] . "(0),";
					}
				}
				$first++;
			}
			#Save all motifs that have different numbers in the different strains
			if($equal == 1) {
				$cands{$pos}{$motif} = 1;
			}
			$first = 0;
		}
	}
	return \%output;
}
my @tmp_array;
my %tmp;
my %tmp_o;
my $rev;
my $count_forward = 0;
my $count_reverse = 0;
my @r;
my @m;
my $motif_start;
my $half;
my $m_from_seq = "";
my $m_from_file = "";

sub print_html_output{
	#Create alignments for the motifs
	foreach my $pos (keys %cands) {
		@chr = split(",", $pos);
		$middle = $chr[1] + (($chr[2] - $chr[1])/2);
		@name = ();
		$num = 0;
		foreach my $motif (keys $cands{$pos}) {
			foreach my $s (keys $save_pos{$pos}{$motif}) {
				$half = int(length($gapless{$s})/2);
				@tmp_array = split(",", $save_pos{$pos}{$motif}{$s});
				%tmp = ();
				for(my $j = 0; $j < @tmp_array; $j++) {
					$tmp{$tmp_array[$j]} = 1;
				}
				$start = 0;
				$motif_start = 0;
				$array[$num] = "";
				$start = 0;
				foreach my $t (sort {$a <=> $b} keys %tmp) {
					$motif_start = $half + $t;
					$m_from_file = $seq{$motif . "_" . $t . "_" . $s};
					if($orientation{$pos}{$motif}{$s}{$t} eq "-" && $t < 0) {
						$motif_start = $motif_start - (length($seq{$motif . "_" . $t . "_" . $s})-1);
					}
					if($orientation{$pos}{$motif}{$s}{$t} eq "-" && $t > 0) {
						$motif_start = $motif_start - (length($seq{$motif . "_" . $t. "_" . $s})-1);
					}
					$m_from_seq = substr($gapless{$s}, $motif_start, length($seq{$motif . "_" . $t . "_" . $s}));
					for(my $l = $start; $l < $motif_start; $l++) {
						$array[$num] .= " ";
					}
					if($motif_start < $start) {
						$array[$num] .= substr($seq{$motif . "_" . $t . "_" . $s}, $start - $motif_start);
						$start = $motif_start + length(substr($seq{$motif . "_" . $t . "_" . $s}, $start - $motif_start)); 
					} else {
						$array[$num] .= $m_from_file;
						$start = $motif_start + length($m_from_file) - 1;
					}
					$start = $motif_start + length($seq{$motif . "_" . $t . "_" . $s});
				}
				for(my $l = $start; $l < length($seqs{$s}); $l++) {
					$array[$num] .= " ";
				}
				#Time to add the gaps
				my @seq_array = split('', $seqs{$s});
				my @motif_array = split('', $array[$num]);
				my $m_pos = 0;
				my $final = "";
				foreach my $sa (@seq_array) {
					if($sa eq "-") {
						$final .= "-";
					} else {
						$final .= $motif_array[$m_pos];
						$m_pos++;
					}
				}
				$array[$num] = $final;
				$name[$num] = $motif . " " . $s;
				$num++;
			}
		}
		
		%ref = %seqs;
		my $name = $pos . ".html";
		print STDERR $name . "\n";
		open OUT, ">seq_$pos.html";
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
			for(my $i = 0; $i < @array; $i++) {
				while(length($array[$i]) > 0 && length($array[$i]) < 100) {
					$array[$i] .= " ";
				}
				print OUT "<pre>";
				for(my $r = 0; $r < 10; $r++) {
					print OUT substr($array[$i], $r * 10, 10) . "|";
				}
				$array[$i] = substr($array[$i], 100);
				print OUT "\t" . $name[$i] . "\n";
			}
		}
		close OUT;
	}
}

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
			$output = &extract_and_analyze_motifs(substr($split[1],3), ($split[2] - $tss_down), ($split[2] + $tss_up));
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

sub analyze_patterns {
	$_ = 0 for my($cand_in_up_are_equal, $cand_in_down_are_equal, $first_cand_in_up, $first_cand_in_down, $occ_in_motif_for_up, $occ_in_motif_for_down, $rest_is_diff_up, $rest_is_diff_down);
	@split = split('\t', $_[0]);
#	print $split[0] . "\n";
	if(@strains != (@split - 8)) {
		print STDERR "Number of strains does not match number of entries in file!\n";
		exit;
	} else {
		$pattern_up = "";
		$pattern_down = "";
		%strain_is_different_down = ();
		%strain_is_different_up = ();
		#Set every gene expression to 16, so we do not include artefacts
		#TODO think about that
		if($split[8] < 16) { $split[8] = 16; }
		@candidate_strains_up = ();
		@candidate_strains_down = ();
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
		}
		print PROMOTER "\tUP: ";
		for(my $i = 0; $i < @candidate_strains_up; $i++) {
			print PROMOTER $candidate_strains_up[$i] . ",";
#			print "up: " .  $candidate_strains_up[$i] . "\n";
			$strain_is_different_up{$candidate_strains_up[$i]} = 1;
		}
		print PROMOTER "\tDOWN: ";
		for(my $i = 0; $i < @candidate_strains_down; $i++) {
			print PROMOTER $candidate_strains_down[$i] . ",";
#			print "down: " . $candidate_strains_down[$i] . "\n";
			$strain_is_different_down{$candidate_strains_down[$i]} = 1;
		}
		$output = &extract_and_analyze_motifs(substr($split[1], 3), ($split[2] - $tss_down), ($split[2] + $tss_up));
		my %save_motif_per_strain;
		foreach my $strains (keys %{$output}) {
			my @a = split('\),', $output->{$strains});
			for(my $j = 0; $j < @a; $j++) {
				$a[$j] =~ s/ //g;
				@b = split(/([^\(]+)$/, $a[$j]);
				chop $b[0];
				$m = $b[0];
				$save_motif_per_strain{$m}{$strains} = $b[-1];	
			}
		}
		my $pot = 0;
		#Try to find a motif that follows the gene expression pattern
		foreach my $motif (keys %save_motif_per_strain) {
#			print $motif . "\n";
			foreach my $strains (keys %{$save_motif_per_strain{$motif}}) {
#				print $strains . "\t" . $save_motif_per_strain{$motif}{$strains} . "\t";
				#The motif occurrence differs between strain and reference and the gene expression in this strain is either up- or downregulated
				$rest_is_diff_up = $rest_is_diff_down = 1;
				$cand_in_up_are_equal = $cand_in_down_are_equal = 0;
				if($save_motif_per_strain{$motif}{$strains} != 0 && (exists $strain_is_different_up{$strains} || exists $strain_is_different_down{$strains})) {
#					print "here\n";
					#We need following pattern
					#Candidates in up and down have to be different
					#Rest of the strains have to be 0
		
					#There are genes that are upregulated
					if(@candidate_strains_up > 0) {
						$first_cand_in_up = $save_motif_per_strain{$motif}{$candidate_strains_up[0]};
						for(my $i = 0; $i < @candidate_strains_up; $i++) {
							#The occurrence of this motif in this strain differs from the occurrence of this motif in another strain that shows the same gene expression pattern 
							#This motif is not considered as a potenital candidate any more
							if($save_motif_per_strain{$motif}{$candidate_strains_up[$i]} != $first_cand_in_up) {
								$cand_in_up_are_equal = 1;
							}
						}
					} else {
						$cand_in_up_are_equal = 1;
					}
#					print "cand in up are equal 1: " . $cand_in_up_are_equal . "\n";
					#Candidate has 0 differences => is not a candidate
					if($first_cand_in_up == 0) { $cand_in_up_are_equal = 1; }	
					if($cand_in_up_are_equal == 0) {
						$occ_in_motif_for_up = $first_cand_in_up;
					}
					
					#Repeat the whole procedure for the downregulated genes/motifs	
					if(@candidate_strains_down > 0) {
						$first_cand_in_down = $save_motif_per_strain{$motif}{$candidate_strains_down[0]};
						for(my $i = 0; $i < @candidate_strains_down; $i++) {
							if($save_motif_per_strain{$motif}{$candidate_strains_down[$i]} != $first_cand_in_down) {
								$cand_in_down_are_equal = 1;
							}
						}
					} else {
						$cand_in_down_are_equal = 1;
					}
#					print "cand in down are equal 1: " . $cand_in_down_are_equal . "\n";
					#Candidate has 0 differences - is not a candidate
					if($first_cand_in_down == 0) { $cand_in_down_are_equal = 1; }	
					if($cand_in_down_are_equal == 0) {
						$occ_in_motif_for_down = $first_cand_in_down;
					}
					$rest_is_diff_up = 0;
					$rest_is_diff_down = 0;
					#Candidates are different - not a candidate
					if($cand_in_down_are_equal == 1) { $rest_is_diff_down = 1; }
					if($cand_in_up_are_equal == 1) { $rest_is_diff_up = 1; }
#					print "do we have a potential candidate?\nrest is diff down: $rest_is_diff_down\nrest is diff up: $rest_is_diff_up\n";
					#Check if we found a potential candidate
					#If so we now have to check that the motif occurrence of the other strains with a different pattern differ from our candidates
					if($cand_in_down_are_equal == 0 || $cand_in_up_are_equal == 0) {
#						print "here2\n";
						if($first_cand_in_up == 0) { $rest_is_diff_up = 1; }
						if($first_cand_in_down == 0) { $rest_is_diff_down = 1; }
						for(my $i = 0; $i < @strains; $i++) {
							if($strains[$i] eq "reference") { next; }
#							print $strains[$i] . "\t" . $save_motif_per_strain{$motif}{$strains[$i]} . "\n";
							#There are strains with up/downregualted gene expression and this strain is not part of this group and the occurrence of this motif for this strain is the same than the occurrence of this motif in the up/downregulated group => This motif can not be considered as a candidate anymore	
							if(keys %strain_is_different_up > 0 && !exists $strain_is_different_up{$strains[$i]} && $save_motif_per_strain{$motif}{$strains[$i]} == $first_cand_in_up) {
								$rest_is_diff_up = 1;
							}
							if(keys %strain_is_different_down > 0 && !exists $strain_is_different_down{$strains[$i]} && $save_motif_per_strain{$motif}{$strains[$i]} == $first_cand_in_down) {
								$rest_is_diff_down = 1;
							}
						}
#						print "rest up: " . $rest_is_diff_up . "\t" . "rest down: " . $rest_is_diff_down . "\n";
											}
#					print "pot: " . $pot . "\n";
				#	last;
				}
			#	last;
			}
#			print "\n";
			if($rest_is_diff_up == 0 && $rest_is_diff_down == 0) {
				#The motif is considered as a candidate for both up and downreguated genes, but the occurrence in the up and downregulated group is either in both positive or in both negative
				#The motif is still considered as a candidate, but it is added to the group of tricky candidates
				if(($first_cand_in_up > 0 && $first_cand_in_down > 0) || ($first_cand_in_up < 0 && $first_cand_in_down < 0)) {
					$cand_tricky .= $motif . "(UP: " . $first_cand_in_up . "; DOWN: " . $first_cand_in_down. "), ";
					$rest_is_diff_up = 1;
					$rest_is_diff_down = 1;	
				} 
			}

			if($rest_is_diff_up == 0) {
				$cand_up .= $motif . "($first_cand_in_up), ";
				$pot = 1;
			}
			if($rest_is_diff_down == 0) {
				$cand_down .= $motif . "($first_cand_in_down), ";
				$pot = 1;
			}
		#	print $rest_is_diff_up . "\t" . $rest_is_diff_down . "\t" . $pot . "\n";
		#	print "cand up: " . $cand_up . "\t\tcand down: " . $cand_down . "\n";
		}
		if($pot == 0) { 
			print PROMOTER "\tno candidate";
		} else {
			print PROMOTER "\tUP candidates: " . $cand_up . "\tDOWN candidates: " . $cand_down . "\tTRICKY candidates: " . $cand_tricky;;
		}
	}
}

sub calculate_stats {
	for(my $i = 0; $i < @strains; $i++) {
		my @a = split('\),', $output->{$strains[$i]});
		for(my $j = 0; $j < @a; $j++) {
			$a[$j] =~ s/ //g;
			@b = split(/([^\(]+)$/, $a[$j]);
			chop $b[0];
			$m = $b[0];
			if(!exists $mean{$strains[$i]}{$m}) { next; }
			if($stddev{$strains[$i]}{$m} == 0) {
				$zscore = $b[-1];
			} else {
				$zscore =  (($b[-1]- $mean{$strains[$i]}{$m})/($stddev{$strains[$i]}{$m}));
			}
			if($zscore > 1 || $zscore < -1) {
				print PROMOTER "sig: " . $m . ",";
				if(!exists $target{$m}{$strains[$i]}{$b[-1]}) {
					$target{$m}{$strains[$i]}{$b[-1]} = 1;
				} else {
					$target{$m}{$strains[$i]}{$b[-1]}++;
				}
			}
		}
	}

}

sub get_lines_in_file {
	my $number_of_lines = `wc -l $_[0]`;
	print $number_of_lines . "\n";
} 

sub check_if_table_exists{
        $sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
        $sth->execute();
        $num = $sth->fetchrow_hashref();
        $sth->finish();
        return $num->{'c'};
}
