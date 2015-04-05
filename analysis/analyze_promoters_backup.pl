#!/usr/bin/perl 

use strict;
use Getopt::Long;
require '../general/config.pm';
#use Data::Dumper;
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
	exit;
}

if(@ARGV < 1) {
        &printCMD();
}

$_ = "" for my ($genome, $table_name, $ann_file, $name, $file, $seqs, $command);
$_ = 0 for my ($homo, $num, $skip_first, $skip_diff, $tss_up, $tss_down, $filter, $html);
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
		"-filter=s" => \$filter)
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
	$sth = $dbh->prepare("SELECT * FROM "  . $genome . "_annotations where annotation = \'P\'");
	$sth->execute();
	while (my $gene = $sth->fetchrow_hashref()) { 
	#	$tss{$gene->{'gene'}} = $gene->{'chr'} . "\t" . $gene->{'start'} . "\t" . $gene->{'stop'};
		$tss{substr($gene->{'gene'}, 0, length($gene->{'gene'}) - 1)}->{'tss'} = $gene->{'chr'} . "\t" . $gene->{'start'} . "\t" . $gene->{'stop'};
		$tss{substr($gene->{'gene'}, 0, length($gene->{'gene'}) - 1)}->{'strand'} = $gene->{'strand'};
	}
} elsif($ann_file ne "") {
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
} else {
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

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
}
if($skip_first == 0) {
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
#Report first differently expressed exon
if($skip_diff == 0) {
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

if($filter > 0) {
	print STDERR "Collect tag counts per gene\n";
	my $command = "analyzeRepeats.pl rna " . $genome . " -count exons -d";
	for(my $i = 0; $i < @dirs; $i++) {
		$command .= " " . $dirs[$i];
	}
	$command .= " > gene_counts.tmp";
	$to_del{'gepne_counts.tmp'} = 1;

	print $command . "\n";
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

%file_genes = ();
if($file ne "") {
	print STDERR "Read in file with genes\n";
	open FH, "<$file";
	foreach my $line (<FH>) {
		@split = split('\t', $line);
		$file_genes{$split[0]} = $line;
	}
	close FH;
} 
else {
	%file_genes = %tss;
}

my $f = 0;
open OUT, ">promoter_analysis.txt";
my %pos;

foreach my $line (keys %tss) {
	chomp $line;
	if($f == 0) {
		print OUT "ID\tchr\tstart\tstop\tstrand\ttss";
		for(my $i = 0; $i < @strains; $i++) {
			print OUT "\ttss mut " . $strains[$i];
		}
		for(my $i = 0; $i < @strains; $i++) {
			print OUT "\t" . "first - " . $strains[$i];
			for(my $j = 0; $j < @strains; $j++) {
				print OUT "\tfirst mut " . $strains[$j];
			}
		}
		for(my $i = 0; $i < @strains - 1; $i++) {
			for(my $j = $i + 1; $j < @strains; $j++) {
				print OUT "\t" . "diff " . $strains[$i] . " vs " . $strains[$j];
				for(my $k = 0; $k < @strains; $k++) {
					print OUT "\tdiff mut " . $strains[$k];
				}
			}
		}
		print OUT "\n";
		$f++;
	}
	if(exists $filter{$line}) { next; }
	chomp $line;
	@split = split('\t', $tss{$line}->{'tss'});
	if(exists $file_genes{$line}) {
		chomp $file_genes{$line};
		print OUT $line . "\t ". "chr" . $split[0] . ":" . ($split[1] - $tss_down) . "-" . ($split[1] + $tss_up) . "\t";
		$pos{$split[0] . "," . ($split[1] - $tss_down) . "," . ($split[1] + $tss_up)} = 1;
		for(my $j = 0; $j < @strains; $j++) {
			if($strains[$j] eq "reference") { next; }
			$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$j] . "_allele_1 WHERE chr = \'" . $split[0] . "\' AND pos >= " . ($split[1] - $tss_down) . " AND pos <= " . ($split[1] + $tss_up));
			$sth->execute();
			while(my $res = $sth->fetchrow_hashref()) {
				print OUT $res->{'pos'} . ":" . $res->{'reference'} . "->" . $res->{'strain'} . ";";
			}
			print OUT "\t";
		}
	} else {
		next;
	}
	if(!exists $first{$line} && !exists $diff{$line}) {
		if(exists $file_genes{$line}) {
			print OUT "\n";
		}
		next;
	}
	for(my $i = 0; $i < @strains; $i++) {
		if(exists $first{$line}->{$strains[$i]}) {
			print OUT "chr" . $first{$line}->{$strains[$i]}->{'exon_chr'} . ":" . ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down) . "-" . ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up) . "\t";
			$pos{$first{$line}->{$strains[$i]}->{'exon_chr'} . "," . ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down) . "," . ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up)} = 1;
			for(my $j = 0; $j < @strains; $j++) {
				if($strains[$j] eq "reference") { next; }
				$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$j] . "_allele_1 WHERE chr = \'" . $first{$line}->{$strains[$i]}->{'exon_chr'} . "\' AND pos >= " . ($first{$line}->{$strains[$i]}->{'exon_start'} - $tss_down) . " AND pos <= " . ($first{$line}->{$strains[$i]}->{'exon_start'} + $tss_up));
				$sth->execute();
				while(my $res = $sth->fetchrow_hashref()) {
					print OUT $res->{'pos'} . ":" . $res->{'reference'} . "->" . $res->{'strain'} . ";";
				}
				print OUT "\t";
			}
		}
	}
	if(!exists $diff{$line}) {
		print OUT "\n";
		next;
	}
	for(my $i = 0; $i < @strains - 1; $i++) {
		if($strains[$i] eq "reference") { next; }
		for(my $j = 1; $j < @strains; $j++) {
			if($strains[$j] eq "reference") { next; }
			if(exists $diff{$line}->{$strains[$i]}->{$strains[$j]}) {
				print OUT "chr" . $diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'} . ":" . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down) . "-" . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up) . "\t";
				$pos{$diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'} . "," . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down) . "," . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up)} = 1;
				for(my $k = 0; $k < @strains; $k++) {
					if($strains[$k] eq "reference") { next; }
					$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$k] . "_allele_1 WHERE chr = \'" . $diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_chr'} . "\' AND pos >= " . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} - $tss_down) . " AND pos <= " . ($diff{$line}->{$strains[$i]}->{$strains[$j]}->{'exon_start'} + $tss_up));
					$sth->execute();
					while(my $res = $sth->fetchrow_hashref()) {
						print OUT $res->{'pos'} . ":" . $res->{'reference'} . "->" . $res->{'strain'} . ";";
					}
					print OUT "\t";
				} 
			}
		}
	}
	print OUT "\n";
}
close OUT;
#Get sequences
my @k;
my %ref;
my %gapless;
my $tmp_seq;
my $length;

print STDERR "Extract sequences for every strain for every region from database\n";
#database_interaction::get_multiple_alignment($split[1], $split[2], $split[0], "+", 0, 0, 0, 0);
foreach my $s (keys %pos) {
	if($s ne "9,45404577,45405577") { next; }
	open OUT, ">seqs.txt";
	@split = split(",", $s);
#	$seqs = database_interaction::get_genomic_seq($split[1], $split[2], $split[0], "+", 0, 0, 0, 0);
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");
	my $ref = database_interaction::get_multiple_alignment($split[1], $split[2], $split[0], "+", 0, 0, 0, 0);
	$length = 0;
	for(my $i = 0; $i < @{$ref}; $i++) {
		if($length < length($ref->[$i])) {
			$length = length($ref->[$i]);
		}
	}
	for(my $i = 0; $i < @{$ref}; $i++) {
		if($i == 0) {
			print "reference\n";
			print $ref->[0] . "\n";
			$seqs{'reference'} = $ref->[0];
			$tmp_seq = $ref->[0];
			$tmp_seq =~ s/-//g;
			print OUT ">reference_" . $s . "\n";
			print OUT $tmp_seq;
			$gapless{'reference'} = $tmp_seq;
		#	for(my $j = length($tmp_seq); $j < $length; $j++) {
		#		print OUT "X";
		#	}
			print OUT "\n";
		} else {
			$seqs{$strains[$i]} = $ref->[$i];
			print $strains[$i] . "\n" . $ref->[$i] . "\n";
			$tmp_seq = $ref->[$i];
			$tmp_seq =~ s/-//g;
			print OUT ">" . $strains[$i] . "_" . $s . "\n";
			print OUT $tmp_seq;
			$gapless{$strains[$i]} = $tmp_seq;
		#	for(my $j = length($tmp_seq); $j < $length; $j++) {
		#		print OUT "X";
		#	}
			print OUT "\n";
		}
	}
	close OUT;

	print STDERR "Search for motifs in these sequences\n";
	#Next step search for motifs
	my $command = "homer2 find -i seqs.txt -m all.motifs > seqs_with_motifs.txt";
	print $command . "\n";
	`$command`;
	#Analyze motifs (work_on_find.pl just smarter)

	print STDERR "readin file\n";
	open FH, "<seqs_with_motifs.txt";

	my %save = ();
	my %save_pos = ();
	my @split;
	my @second;
	my %seq = ();
	my %orientation = ();

	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@second = split("_", $split[0]);
		$seq{$split[3] . "_" . $split[1] . "_" . $second[0]} = $split[2];
		if(exists $save{$second[1]}{$split[3]}{$second[0]}) {
			$save{$second[1]}{$split[3]}{$second[0]}++;
			$save_pos{$second[1]}{$split[3]}{$second[0]} .= "," . $split[1];
			$orientation{$second[1]}{$split[3]}{$second[0]}{$split[1]} = $split[4];
		} else {
			$save{$second[1]}{$split[3]}{$second[0]} = 1;
			$save_pos{$second[1]}{$split[3]}{$second[0]} = $split[1];
			$orientation{$second[1]}{$split[3]}{$second[0]}{$split[1]} = $split[4];
		}
	}
	close FH;

	print STDERR "Extract candidates!\n";
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

	foreach my $pos (keys %save) {
		print "\n\n\n\n";
		$first = 0;
		print STDERR "POS: " . $pos . "\n";
		foreach my $motif (keys $save{$pos}) {
			$equal = 0;
			foreach my $k (sort {$a cmp $b} keys $save{$pos}{$motif}) {
				if($save{$pos}{$motif}{$k} > 10 && $first == 0) { last; }
				#Save motif
				if($first == 0) { 
					$print = $motif; 
				}
				$print .= "\t" . $k . "\t" . $save{$pos}{$motif}{$k};
				#Check if occurrence of motifs is different in the different strains
				if($first == 0) {
					$val = $save{$pos}{$motif}{$k};
				} else {
					if($val != $save{$pos}{$motif}{$k}) { $equal = 1; }
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
	#Create alignments for the motifs
	foreach my $pos (keys %cands) {
		print "pos: " . $pos . "\n";
		@chr = split(",", $pos);
		$middle = $chr[1] + (($chr[2] - $chr[1])/2);
		@name = ();
		$num = 0;
		foreach my $motif (keys $cands{$pos}) {
			print "motif: " . $motif . "\n";
			foreach my $s (keys $save_pos{$pos}{$motif}) {
				print "strain: " . $s . "\n";
				print $gapless{$s} . "\n";	
				$half = int(length($gapless{$s})/2);
				print "half: " . $half . "\n";
				@tmp_array = split(",", $save_pos{$pos}{$motif}{$s});
				%tmp = ();
				for(my $j = 0; $j < @tmp_array; $j++) {
					$tmp{$tmp_array[$j]} = 1;
				}
				$start = 0;
				$motif_start = 0;
				$array[$num] = "";
				$start = 0;
				my $shift;
				foreach my $t (sort {$a <=> $b} keys %tmp) {
					print "pos: " . $pos . "\n";
					print "motif: " . $motif . "\n";
					print "motif postition: " . $t . "\n";
					print "strain: " . $s . "\n";
					print "orientation: " . $orientation{$pos}{$motif}{$s}{$t} . "\n";
					$motif_start = $half + $t;
					print "motif start: " . $motif_start . "\n";
					print "motif from file\n";
					print $seq{$motif . "_" . $t . "_" . $s} . "\n";
					$m_from_file = $seq{$motif . "_" . $t . "_" . $s};
					if($orientation{$pos}{$motif}{$s}{$t} eq "-" && $t < 0) {
						print "now remoge the length of the motif\n";
						$motif_start = $motif_start - (length($seq{$motif . "_" . $t . "_" . $s})-1);
					}
					if($orientation{$pos}{$motif}{$s}{$t} eq "-" && $t > 0) {
						$motif_start = $motif_start - (length($seq{$motif . "_" . $t. "_" . $s})-1);
					}
					print $motif_start . "\n";
					print "MOTIF from seq:\n";
					$m_from_seq = substr($gapless{$s}, $motif_start, length($seq{$motif . "_" . $t . "_" . $s}));
					print substr($gapless{$s}, $motif_start, length($seq{$motif . "_" . $t . "_" . $s})) . "\n";
					print "start: " . $start . "\tmotif_start: " . $motif_start . "\n";
					my $f= 0;
					for(my $l = $start; $l < $motif_start; $l++) {
						$array[$num] .= " ";
						if($f == 0) { print "add spaces\n"; $f++;}
				#		if(substr($seqs{$s}, $l, 1) eq "-") { $motif_start++; $l++ }
					}
					if($motif_start < $start) {
						print "overlapping motifs\n";
						print "just add: " . substr($seq{$motif . "_" . $t . "_" . $s}, $start - $motif_start) . "\n";
						$array[$num] .= substr($seq{$motif . "_" . $t . "_" . $s}, $start - $motif_start);
						$start = $motif_start + length(substr($seq{$motif . "_" . $t . "_" . $s}, $start - $motif_start)); 
					} else {
						$array[$num] .= $m_from_file;
						$start = $motif_start + length($m_from_file) - 1;
					}
					if($m_from_file ne $m_from_seq) {
						print STDERR "ERROR!\n";
						print STDERR "file: " . $m_from_file . "\tseq: " . $m_from_seq . "\n";
						exit;
					}
					print "length of the array: " . length($array[$num]) . "\n";
					$start = $motif_start + length($seq{$motif . "_" . $t . "_" . $s});
					print "SEQ\n" . $gapless{$s} . "\n";
				}
				for(my $l = $start; $l < length($seqs{$s}); $l++) {
					$array[$num] .= " ";
				}
				#Time to add the gaps
				print "\n\ntime to add gaps!\n";
				print $s . "\n";
				print $seqs{$s} . "\n";
				print $gapless{$s} . "\n";
				print $array[$num] . "\n";
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
				print "with gaps\n";
				print $final . "\n";
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
			#	if(length($ref{$strains[$i]}) < 100) {
			#		print OUT "<pre>";
			#		for(my $r = 0; $r < 10; $r++) {
			#			if(length($ref{$strains[$i]}) > $r * 10 + 10) {
			#		#		my @k = split('', substr($ref{'reference'}, $r * 10, 10));
			#		#		my @q = split('', substr($ref{$strains[$i]}, $r * 10, 10));
			#				print OUT substr($ref{$strains[$i]}, $r * 10, 10) . "|";
			#			} else {
			#				print OUT substr($ref{$strains[$i]}, $r * 10);
			#			}
			#		}
			#		$ref{$strains[$i]} = "";
			#	} else {
					print OUT "<pre>";
					for(my $r = 0; $r < 10; $r++) {
					#	my @k = split('', substr($ref_strain, $r * 10, 10));
					#	my @q = split('', substr($ref{$strains[$i]}, $r * 10, 10));
						print OUT substr($ref{$strains[$i]}, $r * 10, 10) . "|";
					}
					print OUT "\t" . $strains[$i] . "\n";
					$ref{$strains[$i]} = substr($ref{$strains[$i]}, 100);
			#	}
			}
			$ref_strain = substr($ref_strain, 100);
			for(my $i = 0; $i < @array; $i++) {
				while(length($array[$i]) > 0 && length($array[$i]) < 100) {
					$array[$i] .= " ";
				}
				#	print OUT "<pre>";
				#	for(my $r = 0; $r < 10; $r++) {
				#		if(length($array[$i]) > $r * 10 + 10) {
				#			print OUT substr($array[$i], $r * 10, 10) . "|";
				#		} else {
				#			print OUT substr($array[$i], $r * 10);
				#		}
				#	}
				#	$array[$i] = "";
				#} else {
					print OUT "<pre>";
					for(my $r = 0; $r < 10; $r++) {
						print OUT substr($array[$i], $r * 10, 10) . "|";
					}
					$array[$i] = substr($array[$i], 100);
			#	}
				print OUT "\t" . $name[$i] . "\n";
			}
		}
		close OUT;
		exit;
	}
	print STDERR "WHAT THE HELL\n";
	exit;
}


#Output in html or different


sub check_if_table_exists{
        $sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
        $sth->execute();
        $num = $sth->fetchrow_hashref();
        $sth->finish();
        return $num->{'c'};
}

