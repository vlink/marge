#!/usr/bin/perl -w
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
use strict;
use Getopt::Long;
use config;
use processing;
use general;
use analysis_tree;
my $config = config::read_config();
use Data::Dumper;

$_= "" for my($file, $output, $motif_file, $data, $tmp_motif_file, $genome, $genome_dir, $data_dir, $tmp_out, $header);
$_ = 0 for my($hetero, $line_number, $allele, $f);
$_= () for my(@strains, @motif, %delete, %peaks, @split, %tree, %lookup_strain, %last_strain, %motif_analysis, @detail, %lines);

sub printCMD{
	print STDERR "Usage:\n";
	print STDERR "General commands:\n";
	print STDERR "\t-file <input file>: coordinates have to be the reference coordinates\n";
	print STDERR "\t-output <output name>: default <file name>_motifs.txt\n";
	print STDERR "\t-strains: one or several strains (comma separated)\n";
	print STDERR "\t-motif: motif to scan for (if not specified: all motifs are scanned for)\n";
	print STDERR "\t-motif_file: motif file to use (if not specified: config file)\n";
	print STDERR "\t-hetero: Strains are heterozygous\n";
        print STDERR "\t-genome_dir <folder to strains genomes>: default defined in config\n";
        print STDERR "\t-data_dir <folder to data directory>: default defined in config\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-file' => 1, '-strains' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "-file=s" => \$file,
	"-output=s" => \$output,
	"-strains=s{,}" => \@strains,
	"-motif=s{,}" => \@motif,
	"-motif_file=s" => \$motif_file,
	"-hetero" => \$hetero,
	"-genome_dir" => \$genome,
	"-data_dir=s" => \$data) or die (&printCMD());

#First set all variables
if($output eq "") {
	$output = substr($file, 0, length($file) - 4) . "_motifs.txt";
}

if($hetero == 1) {
	$allele = 2;
} else {
	$allele = 1;
}

if($motif_file eq "") {
	$motif_file = config::read_config()->{'motif_file'};
}
for(my $i = 0; $i < @motif; $i++) {
	$motif[$i] =~ s/,//g;
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
}

if($genome_dir eq "" && $data eq "") {
        $genome_dir = config::read_config()->{'data_folder'};
        $data = config::read_config()->{'data_folder'};
}
if($genome_dir eq "") {
        $genome_dir = $data;
}
if($data eq "") {
        $data = config::read_config()->{'data_folder'};
}

#Save motifs if only subset is specified
my ($index_motif_ref, $PWM_ref, $score_ref) = analysis::read_motifs($motif_file);
if(@motif > 0) {
	$tmp_motif_file = "tmp_motif_" . rand(5);
	$delete{$tmp_motif_file} = 1;
	open OUT, ">$tmp_motif_file";
	print STDERR "We need to filter out motifs\n";
	#Save motif files
	for(my $i = 0; $i < @motif; $i++) {
		if(!exists $score_ref->{$motif[$i]}) {
			print STDERR "Motif " . $motif[$i] . " does not exist in motif file!\n";
			next;
		}
		print OUT ">" . $motif[$i] . "\t" . $motif[$i] . "\t" . $score_ref->{$motif[$i]} . "\n";
		foreach my $pos (sort {$a <=> $b} keys %{$PWM_ref->{$motif[$i]}}) {
			print OUT $PWM_ref->{$motif[$i]}->{$pos}->{'A'} . "\t" . $PWM_ref->{$motif[$i]}->{$pos}->{'C'} . "\t" . $PWM_ref->{$motif[$i]}->{$pos}->{'G'} . "\t" . $PWM_ref->{$motif[$i]}->{$pos}->{'T'} . "\n";
		}
	}
	close OUT;
} else {
	$tmp_motif_file  = $motif_file; 
	foreach my $keys (keys %{$PWM_ref}) {
		push @motif, $keys;
	}
}


#Save peaks
open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	if($f == 0) {
		$header = $line;
		$f++;
		next;
	}
	@split = split('\t', $line);
        $peaks{substr($split[1], 3)}{$split[2]} = $split[3];
	$lines{$split[1] . "_" . $split[2] . "_" . $split[3]} = $line;
	$line_number++;
}

if($line_number == 1 || $line_number == 0) {
        print STDERR "File was empty!\n";
        print STDERR "No peaks are saved!\n";
        exit;
}

print STDERR "Loading shift vectors\n";
for(my $i = 0; $i < @strains; $i++) {
        my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data, $allele, "ref_to_strain");
        $tree{$strains[$i]} = $tree_ref;
        $lookup_strain{$strains[$i]} = $lookup;
        $last_strain{$strains[$i]} = $last;
}

print STDERR "Get sequences\n";
$tmp_out = "tmp" . rand(15);
$delete{$tmp_out} = 1;
my ($seq_ref, $l_seq, $filter) = analysis::get_seq_for_peaks($tmp_out, \%peaks, \@strains, $genome_dir, $allele, $line_number, 0, 0, \%tree, \%lookup_strain, \%last_strain);

my $tmp_out_main_motif = $tmp_out . "_scanned.txt";
$delete{$tmp_out . "_scanned.txt"} = 1;
analysis::scan_motif_with_homer($tmp_out, $tmp_out_main_motif, $tmp_motif_file);

#Analysis results
open FH, "<$tmp_out_main_motif";
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	@detail = split("_", $split[0]);
	$motif_analysis{$detail[0] . "_" . $detail[1] . "_" . $detail[2]}{$detail[3]}{$split[3]}++;	
}
close FH;

open OUT, ">$output";
print OUT $header;
for(my $i = 0; $i < @motif; $i++) {
	for(my $j = 0; $j < @strains; $j++) {
		print OUT "\t" . $motif[$i] . "_" . $strains[$j];
	}
}
print OUT "\n";

foreach my $key (keys %lines) {
	print OUT $lines{$key};
	for(my $i = 0; $i < @motif; $i++) {
		for(my $j = 0; $j < @strains; $j++) {
			if(!exists $motif_analysis{$key}{$strains[$j]}{$motif[$i]}) {
				print OUT "\t0";
			} else {
				print OUT "\t" . $motif_analysis{$key}{$strains[$j]}{$motif[$i]};
			}
		}
	}
	print OUT "\n";
}
close OUT;

foreach my $d (keys %delete) {
	`rm $d`;
}
