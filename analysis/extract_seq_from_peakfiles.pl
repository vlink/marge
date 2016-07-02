#!/usr/bin/perl
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
use strict;
use Getopt::Long;
use Storable;
use Set::IntervalTree;
use config;
use general;
use analysis_tree;

$_ = "" for my($file, $output, $data_dir, $genome_dir);
$_ = () for my(@strains, %peaks, %strand, @split, %tree, %lookup_strain, %last_strain, @tmp_split, %save_id);
$_ = 0 for my($hetero, $allele, $line_number, $id);

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - Order must overlay with order in annotated peak file\n";
        print STDERR "\t-file <file>: File with genomic coordinates to pull the sequences\n";
	print STDERR "\t-data_dir <path to strain mutation data>: default defined in config\n";
	print STDERR "\t-genome_dir <path to strain genomes>: default defined in config\n";
	print STDERR "\t-output: Name of the output files (default: sequences.txt)\n";
	print STDERR "\t-hetero: Data is heterozygous\n";
	print STDERR "\t-id: Use gene ID: can only be used when only one strain was specified\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my %mandatory = ('-file' => 1, '-strains' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(   	"file=s" => \$file,
		"strains=s{,}" => \@strains,
		"genome_dir=s" => \$genome_dir,
		"data_dir=s" => \$data_dir,
		"hetero" => \$hetero,
		"id" => \$id,
		"output=s" => \$output)
	or die(&printCMD());
#First step: Get the sequences for the peaks

#Set variables
if($hetero == 1) {
	$allele = 2;
} else {
	$allele = 1;
}
if($data_dir eq "") {
	$data_dir = config::read_config()->{'data_folder'};
}
if($genome_dir eq "") {
	$genome_dir = $data_dir;
}
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}
if($output eq "") {
	$output = "sequences.txt";
}

if(@strains > 1 && $id == 1) {
	$id = 0;
}

print STDERR "Saving peaks\n";
open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#" || substr($line, 0, 6) eq "PeakID") {
		next;
	}
	@split = split('\t', $line);
	@tmp_split = split("_", $split[1]);
	$peaks{substr($tmp_split[0], 3)}{$split[2]} = $split[3];
	$strand{substr($tmp_split[0], 3)}{$split[2]} = $split[4];
	$save_id{$split[1] . "_" . $split[2] . "_" . $split[3]} = $split[0];
	$line_number++;
}
close FH;

for(my $a = 1; $a <= $allele; $a++) {
	#Read in strains data
	print STDERR "Loading shift vectors\n";
	for(my $i = 0; $i < @strains; $i++) {
		my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data_dir, $allele, "ref_to_strain");
		$tree{$strains[$i]} = $tree_ref;
		$lookup_strain{$strains[$i]} = $lookup;
		$last_strain{$strains[$i]} = $last;
	}
	#Get sequnecs for every peak
	if($id == 0) {
		analysis::get_seq_for_peaks($output, \%peaks, \@strains, $genome_dir, $a, $line_number, 0, 0, \%tree, \%lookup_strain, \%last_strain, \%strand);
	} else {
		my $tmp = "tmp_" . rand(10) . ".txt";
		my ($seq, $long, $mut) = analysis::get_seq_for_peaks($tmp, \%peaks, \@strains, $genome_dir, $a, $line_number, 0, 0, \%tree, \%lookup_strain, \%last_strain, \%strand);
		`rm $tmp`;
		open OUT, ">$output";
		foreach my $header (keys %{$seq}) {
			@split = split("_", $header);
			print OUT $save_id{$split[0] . "_" . $split[1] . "_" . $split[2]} . "\t" . $seq->{$header} . "\n";
		}
		close OUT;
	}
}
