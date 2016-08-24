#!/usr/bin/perl -w

use strict;
use Getopt::Long;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;

$_ = () for my($genome, @strains, @genome_files, @file, %strains);
$_ = "" for my($chr, $data, $output, $folder, $allele, $ref_genome);
$_ = 0 for my ($hetero);

sub printCMD{ 
	print STDERR "Usage:\n";
	print STDERR "\t-genome <genome>: Path to folder with fastq files per chromosome\n";
	print STDERR "\t-strains <strains>: one or several strains (comma separated)\n";
	print STDERR "\t-data <path to folder with strains mutations> (default defined in config)\n";
	print STDERR "\t-out <name of the output directory>: Script automatically creates a directory for each strain (default in config)\n";
	print STDERR "\t-ref <name for reference>: creates folder for reference genome (if nothing is specified folder is called REFERENCE\n";
	print STDERR "\t-hetero: Genome is heterozygous\n";
	exit;	
}

if(@ARGV < 1) {
	&printCMD();
}

#Check mandatory command line arguments
my %mandatory = ('-genome' => 1, '-strains' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(	"genome=s" => \$genome,
		"strains=s{,}" => \@strains,
		"out=s" => \$output,
		"data=s" => \$data,
		"ref=s" => \$ref_genome,
		"hetero" => \$hetero)
or die (&printCMD());

#Set variables
if($ref_genome eq "") {
	$ref_genome = "REFERENCE";
}
if($hetero == 1) {
	$allele = 2;
} else {
	$allele = 1;
}
if($output eq "") {
	$output = config::read_config()->{'data_folder'};
}
if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}
push(@strains, $ref_genome);

#Generate folder for each strain 
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
	$strains{$strains[$i]} = 1;
	if(!-e $output . "/" . $strains[$i]) {
		$folder = $output . "/" . $strains[$i];
		`mkdir -p $folder`;
	}
}

#Generate genome for each chromsome from the reference for all the strains
@genome_files = `ls $genome/*fa`;
foreach my $g_file (@genome_files) {
	chomp $g_file;
	@file = split("/", $g_file);
	$chr = substr($file[-1], 3, length($file[-1]) - 6);
	processing::create_genome($chr, \%strains, $data, $output, $g_file, $allele);
}
