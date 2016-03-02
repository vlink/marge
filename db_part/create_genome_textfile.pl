#!/usr/bin/perl -w

use strict;
use Getopt::Long;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'};
use config;
use Storable;

if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "\t-genome <genome>: Path to folder with fastq files per chromosome\n";
	print STDERR "\t-strains <strains>: one or several strains (comma separated)\n";
	print STDERR "\t-in <path to folder with strains mutations\n";
	print STDERR "\t-out <name of the output directory>: Script automatically creates for each strain a directory with this name + strain\n";
	print STDERR "\t-homo: Genome is homozygous\n";
	exit;	
}

$_ = () for my($genome, @strains, @genome_files, @file, %lookup, %lookup_strain);
$_ = "" for my($chr, $filename, $mut_file, $homo, $input, $output, $folder, $allele);

GetOptions(	"genome=s" => \$genome,
		"strains=s{,}" => \@strains,
		"out=s" => \$output,
		"in=s" => \$input,
		"homo" => \$homo)
or die ("Error in command line arguments!\n");

if($homo == 1) {
	$allele = 1;
} else {
	$allele = 2;
}


if($output eq "") {
	$output = config::read_config()->{'data_folder'};
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
	if(!-e $output . "/" . $strains[$i]) {
		$folder = $output . "/" . $strains[$i];
		`mkdir -p $folder`;
	}
}

@genome_files = `ls $genome/*fa`;
foreach my $g_file (@genome_files) {
	chomp $g_file;
	@file = split("/", $g_file);
	$chr = substr($file[-1], 3, length($file[-1]) - 6);
	for(my $i = 1; $i <= $allele; $i++) {
		processing::create_genome($chr, \@strains, $input, $output, $g_file, $allele);
	}
}
