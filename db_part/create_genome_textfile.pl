#!/usr/bin/perl -w

use strict;
use Getopt::Long;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'};
use config;

if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "\t-genome <genome>: Path to folder with fastq files per chromosome\n";
	print STDERR "\t-strains <strains>: one or several strains (comma separated)\n";
	print STDERR "\t-in <path to folder with strains mutations\n";
	print STDERR "\t-out <name of the output directory>: Script automatically creates for each strain a directory with this name + strain\n";
	print STDERR "\t-homo: Genome is homozygous\n";
	exit;	
}

$_ = () for my($genome, @strains, @genome_files);
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
}

@genome_files = `ls $genome`;
foreach my $g_file (@genome_files) {
	print $g_file;
	$chr = substr($g_file, 3, length($g_file) - 7);
	for(my $s = 0; $s < @strains; $s++) {
		for(my $i = 0; $i < $allele; $i++) {
			#make sure output folder exists
			if(!-e $output . "/" . $strains[$s]) {
				$folder = $output . "/" . $strains[$s];
				`mkdir -p $folder`;
			}
			$filename = $output . "/" . $strains[$s] . "/chr" . $chr . "_allele_" . ($i + 1) . ".fa";
			$mut_file = $input . "/" . $strains[$s] . "/chr" . $chr . "_allele_" . ($i + 1) . ".mut";
			print $filename . "\n";
			print $mut_file. "\n";
			processing::create_genome($chr, $filename, $genome . "/" . $g_file, $mut_file);
		}

	}
}
