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

$_ = "" for my($output, $data_dir, $genome_dir, $strain);
$_ = () for my(%peaks, %strand, @split, %tree, %lookup_strain, %last_strain, @tmp_split, %save_id, $tree, @files, $tree_tmp);
$_ = 0 for my($hetero, $allele, $line_number, $id);

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-strain <strain>: Strain we look for muts in\n";
        print STDERR "\t-files <files>: Comma seperated list of files\n";
	print STDERR "\t-data_dir <path to strain mutation data>: default defined in config\n";
	print STDERR "\t-genome_dir <path to strain genomes>: default defined in config\n";
	print STDERR "\t-hetero: Data is heterozygous\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my %mandatory = ('-files' => 1, '-strain' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(   	"files=s{,}" => \@files,
		"strain=s" => \$strain,
		"genome_dir=s" => \$genome_dir,
		"data_dir=s" => \$data_dir,
		"hetero" => \$hetero)
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
$strain = uc($strain);

for(my $i = 0; $i < @files; $i++) {
	$files[$i] =~ s/,//g;
}

print STDERR "Read in mutations\n";
my ($tree, $last) = general::read_strains_mut($strain, $data_dir, $allele);

print STDERR "Processing sam file\n";
foreach my $file (@files) {
	open FH, "<$file";
	@split = split("/", $file);
	$output = "";
	for(my $i = 0; $i < @split - 1; $i++) {
		$output .= $split[$i] . "/";
	}
	$output .= "only_muts_" . $split[-1];
	open OUT, ">$output";

	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq "@") {
			print OUT $line . "\n";
		} else {
			@split = split('\t', $line);
			if($split[3] + length($split[9]) > $last->{substr($split[2], 3)}) { next; }
			$tree_tmp = $tree->{substr($split[2], 3)}->fetch($split[3], $split[3] + length($split[9]));	
			if(exists $tree_tmp->[0]->{'mut'}) {
				print OUT $line . "\n";
			}	
		}
	}
	close FH;
	close OUT;
}
