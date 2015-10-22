#!/usr/bin/perl -w

#Script to interact with the database that saves all the mutation data of the mouse strains

use strict;
use DBI;
use Getopt::Long;
use Data::Dumper;
require '../general/config.pm';
require '../general/system_interaction.pm';

#Create database connection
my $dbh;
my $sth;

if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "\t-genome <genome>: Path to folder with fastq files per chromosome\n";
	print STDERR "\t-strains <strains>: one or several strains (comma separated) - all or not specified: all strains\n";
	print STDERR "\tif this parameter is empty the genome is created for all possible strains!\n";
	print STDERR "\t-out <name of the output directory>: Script automatically creates for each strain a directory with this name + strain\n";
	print STDERR "\tIf nothing is specified this script uses multithreading with 1/4 of all cores if possible\n";
	print STDERR "\t-resume: Starts from last successful log point!\n";
	print STDERR "\t-homo: Genome is homozygous\n";
	exit;	
}

$_ = () for my($num, %hash, $genome, $mutations, $mut_ref, $end, $seq, $c, $last_pos, $strain, $name, @strains, @split, @cache, @cache_save);
$_ = "" for my($organism, $strains, $output_dir, $table_name, $dir, $chr, $strain_seq, $output, $file);
$_ = 0 for my ($multithreading, $user_cores, $resume, $first, $homo);
my $allele = 1;
if($homo == 1) {
	$allele = 1;
}

my %mandatory = ('-genome' => 1, '-out' => 1);
my %convert = map { $_ => 1 } @ARGV;

config::check_parameters(\%mandatory, \%convert);


GetOptions(	"genome=s" => \$organism,
		"strains=s" => \@strains,
		"out=s" => \$output_dir,
		"homo" => \$homo)
or die ("Error in command line arguments!\n");

my $data = config::read_config()->{'data_folder'};

if($homo == 0) {
	print STDERR "Heterozygousity is not supported at the moment!\n";
	exit;
}

#Check if genome is a folder:
if(-d $organism) {
	print STDERR "Genome is a directory, so there should be fa files in there\n";
	my @files = `ls $organism/*`;
	foreach my $f (@files) {
		chomp $f;
	&cache_chromosome($f);
		if(@strains > 1) {
			@cache_save = @cache;
		}
		for(my $i = 0; $i < @strains; $i++) {
			if($i > 0) {
				@cache = @cache_save;
			}
			print $f . "\n";
			$file = $data . "/" . uc($strains[$i]) . "/chr" . $chr . "_allele_" . $allele . ".mut";
			if(-e $file) {
				print $file . "\n";
				print "Create vector\n";
				&create_vectors($file);
			}
			print "Done with incooperating the mutations\n";
			#Flaten the array
			$strain_seq = join("", @cache);
			$output = $output_dir . "_" . $strains[$i]; 
			if(!-e $output) {
				`mkdir $output`;
			}
			&write_output($chr);
		}
	}
} else {
	print STDERR "We do not support only one fastq file at this moment!\n";
	print STDERR "We also do not automatically check for HOMER genomes at the moment!\n";
	exit;
}

sub write_output{
	my $output_file = $output . "/chr" . $chr . ".fa";
	open OUT, ">$output_file";
	print OUT ">chr" . $chr . "\n";
	while(length($strain_seq) > 50) {
		print OUT substr($strain_seq, 0, 50) . "\n";
		$strain_seq =~ s/^.{50}//s; 
	}
	print OUT $strain_seq . "\n";
	close OUT;
}

sub cache_chromosome{
	my $file = $_[0];
	open FH, "<$file";
	@cache = ();
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			$chr = substr($line, 4);
			print STDERR "Saving chromsome " . $chr . "\n";
			next;
		}
		@split = split('', $line);
		foreach my $s (@split) {
			push(@cache, $s);
		}
	}
}

sub create_vectors{
	my $file = $_[0];	
	my @split_mut;
	open FH, "<$file";
	foreach my $line (<FH>) {
		chomp $line;
		@split_mut = split('\t', $line);
		if($split_mut[0] > @cache) { last; }
		if(length($split_mut[1]) > 1) {
			for(my $i = 1; $i < length($split_mut[1]); $i++) {
				$cache[$split_mut[0] + $i - 1] = "";
			}
		} else {
			$cache[$split_mut[0] - 1] = $split_mut[2];
		}
	}
	close FH;
}


