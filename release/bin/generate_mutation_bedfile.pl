#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use config;
use general;

$_ = "" for my($path, $output, $chr, $data);
$_ = () for my(@strains, @split, @files);
$_ = 0 for my($exists);

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-ind <individuals>: Comma-separated list of individuals\n";
	print STDERR "\t-data_dir <data folder>: Default specified in config\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
	exit;
}

my %mandatory = ('-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "data_dir=s" => \$data,
                "ind=s{,}" => \@strains)
        or die("Error in command line options!\n");

if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	print STDERR "Collection mutations for " . $strains[$i] . "\n";
	$strains[$i] = uc($strains[$i]);
	$path = $data . "/" . uc($strains[$i]);
	for(my $al = 1; $al <= 2; $al++) {
		@files = `ls $path/*_allele_$al.mut`;
		if(@files > 1) {
			$output = "output_mut_" . uc($strains[$i]) . "_allele_$al.bed";
			$exists = 0;
			if(-e $output) {
				general::wait_10_secs($output);
				$exists = 1;
			}
			if(-e $output . ".gz" && $exists == 0) {
				general::wait_10_secs($output . ".gz");
				`rm $output.gz`;
			}		
			open OUT, ">$output";
			#Iterate through all mutation files 
			print OUT "track name=\"" . $strains[$i] . " allele " . $al . " mutations\"\n";
			foreach my $file (@files) {
				chomp $file;
				if(!-e $file) { next; }
				@split = split("/", $file);
				@split = split("_", $split[-1]);
				$chr = $split[0];
				print STDERR "\t\t\t\t$chr - allele $al\n";
				open FH, "<$file";
				foreach my $line (<FH>) {
					chomp $line;
					@split = split('\t', $line);
					print OUT $chr . "\t" . ($split[0] - 1) . "\t" . (($split[0] + length($split[1])) - 1) ."\t" . $split[1] . "->" . $split[2] . "\n";
				}
			}
			close OUT;
			print STDERR "Compressing file\n";
			`gzip $output`;
			print STDERR "Bed graph file: " . $output . ".gz\n";
		}
	}
}
