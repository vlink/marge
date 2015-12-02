#!/usr/bin/perl


use strict;
use Getopt::Long;

BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
#require '/Users/verenalink/workspace/strains/general/config.pm';
use analysis;
$_ = "" for my($genome, $strain, $file, $output);
$_ = () for my(%peaks, @split, @strains);
$_ = 0 for my($line_number);

$output = "output_seq_strains.txt";
sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
        print STDERR "\t-strain <strain>: Name of strain for which sequences should be gather\n";
        print STDERR "\t-file <file>: File with genomic coordinates relative to reference\n";
	print STDERR "\t-output <file>: Output file (default: output_seq_strains.txt)\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $param = config::read_config();

GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
                "strains=s" => \$strain,
		"-output=s" => \$output)
        or die("Error in command line options!\n");
#First step: Get the sequences for the peaks
my $allele = 1;
my $data = config::read_config()->{'data_folder'};

open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#" || substr($line, 0, 6) eq "PeakID") {
		next;
	}
	@split = split('\t', $line);
	$peaks{substr($split[1], 3)}{$split[2]} = $split[3];
	$line_number++;
}
close FH;

print STDERR "Extracting sequences from strain genomes\n";
$strains[0] = $strain;
my ($seq_ref, $save_local_shift_ref) = analysis::get_seq_for_peaks($output, \%peaks, \@strains, $data, $allele, $line_number, 0);
print STDERR "sequences saved in " . $output. "\n";
