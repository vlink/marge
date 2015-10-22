#!/usr/bin/perl -w

use strict;
require '/data/home/vlink/mouse_strains/marge/db_part/database_interaction.pm';
use Getopt::Long;
require '../general/config.pm';
use Data::Dumper;

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
	print STDERR "\t-peakfile: Peakfile\n";
	print STDERR "\t-output <output> (Default: peakfile_sequences.txt)\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - if not defined only the reference is used - if all is specified all strains are used\n";
	print STDERR "\t-homo: homozygouse (default: heterzygouse) TODO\n";
	print STDERR "\t-format: <peak|bed> (default: peak)\n"; 
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

$_ = "" for my ($genome, $peakfile, $output, $return, $format);
$_ = 0 for my($homo);
$_ = () for my(@strains);

my %mandatory = ('-peakfile' => 1, '-genome' => 1);
my %convert = map { $_ => 1 } @ARGV;

config::check_parameters(\%mandatory, \%convert);

GetOptions(	"genome=s" => \$genome,
		"peakfile=s" => \$peakfile,
		"strains=s{,}" => \@strains,
		"output=s"=> \$output,
		"homo" => \$homo,
		"format=s" => \$format)
	or die("Error in command line options!\n");



database_interaction::set_global_variables(\@strains, $genome, $homo, 0, "genomic");

#Strains 1
my $ref_as_strain = 0;
for(my $i = 0; $i < @strains; $i++) {
	chomp $strains[$i];
	$strains[$i] =~ s/,//g;
	if($strains[$i] eq "reference") {
		$ref_as_strain = 1;
	}
}
if(@strains > 0 && $strains[0] eq "all") {
	my $hash = config::get_strains($genome);
	my $i = 0;
	foreach my $key (keys %{$hash}) {
		$strains[$i] = $key;
		$i++;
		if($key eq "reference") {
			$ref_as_strain = 1;
		}
	}
}
if($ref_as_strain == 0) {
	push(@strains, "reference");
}

open FH, "<$peakfile";
my $total = `wc -l $peakfile`;
$total = (split('\s', $total))[0];

if($output eq "") {
	if(substr($peakfile, length($peakfile) - 4) eq ".txt" || substr($peakfile, length($peakfile) - 4) eq ".bed") {
		$peakfile = substr($peakfile, 0, length($peakfile) - 4);
	}
	$output = $peakfile . "_sequences.txt";
}

open OUT, ">$output";

my @split;
my $count = 0;
if($format eq "peak") {
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq "#") {
			next;
		}
		@split = split('\t', $line);
		print OUT ">$split[1]_$split[2]_$split[3]_$split[4]\n";
		$return = database_interaction::get_genomic_seq($split[2], $split[3], substr($split[1], 3), $split[4], 0, 0, 0, 0);
		@split = split('\n', $return);
		print OUT $split[1] . "\n";
		$count++;
		print "Status: " . int(($count/$total)*100) . "% Completed \r";
	}
} elsif($format eq "bed") {
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		print "Status: " . int(($count/$total)*100) . "% Completed \r";	
		print OUT ">$split[0]_$split[1]_$split[2]_+\n";
		$return = database_interaction::get_genomic_seq($split[1], $split[2], substr($split[0], 3), "+", 0, 0, 0, 0);
		@split = split('\n', $return);
		print OUT $split[1] . "\n";
		$count++;
	}
} else {
	print STDERR "Wrong format!\n";
	exit;
}
close FH;
close OUT;

