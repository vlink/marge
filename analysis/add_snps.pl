#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use DBI;
require '../general/config.pm';
require '../general/system_interaction.pm';
use Data::Dumper;
	
$|=1;

if(@ARGV < 1) {
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
        print STDERR "\t-genome <genome>\n";
        print STDERR "\t-file <input file>: coordinates have to be the reference coordinates\n";
	print STDERR "\t-output <output name>: default <file name>_snps.txt\n";
        print STDERR "\t-strains: one or several strains (comma separated) - all or not specified: all strains\n";
	print STDERR "\n\nOptions for RNA-Seq:\n";
	print STDERR "\t-exons: Just looks for mutations within the exons\n";
	print STDERR "\t-all: Looks for mutations in exons and introns\n";
	print STDERR "Annotation RNA-Seq esp. with exons takes a while\n";
        print STDERR "\n\nFormat of input file:\n";
        print STDERR "\t-homer: shifts HOMER peak and RNA files (default)\n";
        print STDERR "\t-format: checks formats.txt to find the user defined format called annotation_format\n";
	print STDERR "\t-homo: Strains are homozygous\n";
        exit;
}


$_ = "" for my($genome, $file, $output, $split_sign, $chr_pos, $shift_order, $start, $end, $tmp_start, $transcript);
$_ = 0 for my($homer, $format, $chr, $print_header, $exons, $all, $found, $total, $count, $homo);
$_ = () for my(@strains, @split, $dbh, $sth, $sth_exon, $header, %lookup_header, $mutations);
my $allele = 2;

my %mandatory = ('-genome' => 1, '-file' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
		"output=s" => \$output,
                "strains=s{,}" => \@strains,
		"exons" => \$exons,
		"all" => \$all,
                "homer" => \$homer,
		"homo" => \$homo,
                "format" => \$format)
        or die ("Error in command line arguments!\n");


my $a = `wc -l $file`;
$total = ((split('\s+', $a))[0]);

if($homo == 1) {
	$allele = 1;
}

my $ref_as_strain = 0;
for(my $i = 0; $i < @strains; $i++) {
        chomp $strains[$i];
        $strains[$i] =~ s/,//g;
        if($strains[$i] eq "reference") {
                $ref_as_strain = 1;
        }
}

if(@strains == 0 || $strains[0] eq "all") {
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

if($homer == 0 && $format == 0) {
	$homer = 1;
}

#Get format
open FH, "<../config/formats.txt";

#Read format
foreach my $line (<FH>) {
        chomp $line;
	if($line eq "") { next; }
        @split = split(':', $line);
        if($homer == 1 && $split[0] eq "annotation_homer") {
                ($split_sign, $chr_pos, $shift_order) = config::read_format($split[1]);
        } elsif($homer == 1 && $split[0] eq "header_homer") {
		($header) = config::read_header($split[1]);
	} elsif($format == 1 && $split[0] eq "annotation_format") {
                ($split_sign, $chr_pos, $shift_order) = config::read_format($split[1]);
        } elsif($format == 1 && $split[0] eq "header_format") {
		($header) = config::read_header($split[1]);
	}
}

for(my $i = 0; $i < @{$header}; $i++) {
	if($header->[$i] eq "") {
		print STDERR "Header is not defined for this file!\n";
	}
	$lookup_header{$header->[$i]} = 1;
}

#Read in file
open FH, "<$file";

if($output eq "") {
	#Check if input filename has a file ending
	@split = split('\.', $file);
	$output = $split[0] . "_snps.txt";
}
open OUT, ">$output";

my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

print STDERR "Parsing through file and extracting mutations from database\n";
foreach my $line (<FH>) {
	chomp $line;
	@split = split($split_sign, $line);
	if($print_header == 0) {
		if(exists $lookup_header{$split[0]}) {
			$print_header = 1;
		} else {
			my @a = split('\s+', $line);
			if(exists $lookup_header{$a[0]}) {
				print OUT $line;
			}
				$print_header = 1;
		}
		if($print_header == 1) {
			print OUT $line;
			for(my $i = 0; $i < @strains; $i++) {
				for(my $a = 0; $a < $allele; $a++) {
					print OUT "\t" . $strains[$i] . " (allele " . ($a + 1) . ")";
				}
			}
			print OUT "\n";
		}
		next;
	}
	if(@split == 0) {
		print OUT $line . "\n";
	} elsif(substr($line, 0, 1) eq "#") {
		print OUT $line . "\n";
	} else {
		print OUT $line;
		#We start to look the mutations up
		for(my $i = 0; $i < @{$shift_order}; $i++) {
			if($shift_order->[$i] eq "chr") {
				$chr = substr($split[$i], 3);
			} elsif($shift_order->[$i] eq "start") {
				$start = $split[$i];
			} elsif($shift_order->[$i] eq "end") {
				$end = $split[$i];
			} elsif($shift_order->[$i] eq "transcript") {
				$transcript = $split[$i];
			}
		}
		if($end < $start) {
			$tmp_start = $end;
			$end = $start;
			$start = $tmp_start;
		}
		if($exons == 1) {
			if($transcript eq "") {
				print STDERR "Your file format does not contain transcript information!\n";
				print STDERR "Please revise the format config file!\n";
				exit;
			}
			$sth_exon = $dbh->prepare("SELECT * FROM " . $genome . "_genes WHERE transcript =\'" . $transcript . "\'");
			$sth_exon->execute();
			while(my $e = $sth_exon->fetchrow_hashref()) {
				$found = 1;
				&add_mutations($e->{'start'}, $e->{'stop'});
			}
			$sth_exon->finish();
			if($found == 0) {
				print STDERR $transcript . " not found!\n";
				$found = 0;
				exit;
			}
		} else {
			&add_mutations($start, $end);
		}
		print OUT $mutations . "\n";
		$count++;
		print "Status: " . int(($count/$total)*100) . "% Completed \r";
	}
}

sub add_mutations{
	my $start = $_[0];
	my $end = $_[1];
	$mutations = "";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $a = 0; $a < $allele; $a++) {
			$mutations .= "\t";
			if($strains[$i] eq "reference") { next; }
			$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$i] . "_allele_" . ($a + 1) . " WHERE chr = \'" . $chr . "\' AND pos >= " . $start . " AND pos <= " . $end);
#			print "SELECT * FROM " . $genome . "_mutations_" . $strains[$i] . "_allele_" . ($a + 1) . " WHERE chr = \'" . $chr . "\' AND pos >= " . $start . " AND pos <= " . $end . "\n";
			$sth->execute();
			while(my $mut = $sth->fetchrow_hashref()) {
				$mutations .= $mut->{'pos'} . ":" . $mut->{'reference'} . "->" . $mut->{'strain'} . ";";
			}
		}
	}
	$sth->finish();
}
$dbh->disconnect();

