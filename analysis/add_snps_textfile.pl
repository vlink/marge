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
$_ = 0 for my($homer, $format, $chr, $print_header, $exons, $aall, $found, $total, $count, $homo);
$_ = () for my(@strains, @split, $dbh, $sth, $sth_exon, $header, %lookup_header, $mutations);
my $aallele = 2;

my %mandatory = ('-genome' => 1, '-file' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);
my $data = config::read_config()->{'data_folder'};


GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
		"output=s" => \$output,
                "strains=s{,}" => \@strains,
		"exons" => \$exons,
		"all" => \$aall,
                "homer" => \$homer,
		"homo" => \$homo,
                "format" => \$format)
        or die ("Error in command line arguments!\n");


my $aa = `wc -l $file`;
$total = ((split('\s+', $aa))[0]);

if($homo == 1) {
	$aallele = 1;
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

if($output eq "") {
	#Check if input filename has a file ending
	@split = split('\.', $file);
	$output = $split[0] . "_snps.txt";
}
open OUT, ">$output";

my @fileHandles;
my $filename;
my $outfile_open = 0;

sub open_filehandles{
	for(my $i = 0; $i < @strains - 1; $i++) {
		$filename = $data . "/" . uc($strains[$i]) . "/chr" . $_[0] . "_allele_" . $_[1] . ".mut";
		print STDERR $filename . "\n";
		print $filename . "\n";
		open my $fh_mut, "<", "$filename" or die "Can't open $filename: $!\n";
		push @fileHandles, $fh_mut;
	}
        $outfile_open = 1;
}

sub close_filehandles{
	for(my $i = 0; $i < @strains - 1; $i++) {
                close $fileHandles[$i];
	}
        undef @fileHandles;
}
my @current_line;
print STDERR "Add automated chromosome query\n";
my %save_peaks = ();
for(my $chr = 1; $chr < 19; $chr++) {
	open FH, "<$file";
	if($chr > 1) {
		&close_filehandles($chr, 1);
	}
	&open_filehandles($chr, 1);
	print STDERR "Parsing through file and extracting mutations from database for chromosome " . $chr . "\n";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq "#") {
			if($chr == 1) {
				print OUT $line . "\n";
			}
			next;
		}
		@split = split($split_sign, $line);
		if($split[1] eq "chr" . $chr) {
			my $start = int($split[2]) * 1;
			my $stop = int($split[3]) * 1;
			$save_peaks{$start}->{'stop'} = $stop;
			$save_peaks{$start}->{'line'} = $line;
		}
	}
	foreach my $start_pos (sort {$a <=> $b} keys %save_peaks) {
		for(my $run = 0; $run < @strains - 1; $run++) {
			$save_peaks{$start_pos}->{'line'} .= "\t";
			if(defined read_file_line($fileHandles[$run])) {
				@current_line = split('\t', read_file_line($fileHandles[$run])) if defined read_file_line($fileHandles[$run]);

			} else {
				next;
			} 
			while($current_line[1] < $save_peaks{$start_pos}->{'stop'}) {
				if($current_line[1] > $start_pos) {
					$save_peaks{$start_pos}->{'line'} .= $current_line[1] . ",";
				}
				if(defined read_file_line($fileHandles[$run])) {
					@current_line = split('\t', read_file_line($fileHandles[$run])) if defined read_file_line($fileHandles[$run]);
				} else {
					last;
				}
			}
		}
	}
	foreach my $start_pos (keys %save_peaks) {
		print OUT $save_peaks{$start_pos}->{'line'} . "\n";		
	}
	%save_peaks = ();
}

sub read_file_line {
        my $fh = shift;
        my @split;

        if ($fh and my $line = <$fh>) {
                chomp $line;
#               return [ split('\t', $line) ];
                return $line;
        }
        return;
}

sub add_mutations{
	my $start = $_[0];
	my $end = $_[1];
	$mutations = "";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $aa = 0; $aa < $aallele; $aa++) {
			$mutations .= "\t";
			if($strains[$i] eq "reference") { next; }
			$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$i] . "_allele_" . ($aa + 1) . " WHERE chr = \'" . $chr . "\' AND pos >= " . $start . " AND pos <= " . $end);
#			print "SELECT * FROM " . $genome . "_mutations_" . $strains[$i] . "_allele_" . ($aa + 1) . " WHERE chr = \'" . $chr . "\' AND pos >= " . $start . " AND pos <= " . $end . "\n";
			$sth->execute();
			while(my $mut = $sth->fetchrow_hashref()) {
				$mutations .= $mut->{'pos'} . ":" . $mut->{'reference'} . "->" . $mut->{'strain'} . ";";
			}
		}
	}
	$sth->finish();
}
$dbh->disconnect();

