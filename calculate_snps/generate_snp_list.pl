#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use DBI;
require '/data/home/vlink/mouse_strains/software/general/config.pm';
require '/data/home/vlink/mouse_strains/software/general/system_interaction.pm';

$_ = "" for my ($chain, $org1, $org2, $org1_start, $org1_end, $org2_start, $org2_end, $org1_chr, $org2_chr, $org1_ws, $org2_ws, $org1_mut, $org2_mut, $seq, $org1_strand, $org2_strand);
$_ = 0 for my ($help, $del, $ins, $snp, $conversion, $org1_index, $org2_index, $count, $diff, $run, $count_seq, $skip, $line_end);
$_ = () for my (@split, @org1, @org2, $dbh, $sth, $seen);

my %mandatory = ('-chain' => 1, '-org1' => 1, '-org2' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

sub printCMD() {
	print STDERR "\nUsage:\n";
	print STDERR "\t-chain <liftover file from UCSC>\n";
	print STDERR "\t-org1 <First species in liftover file>\n";
	print STDERR "\t-org2 <Second species in liftover file>\n";
	print STDERR "\n\t-h | --help: Shows this help\n";
	exit;
}

GetOptions(     "chain=s" => \$chain,
                "help" => \$help,
		"org1=s" => \$org1,
		"org2=s" => \$org2,
		"h" => \$help)
        or die ("Error in command line arguments!\n");

if($help == 1) { &printCMD(); }


my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

my $total = `grep chain $chain | wc -l`;
$total =~ s/\s//g;
my $s;

use Data::Dumper;

open OUT, ">vcf_test_file.txt";
$org1 = &get_seq($org1);
#print Dumper $org1;
#print "\n\n";
$org2 = &get_seq($org2);
#print Dumper $org2;
#print "\n\n";

open FH, "<$chain";
print OUT "#org1_chr\torg1_pos\torg1_seq\torg2_chr\torg2_pos\torg2_seq\n";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 5) eq "chain") {
		$count_seq++;
		@split = split('\s+', $line);
		$s = $split[1];
		#Extract information from chain line
		$org1_chr = substr($split[2], 3);
		$org2_chr = substr($split[7], 3);
		$org1_strand = $split[4];
		$org2_strand = $split[9]; 
		$org1_start = $split[5];
		$org2_start = $split[10];
		$org1_index = $org1_start;
		$org2_index = $org2_start;
		#Generate sequence for this alignment
		if(!defined $org1->{$org1_chr}) { $skip = 1; }
		if(!defined $org2->{$org2_chr}) { $skip = 1; }
		$run = -1;
		print STDERR "Sequence are fetched! ($count_seq out of $total)\n";
		$line_end = $org1_start;
	} elsif($line eq "") {
		next;
	} else {
	#	print OUT $line . "\n";
		if($skip == 1) { next; }
		@split = split('\t', $line);
		$line_end += $split[0];
	#	print $line_end . "\n";
		if(@split < 3) { next; }
	#	my $a = "";
	#	my $b = "";
#		print $org1_start . "\n";
#		print $org2_start . "\n";
#		print $org1->{$org1_chr}->[$org1_start] . "\n";
#		print $org2->{$org2_chr}->[$org2_start] . "\n";
		#Start with looping over sequence before alignment
		while($org1_index < $line_end - 1) {
			if($org1->{$org1_chr}->[$org1_index] ne $org2->{$org2_chr}->[$org2_index] && $org1->{$org1_chr}->[$org1_index] ne "N" && $org2->{$org2_chr}->[$org2_index] ne "N") {
				print OUT $org1_chr . "\t" . $org1_strand . "\t" . $org1_index . "\t" . $org1->{$org1_chr}->[$org1_index] . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $org2_index . "\t" . $org2->{$org2_chr}->[$org2_index] . "\n";
				$snp++;
			}
		#	$a .= $org1->{$org1_chr}->[$org1_index];
		#	$b .= $org2->{$org2_chr}->[$org2_index];
			$org1_index++;
			$org2_index++;
		}
	#	print OUT "DONE\n";
	#	print OUT $org1_index . "\n";
	#	print OUT $org2_index . "\n";
	#	print OUT $a . "\n";
	#	print OUT $b . "\n";
		#Set index to current position
		$org1_mut = "";
		$org2_mut = "";
		#Create alignment
		$org1_mut .= $org1->{$org1_chr}->[$org1_index];
		$org1_index++;
		$org2_mut .= $org2->{$org2_chr}->[$org2_index];
		$org2_index++;
		if($org2_index + $split[2] > @{$org2->{$org2_chr}}) {
			print STDERR "We are running out of the string\n";
			print STDERR "Here for org2 \n";
			print STDERR "mutation is " . $line . "\n";
		}
		if($org1_index + $split[1] > @{$org1->{$org1_chr}}) {
			print STDERR "We are running out of the string\n";
			print STDERR "Here for org1 \n";
			print STDERR "mutation is " . $line . "\n";

		}
		for(my $gap = 0; $gap < $split[2]; $gap++) {
			$org2_mut .= $org2->{$org2_chr}->[$org2_index];
			$org2_index++;
	#		print OUT "gap in second\n";
		}
		
		if($split[1] > 0 && $split[2] > 0) {
			print OUT $org1_chr . "\t" . $org1_strand . "\t" . $org1_index . "\t" . $org1_mut . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $org2_index . "\t" . $org2_mut . "\n";
			$org1_index++;
			$org2_index++;
			$org1_mut = "";
			$org2_mut = "";
		}
		for(my $gap = 0; $gap < $split[1]; $gap++) {
			$org1_mut .= $org1->{$org1_chr}->[$org1_index];
			$org1_index++;
	#		print OUT "gap in first\n";
		}
		if($split[1] > 0 && $split[2] > 0) {
			$org1_mut .= $org1->{$org1_chr}->[$org1_index];
			$org2_mut .= $org2->{$org2_chr}->[$org2_index];
			$org1_index++;
			$org2_index++;
		}

		print OUT $org1_chr . "\t" . $org1_strand . "\t" . $org1_index . "\t" . $org1_mut . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $org2_index . "\t" . $org2_mut . "\n";
		if($split[1] > 0 && $split[2] > 0) {
			$org1_index--;
			$org2_index--;
		}
		if($split[1] > 0) { $del++; }
		if($split[2] > 0) { $ins++; }
	}
}
print STDERR "SNPS:\t" . $snp . "\n";
print STDERR "Insertions:\t" . $ins . "\n";
print STDERR "Deletions:\t" . $del . "\n";
close OUT;

sub get_seq {
	my $org = $_[0];
	my $total_chr = config::chromosome_number($org);
	my $save;
	print STDERR "Organism " . $org . "\n";
	foreach my $chr (keys %{$total_chr}) {
		print $chr . "\n";
		$sth = $dbh->prepare("SELECT * from " . $org . "_ref_genome where chr = \'" . $chr . "\' ORDER BY pos");
		$sth->execute();
		$count = 0;
		my $comp_seq = "";
		while(my $fetch = $sth->fetchrow_hashref()) {
			$seq = "";
			while($count + 1 < $fetch->{'pos'}) {
				for(my $i = 0; $i < 100; $i++) {
					$seq .= "N";
				}
				$count++;
			}
			$count++;
			$comp_seq .= $seq;
			$comp_seq .= $fetch->{'seq'};
		}
		print STDERR "Add chr " . $chr . "\n";
		@{$save->{$chr}} = split('', $comp_seq);	
	#	print OUT $comp_seq . "\n";
	}
	return $save;
}
