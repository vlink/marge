#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use DBI;
require '/data/home/vlink/mouse_strains/software/general/config.pm';
require '/data/home/vlink/mouse_strains/software/general/system_interaction.pm';

#my $chr_file2 = "/data/home/vlink/mouse_strains/software/calculate_snps/panTro4_chr1.fa";
#my $chr_file1 = "/data/home/vlink/mouse_strains/software/db_part/hg19/chr1.fa"; 
my $c;
$_ = "" for my ($chain, $org1, $org2, $org1_start, $org1_end, $org2_start, $org2_end, $org1_chr, $org2_chr, $org1_ws, $org2_ws, $org1_mut, $org2_mut, $seq, $org1_strand, $org2_strand, $organism1, $organism2, $ali1, $ali2);
$_ = 0 for my ($help, $del, $ins, $snp, $conversion, $org1_index, $org2_index, $count, $diff, $run, $count_seq, $skip, $line_end, $pos1, $pos2);
$_ = () for my (@split, $dbh, $sth, $seen);

my %mandatory = ('-chain' => 1, '-org1' => 1, '-org2' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

my %rev;
$rev{'A'} = 'T';
$rev{'T'} = 'A';
$rev{'C'} = 'G';
$rev{'G'} = 'C';
$rev{'N'} = 'N';

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
		"org1=s" => \$organism1,
		"org2=s" => \$organism2,
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

print OUT "#org1_chr\torg1_pos\torg1_seq\torg2_chr\torg2_pos\torg2_seq\n";
my $all_chr_org1 = config::chromosome_number($organism1);
my $all_chr_org2 = config::chromosome_number($organism2);
foreach my $key_org1 (sort {$a cmp $b } keys %{$all_chr_org1}) {
#	$org1 = &get_seq_file($chr_file1);
	$org1 = &get_seq($organism1, $key_org1);
	foreach my $key_org2 (sort {$a cmp $b} keys %{$all_chr_org2}) {
		print STDERR "Organism " . $organism1 . " is currently using chromosome " . $key_org1 . "\n";
		$org2 = &get_seq($organism2, $key_org2);
	#	$org2 = &get_seq_file($chr_file2);
		open FH, "<$chain";
		foreach my $line (<FH>) {
			chomp $line;
			if(substr($line, 0, 5) eq "chain") {
				@split = split('\s+', $line);
				$s = $split[1];
				$skip = 0;
				#Extract information from chain line
				$org1_chr = substr($split[2], 3);
				$org2_chr = substr($split[7], 3);
				if($org1_chr ne $key_org1 || $org2_chr ne $key_org2) {
					$skip = 1;
					next;
				}
				$count_seq++;
				print STDERR "Mutations are generated! (chain $count_seq out of $total)\n";
				$org1_strand = $split[4];
				$org1_start = $split[5];
				$org2_strand = $split[9]; 
				$org2_start = $split[10];
				if($org1_strand eq "-") {
					$org1_start = $split[3] - $split[6] + ($split[6] - $split[5]) - 1; 
				}
				if($org2_strand eq "-") {
					$org2_start = $split[8] - $split[11] + ($split[11] - $split[10]) - 1; 
				#	$org2_start = $split[8] - $split[10] - 1;
				}
				$org1_index = $org1_start;
				$org2_index = $org2_start;

				#Generate sequence for this alignment
				if(@{$org1} == 0) { $skip = 1; }
				if(@{$org2} == 0) { $skip = 1; }
				$run = -1;
				$line_end = $org1_start;
			} elsif($line eq "") {
				next;
			} else {
				if($skip == 1) { next; }
				@split = split('\t', $line);
				if($org1_strand eq "+" && $org2_strand eq "+") {
					for(my $i = 0; $i < $split[0] - 1; $i++) {
						if($org1->[$org1_index] ne $org2->[$org2_index] && $org1->[$org1_index] ne "N" && $org2->[$org2_index] ne "N") {
							print OUT $org1_chr . "," . $org1_strand . "," . $org1_index . "," . $org1->[$org1_index] . "," . $org2_chr . "," . $org2_strand . "," . $org2_index . "," . $org2->[$org2_index] . "\n";
							$snp++;
						}
						$org1_index++;
						$org2_index++;
					}
					if(@split < 3) { next; }
					$org1_mut = "";
					$org2_mut = "";
					#Create alignment
					$pos2 = $org2_index;
					$pos1 = $org1_index;
					$org1_mut .= $org1->[$org1_index];
					$org1_index++;
					$org2_mut .= $org2->[$org2_index];
					$org2_index++;
	
					for(my $gap = 0; $gap < $split[1]; $gap++) {
						$org1_mut .= $org1->[$org1_index];
						$org1_index++;
					}

					if($split[1] > 0 && $split[2] > 0) {
						print OUT $org1_chr . "," . $org1_strand . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $org2_strand . "," . $pos2 . "," . $org2_mut . "\n";
						$pos1++;
						$pos2++;
						$org1_mut = "";
						$org2_mut = "";
					}

					if($split[1] > 0 && $split[2] > 0) {
						$org1_mut .= $org1->[$org1_index];
						$org2_mut .= $org2->[$org2_index];
						$org1_index++;
						$org2_index++;
					}
			
					for(my $gap = 0; $gap < $split[2]; $gap++) {
						$org2_mut .= $org2->[$org2_index];
						$org2_index++;
					}

					print OUT $org1_chr . "," . $org1_strand . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $org2_strand . "," . $pos2 . "," . $org2_mut . "\n";
					if($split[1] > 0 && $split[2] > 0) {
						$org1_index--;
						$org2_index--;
					}

#				} elsif($org1_strand eq "-" && $org2_strand eq "-") {
#					for(my $i = 0; $i < $split[0]; $i++) {
#						if($org1->[$org1_index] ne $org2->[$org2_index] && $org1->[$org1_index] ne "N" && $org2->[$org2_index] ne "N") {
#							print OUT $org1_chr . "\t" . $org1_strand . "\t" . $org1_index . "\t" . $rev{$org1->[$org1_index]} . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $org2_index . "\t" . $rev{$org2->[$org2_index]} . "\n";
#							$snp++;
#						}
#						$org1_index--;
#						$org2_index--;
#					}
#					if(@split < 3) { next; }
#					$org1_mut = "";
#					$org2_mut = "";
#				#	$ali1 = "";
#				#	$ali2 = "";
#					for(my $i = $org1_index + 10; $i > $org1_index + 1; $i--) {
#						$ali1 .= $rev{$org1->[$i]};
#					}
#					for(my $i = $org2_index + 10; $i > $org2_index + 1; $i--) {
#						$ali2 .= $rev{$org2->[$i]};
#					}
#					#Create alignment
#					$pos2 = $org2_index;
#					$pos1 = $org1_index;
#					$org1_mut .= $rev{$org1->[$org1_index]};
#					$ali1 .= $rev{$org1->[$org1_index]};
#					$org1_index--;
#					$org2_mut .= $rev{$org2->[$org2_index]};
#					$ali2 .= $rev{$org2->[$org2_index]};
#					$org2_index--;
#					for(my $gap = 0; $gap < $split[2]; $gap++) {
#						$org2_mut .= $rev{$org2->[$org2_index]};
#						$ali2 .= $rev{$org2->[$org2_index]};
#						$ali1 .= "-";
#						$org2_index--;
#					}
#					if($split[1] > 0 && $split[2] > 0) {
#						print OUT $org1_chr . "\t" . $org1_strand . "\t" . $pos1 . "\t" . $org1_mut . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $pos2 . "\t" . $org2_mut . "\n";
#						$pos1--;
#						$pos2--;
#						$org1_mut = "";
#						$org2_mut = "";
#					}
#					for(my $gap = 0; $gap < $split[1]; $gap++) {
#						$org1_mut .= $rev{$org1->[$org1_index]};
#						$ali1 .= $rev{$org1->[$org1_index]};
#						$ali2 .= "-";
#						$org1_index++;
#					}
#					if($split[1] > 0 && $split[2] > 0) {
#						$org1_mut .= $rev{$org1->[$org1_index]};
#						$ali1 .= $rev{$org1->[$org1_index]};
#						$org2_mut .= $rev{$org2->[$org2_index]};
#						$ali2 .= $rev{$org2->[$org2_index]};
#						$org1_index--;
#						$org2_index--;
#					}
#					for(my $i = $org1_index; $i > $org1_index - 20; $i--) {
#						$ali1 .= $rev{$org1->[$i]};
#					}
#					for(my $i = $org2_index; $i > $org2_index - 20; $i--) {
#						$ali2 .= $rev{$org2->[$i]};
#					}
#
#					print OUT $org1_chr . "\t" . $org1_strand . "\t" . $pos1 . "\t" . $org1_mut . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $pos2 . "\t" . $org2_mut . "\n";
#					print OUT $ali1 . "\n" . $ali2 . "\n";
#					if($split[1] > 0 && $split[2] > 0) {
#						$org1_index++;
#						$org2_index++;
#					}
				} elsif($org1_strand eq "+" && $org2_strand eq "-") {
					for(my $i = 0; $i < $split[0] - 1; $i++) {
						if($org1->[$org1_index] ne $rev{$org2->[$org2_index]} && $org1->[$org1_index] ne "N" && $org2->[$org2_index] ne "N") {
							print OUT $org1_chr . "," . $org1_strand . "," . $org1_index . "," . $org1->[$org1_index] . "," . $org2_chr . "," . $org2_strand . "," . $org2_index . "," . $rev{$org2->[$org2_index]} . "\n";
							$snp++;
						}

						$ali1 .= $org1->[$org1_index];
						$ali2 .= $rev{$org2->[$org2_index]};
						$org1_index++;
						$org2_index--;
					}
					if(@split < 3) { next; }
					$org1_mut = "";
					$org2_mut = "";
					#Create alignment
					$pos2 = $org2_index;
					$pos1 = $org1_index;
					$org1_mut .= $org1->[$org1_index];
					$ali1 .= $org1->[$org1_index];
					$org1_index++;
					$org2_mut .= $rev{$org2->[$org2_index]};
					$ali2 .= $rev{$org2->[$org2_index]};
					$org2_index--;
					for(my $gap = 0; $gap < $split[1]; $gap++) {
						$org1_mut .= $org1->[$org1_index];
						$ali1 .= $org1->[$org1_index];
						$ali2 .= "-";
						$org1_index++;
					}

					if($split[1] > 0 && $split[2] > 0) {
						print OUT $org1_chr . "," . $org1_strand . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $org2_strand . "," . $pos2 . "," . $org2_mut . "\n";
						$pos1++;
						$pos2--;
					#	$org1_index++;
					#	$org2_index--;
						$org1_mut = "";
						$org2_mut = "";
					}
					for(my $gap = 0; $gap < $split[2]; $gap++) {
						$org2_mut .= $rev{$org2->[$org2_index]};
						$ali2 .= $rev{$org2->[$org2_index]};
						$ali1 .= "-";
						$org2_index--;
					}

					if($split[1] > 0 && $split[2] > 0) {
						$org1_mut .= $org1->[$org1_index];
						$org2_mut .= $rev{$org2->[$org2_index]};
						$ali1 .= $org1->[$org1_index];
						$ali2 .= $rev{$org2->[$org2_index]};
					}
					print OUT $org1_chr . "," . $org1_strand . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $org2_strand . "," . $pos2 . "," . $org2_mut . "\n";
				#	if($split[1] > 0 && $split[2] > 0) {
				#		$org1_index--;
				#		$org2_index++;
				#	}
#				} elsif($org1_strand eq "-" && $org2_strand eq "+") {
#					for(my $i = 0; $i < $split[0] - 1; $i++) {
#						if($rev{$org1->[$org1_index]} ne $org2->[$org2_index] && $org1->[$org1_index] ne "N" && $org2->[$org2_index] ne "N") {
#							print OUT $org1_chr . "\t" . $org1_strand . "\t" . $org1_index . "\t" . $rev{$org1->[$org1_index]}. "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $org2_index . "\t" . $org2->[$org2_index] . "\n";
#							$snp++;
#						}
#						$org1_index--;
#						$org2_index++;
#					}
#					if(@split < 3) { next; }
#					$org1_mut = "";
#					$org2_mut = "";
#					#Create alignment
#					$pos2 = $org2_index;
#					$pos1 = $org1_index;
#			#		$ali1 = "";
#			#		$ali2 = "";
#					for(my $i = $org1_index + 10; $i > $org1_index + 1; $i++) {
#						$ali1 .= $rev{$org1->[$i]};
#					}
#					for(my $i = $org2_index - 10; $i > $org2_index - 1; $i--) {
#						$ali2 .= $org2->[$i];
#					}
#					$org1_mut .= $rev{$org1->[$org1_index]};
#					$ali1 .= $rev{$org1->[$org1_index]};
#					$org1_index--;
#					$org2_mut .= $org2->[$org2_index];
#					$ali2 .= $org2->[$org2_index];
#					$org2_index++;
#					for(my $gap = 0; $gap < $split[2]; $gap++) {
#						$org2_mut .= $org2->[$org2_index];
#						$ali2 .= $org2->[$org2_index];
#						$ali1 .= "-";
#						$org2_index++;
#					}
#					print STDERR $line . "\n";
#					print STDERR $split[1] . "\t" . $split[2] . "\n";
#					if($split[1] > 0 && $split[2] > 0) {
#						print OUT $org1_chr . "\t" . $org1_strand . "\t" . $pos1 . "\t" . $org1_mut . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $pos2 . "\t" . $org2_mut . "\n";
#						$pos1--;
#						$pos2++;
#						$org1_mut = "";
#						$org2_mut = "";
#					}
#					print OUT "so far 1: " . $ali1 . "\n";
#					for(my $gap = 0; $gap < $split[1]; $gap++) {
#						$org1_mut .= $rev{$org1->[$org1_index]};
#						$ali1 .= $rev{$org1->[$org1_index]};
#						$ali2 .= "-";
#						$org1_index--;
#					}
#					print OUT "mut 1: " . $org1_mut. "\n";
#					if($split[1] > 0 && $split[2] > 0) {
#						$org1_mut .= $rev{$org1->[$org1_index]};
#						$org2_mut .= $org2->[$org2_index];
#						$ali1 .= $rev{$org1->[$org1_index]};
#						$ali2 .= $org2->[$org2_index];
#						$org1_index--;
#						$org2_index++;
#					}
#					for(my $i = $org1_index; $i > $org1_index - 20; $i++) {
#						$ali1 .= $rev{$org1->[$i]};
#					}
#					for(my $i = $org2_index; $i < $org2_index + 20; $i--) {
#						$ali2 .= $org2->[$i];
#					}
#					print OUT $org1_chr . "\t" . $org1_strand . "\t" . $pos1 . "\t" . $org1_mut . "\t" . $org2_chr . "\t" . $org2_strand . "\t" . $pos2 . "\t" . $org2_mut . "\n";
#		#			print OUT $ali1 . "\n" . $ali2 . "\n";
#					if($split[1] > 0 && $split[2] > 0) {
#						$org1_index--;
#						$org2_index++;
#					}
				} else {
					print STDERR "SOMETHING WEIRD HAPPENED!\n";	
					print STDERR $line . "\n";
					exit;
				}
				if($split[1] > 0) { $del++; }
				if($split[2] > 0) { $ins++; }
			}
		}
		close FH;
	}
}

#my $gap_up = 0;
#my $gap_down = 0;
#for(my $i = 0; $i < length($ali1)/100; $i++) {
#	print OUT "HERE\n";
#	print OUT substr($ali1, $i * 100, 100) . "\n";
#	print OUT substr($ali2, $i * 100, 100) . "\n";
#	if(substr($ali1, $i * 100, 100) eq "----------------------------------------------------------------------------------------------------") {
#		$gap_up++;
#	} else {
#		print OUT "Gap up: " . $gap_up . "\n";
#		$gap_up = 0;
#	}
#	if(substr($ali2, $i * 100, 100) eq "----------------------------------------------------------------------------------------------------") {
#		$gap_down++;
#	} else {
#		print OUT "Gap down: " . $gap_down . "\n";
#		$gap_down = 0;
#	}
#	print OUT "\n\n";
#}

print STDERR "SNPS:\t" . $snp . "\n";
print STDERR "Insertions:\t" . $ins . "\n";
print STDERR "Deletions:\t" . $del . "\n";
close OUT;

sub get_seq {
	my $org = $_[0];
	my @save;
	print STDERR "Organism " . $org . "\n";
	my $chr = $_[1];
	print STDERR "Fetch chromosome " . $_[1] . "\n";
	$sth = $dbh->prepare("SELECT * from " . $org . "_ref_genome where chr = \'" . $chr . "\' ORDER BY pos");
	$sth->execute();
	$count = -1;
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
	@save = split('', $comp_seq);	
	return \@save;
}

sub get_seq_file {
	open FH, "<$_[0]";
	print STDERR "File: " . $_[0] . "\n";
	my @save;
	my $seq = "";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			next;
		}
		$seq .= $line;
	}
	@save = split('', $seq);
	return \@save;
}
