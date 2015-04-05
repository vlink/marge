#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use DBI;
require '/data/home/vlink/mouse_strains/software/general/interspecies_methods.pm';
require '/data/home/vlink/mouse_strains/software/general/config.pm';
use Cwd;

my $dir = getcwd;

#my $chr_file2 = "/data/home/vlink/mouse_strains/software/calculate_snps/panTro4_chr1.fa";
#my $chr_file1 = "/data/home/vlink/mouse_strains/software/db_part/hg19/chr1.fa"; 
my $c;
$_ = "" for my ($chain, $org1, $org2, $org1_start, $org1_end, $org2_start, $org2_end, $org1_chr, $org2_chr, $org1_ws, $org2_ws, $org1_mut, $org2_mut, $seq, $org1_strand, $org2_strand, $organism1, $organism2, $ali1, $ali2, $table_mut, $table_lift, $command, $table_shift);
$_ = 0 for my ($help, $del, $ins, $snp, $conversion, $org1_index, $org2_index, $count, $diff, $run, $count_seq, $skip, $line_end, $pos1, $pos2, $force, $local_shift, $global_shift);
$_ = () for my (@split, $dbh, $sth, $seen, @add_mut, @add_lift);

my $gap = 50;
my %mandatory = ('-chain' => 1, '-org1' => 1, '-org2' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

my %rev;
$rev{'A'} = 'T';
$rev{'T'} = 'A';
$rev{'C'} = 'G';
$rev{'G'} = 'C';
$rev{'N'} = 'N';

sub printCMD{
	print STDERR "\nUsage:\n";
	print STDERR "\t-chain <liftover file from UCSC>\n";
	print STDERR "\t-org1 <First species in liftover file>\n";
	print STDERR "\t-org2 <Second species in liftover file>\n";
	print STDERR "\t-gap: Specifies the maximal length of a gap that is allowed in interspecies alignments (default: 50)\n";
	print STDERR "\t-force: Overwrites table if it already exists! (Default: exit)\n";
	print STDERR "\n\t-h | --help: Shows this help\n";
	exit;
}

GetOptions(     "chain=s" => \$chain,
                "help" => \$help,
		"org1=s" => \$organism1,
		"org2=s" => \$organism2,
		"gap=s" => \$gap,
		"force" => \$force,
		"h" => \$help)
        or die ("Error in command line arguments!\n");

if($help == 1) { &printCMD(); }

print STDERR "Split chain file in smaller pieces.\nGaps greater than $gap are not allowed!\n";
#Split chain file in smaller file so that there are no long gaps
interspecies_methods::split_chain_file($chain, $gap, $dir . "/chain.tmp");
print STDERR "done!\n\nStart reading in genomes!\n";

my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

my $total = `grep chain $chain | wc -l`;
$total =~ s/\s//g;
my $s;

my $all_chr_org1 = config::chromosome_number($organism1);
my $all_chr_org2 = config::chromosome_number($organism2);

my $chr_length_org1 = 0;
my $chr_length_org2 = 0;

foreach my $keys (keys %{$all_chr_org1}) { 
	if(length($keys) > $chr_length_org1) {
		$chr_length_org1 = length($keys);
	}
}
foreach my $keys (keys %{$all_chr_org2}) {
	if(length($keys) > $chr_length_org2) {
		$chr_length_org2 = length($keys);
	}
}

#Generate tables
$table_mut = "mut_" . $organism1 . "_" . $organism2;
if(&check_if_table_exists($table_mut) == 1) {
	print STDERR "Table " . $table_mut . " already exists!\n";
	if($force == 0) {
		print STDERR "Aborting...\n";
#		exit;
	} else {
		print STDERR "Table will be overwritten!\n";
		$dbh->do("DROP TABLE " . $table_mut);
	}
}
#Create mutation table
#$dbh->{ Warn } = 0;
#$command = "CREATE TABLE " . $table_mut . " (chr1 varchar($chr_length_org1) not null, pos1 bigint not null, mut1 varchar(" . ($gap + 1) . "), chr2 varchar($chr_length_org2) not null, pos2 bigint not null, mut2 varchar(" . ($gap+1) . "), PRIMARY KEY (chr1, pos1, chr2, pos2))";
#$command = "CREATE TABLE " . $table_mut . " (chr1 varchar($chr_length_org1) not null, pos1 bigint not null, mut1 varchar(" . ($gap + 1) . "), chr2 varchar($chr_length_org2) not null, pos2 bigint not null, mut2 varchar(" . ($gap+1) . "))";
#$dbh->do($command);
#$command = "CREATE INDEX ON " . $table_mut . "(chr2, pos2)";
#$dbh->do($command);
#$command = "CREATE INDEX ON " . $table_mut . "(chr1, pos1)";
#$dbh->do($command);
#$dbh->{ Warn } = 1;

$table_lift = "lift_" . $organism1 . "_" . $organism2;
if(&check_if_table_exists($table_lift) == 1) {
	print STDERR "Table " . $table_lift . " already exists!\n";
	if($force == 0) {
		print STDERR "Aborting...\n";
#		exit;
	} else {
		print STDERR "Table will be overwritten!\n";
		$dbh->do("DROP TABLE " . $table_lift);
	}
}

#Create lift table
#$dbh->{ Warn } = 0;
#$command = "CREATE TABLE " . $table_lift . " (chr1 varchar($chr_length_org1) not null, start1 bigint not null, end1 bigint not null, chr2 varchar($chr_length_org2) not null, start2 bigint not null, end2 bigint not null, inversion char(1))";
#$dbh->do($command);
#$command = "CREATE INDEX ON " . $table_lift . "(chr1, start1, start1)";
#$dbh->do($command);
#$command = "CREATE INDEX ON " . $table_lift . "(chr1, start1, start2)";
#$dbh->do($command);
#$dbh->{ Warn } = 1;

$table_shift = "offset_" . $organism1 . "_" . $organism2;
if(&check_if_table_exists($table_shift) == 1) {
	print STDERR "Table " . $table_shift . " already exists!\n";
	if($force == 0) {
		print STDERR "Aborting...\n";
#		exit;
	} else {
		print STDERR "Table will be overwritten!\n";
		$dbh->do("DROP TABLE " . $table_shift);
	}
}


#Create offset table
#$dbh->{ Warn } = 0;
#$command = "CREATE TABLE " . $table_shift . " (chr1 varchar($chr_length_org1) not null, start1 bigint not null, chr2 varchar($chr_length_org2), start2 bigint not null, shift bigint not null)";
#$dbh->do($command);
#$command = "CREATE INDEX ON " . $table_shift . " (chr1, start1)";
#$dbh->do($command);
#$command = "CREATE INDEX ON " . $table_shift . " (chr1, start2)";
#$dbh->do($command);
#$dbh->{ Warn } = 1;

my $mut = $dbh->prepare("INSERT INTO $table_mut (chr1, pos1, mut1, chr2, pos2, mut2) VALUES (?, ?, ?, ?, ?, ?)");
my $lift = $dbh->prepare("INSERT INTO $table_lift (chr1, start1, end1, chr2, start2, end2, inversion) VALUES (?, ?, ?, ?, ?, ?, ?)");
my $shift = $dbh->prepare("INSERT INTO $table_shift (chr1, start1, chr2, start2, shift) VALUES (?, ?, ?, ?, ?)");

my $tmp;
my $genome_org2;
#Save chromosomes from organism2
foreach my $key_org2 (sort {$a cmp $b} keys %{$all_chr_org2}) {
	$tmp = &get_seq($organism2, $key_org2);
	$genome_org2->{$key_org2} = $tmp;
}

print STDERR "Organism " . $organism2 . " is fetched\n";
open OUT, ">out.txt";

$dbh->{ AutoCommit } = 0;
#Generate database tables
foreach my $key_org1 (sort {$a cmp $b } keys %{$all_chr_org1}) {
#	$org1 = &get_seq_file($chr_file1);
	$dbh->commit();
	$org1 = &get_seq($organism1, $key_org1);
	foreach my $key_org2 (sort {$a cmp $b} keys %{$all_chr_org2}) {
		print STDERR "Organism " . $organism1 . " is currently using chromosome " . $key_org1 . "\n";
		print STDERR "Organism " . $organism2 . " is currently using chromosome " . $key_org2 . "\n";
	#	$org2 = &get_seq($organism2, $key_org2);
	#	$org2 = &get_seq_file($chr_file2);
		$org2 = $genome_org2->{$key_org2};
		open FH, "<$dir/chain.tmp";
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
				}
				$org1_index = $org1_start;
				$org2_index = $org2_start;

				#Generate sequence for this alignment
				if(@{$org1} == 0) { $skip = 1; }
				if(@{$org2} == 0) { $skip = 1; }
				$local_shift = 0;
				$global_shift = 0;	
				#Generate liftover file
				#1 means negative strand, 0 means positive strand
				print $org1_chr . "," . $split[5] . "," . $split[6] . "," . $org2_chr . "," . $split[10] . "," . $split[11] . ",";
			#	if($org2_strand eq "-") {
			#		$lift->execute($org1_chr, $split[5], $split[6], $org2_chr, $split[10], $split[11], 1);
			#	} else {
			#		$lift->execute($org1_chr, $split[5], $split[6], $org2_chr, $split[10], $split[11], 0);
			#	}
			#	$lift->finish();
				$run = -1;
				$line_end = $org1_start;

			} elsif($line eq "") {
				next;
			} else {
				if($skip == 1) { next; }
				@split = split('\t', $line);
				if($org1_strand eq "+" && $org2_strand eq "+") {
					print OUT "org1: " . $org1_index . "\t"  . "org2: " . $org2_index . "\n";
					for(my $i = $org1_index; $i < $org1_index + 10; $i++) {
						print OUT $org1->[$i];
					}
					print OUT "\t";
					for(my $i = $org2_index; $i < $org2_index + 10; $i++) {
						print OUT $org2->[$i];
					}
					print OUT "\n";
					for(my $i = 0; $i < $split[0] - 1; $i++) {
						if($org1->[$org1_index] ne $org2->[$org2_index] && $org1->[$org1_index] ne "N" && $org2->[$org2_index] ne "N") {
						#	$mut->execute($org1_chr, $org1_index, $org1->[$org1_index], $org2_chr, $org2_index, $org2->[$org2_index]);
						#	$mut->finish();
							print OUT $org1_chr . "," . $org1_index . "," . $org1->[$org1_index] . "," . $org2_chr . "," . $org2_index . "," . $org2->[$org2_index] . "\n";
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
						$local_shift--;
					}

					if($split[1] > 0 && $split[2] > 0) {
						print OUT $org1_chr . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . ","  . $pos2 . "," . $org2_mut . "\n";
				#		$mut->execute($org1_chr, $pos1, $org1_mut, $org2_chr, $pos2, $org2_mut);
				#		$mut->finish();
						$pos1++;
						$pos2++;
						$org1_mut = "";
						$org2_mut = "";
			#		}

			#		if($split[1] > 0 && $split[2] > 0) {
						$org1_mut .= $org1->[$org1_index];
						$org2_mut .= $org2->[$org2_index];
						$org1_index++;
						$org2_index++;
					}
			
					for(my $gap = 0; $gap < $split[2]; $gap++) {
						$org2_mut .= $org2->[$org2_index];
						$org2_index++;
						$local_shift++;
					}
					$global_shift += $local_shift;
					print OUT $org1_chr . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $pos2 . "," . $org2_mut . "\n";
					if($local_shift != 0) {
						print OUT "INDEL\n";
					#	$shift->execute($org1_chr, $pos1, $org2_chr, $pos2, $global_shift);
					#	$shift->finish();
						print OUT $org1_chr . "," . $pos1 . "," . $org2_chr . "," . $pos2 . "," . $global_shift . "\n";
					}
					$local_shift = 0;
			#		$mut->execute($org1_chr, $pos1, $org1_mut, $org2_chr, $pos2, $org2_mut);
			#		$mut->finish();
					if($split[1] > 0 && $split[2] > 0) {
						$org1_index--;
						$org2_index--;
					}
				} elsif($org1_strand eq "+" && $org2_strand eq "-") {
					for(my $i = 0; $i < $split[0] - 1; $i++) {
						if($org1->[$org1_index] ne $rev{$org2->[$org2_index]} && $org1->[$org1_index] ne "N" && $org2->[$org2_index] ne "N") {
							print OUT $org1_chr . "," .  $org1_index . "," . $org1->[$org1_index] . "," . $org2_chr . "," . $org2_index . "," . $rev{$org2->[$org2_index]} . "\n";
					#		$mut->execute($org1_chr, $org1_index, $org1_mut, $org2_chr, $org2_index, $org2_mut);
					#		$mut->finish();
							$snp++;
						}

				#		$ali1 .= $org1->[$org1_index];
				#		$ali2 .= $rev{$org2->[$org2_index]};
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
				#	$ali1 .= $org1->[$org1_index];
					$org1_index++;
					$org2_mut .= $rev{$org2->[$org2_index]};
				#	$ali2 .= $rev{$org2->[$org2_index]};
					$org2_index--;
					for(my $gap = 0; $gap < $split[1]; $gap++) {
						$org1_mut .= $org1->[$org1_index];
				#		$ali1 .= $org1->[$org1_index];
				#		$ali2 .= "-";
						$org1_index++;
					}

					if($split[1] > 0 && $split[2] > 0) {
						print OUT $org1_chr . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $pos2 . "," . $org2_mut . "\n";
				#		$mut->execute($org1_chr, $pos1, $org1_mut, $org2_chr, $pos2, $org2_mut);
				#		$mut->finish();
						$pos1++;
						$pos2--;
					#	$org1_index++;
					#	$org2_index--;
						$org1_mut = "";
						$org2_mut = "";
					}
					for(my $gap = 0; $gap < $split[2]; $gap++) {
						$org2_mut .= $rev{$org2->[$org2_index]};
					#	$ali2 .= $rev{$org2->[$org2_index]};
					#	$ali1 .= "-";
						$org2_index--;
					}

					if($split[1] > 0 && $split[2] > 0) {
						$org1_mut .= $org1->[$org1_index];
						$org2_mut .= $rev{$org2->[$org2_index]};
					#	$ali1 .= $org1->[$org1_index];
					#	$ali2 .= $rev{$org2->[$org2_index]};
					}
			#		$mut->execute($org1_chr, $pos1, $org1_mut, $org2_chr, $pos2, $org2_mut);
			#		$mut->finish();
					print OUT $org1_chr . "," . $pos1 . "," . $org1_mut . "," . $org2_chr . "," . $pos2 . "," . $org2_mut . "\n";

		
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
	print OUT $key_org1 . "\n";
	print OUT "SNPS:\t" . $snp . "\n";
	print OUT "Insertions:\t" . $ins . "\n";
	print OUT "Deletions:\t" . $del . "\n";

}

$dbh->{ AutoCommit } = 1;
close OUT;
#close LIFT;
$dbh->disconnect();

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

sub check_if_table_exists{
	my $table_name = $_[0];
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
	$sth->execute();
	my $num = $sth->fetchrow_hashref();
	$sth->finish();
	return $num->{'c'};
}
