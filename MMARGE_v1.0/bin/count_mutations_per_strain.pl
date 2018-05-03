#!/usr/bin/env perl
use strict;
use warnings;

# Copyright Verena M. Link <vlink@ucsd.edu>
# 
# This file is part of MMARGE
#
# MMARGE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MMARGE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


use Getopt::Long;
use config;
use Data::Dumper;

$_ = () for my(@split, @strains, %save, @name, @filename);
$_ = 0 for my($snps, $indels, $private, $hetero, $snps_homo_homo, $snps_homo_hetero, $snps_hetero_homo, $snps_hetero_hetero, $indels_homo_homo, $indels_homo_hetero, $indels_hetero_homo, $indels_hetero_hetero);
$_ = "" for my ($data, $path);

sub printCMD {
        print STDERR "\nUsage:\n";
        print STDERR "\t-ind <individuals>: Comma separated list to count mutations\n";
	print STDERR "\t\tinput at least two strains to do a comparion\n";
	print STDERR "\t-data_dir: data folder: default specified in config file\n";
	print STDERR "\t-private: Counts number of private mutations\n\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}


my %mandatory = ('-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);


GetOptions(	"ind=s{,}" => \@strains, 
		"data_dir=s" => \$data, 
		"private" => \$private)
        or die(&printCMD());

#Set variables
if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}
if(@strains < 2) {
	print STDERR "\n\nERROR: Not enough individuals specified!\n\n";
	&printCMD();
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

#Save mutations per strain
foreach my $strains (@strains) {
	if(!-e $data . "/" . uc($strains)) {
		print STDERR "No mutation data found for " . $strains . "\n";
		next;
	}
	$path = $data . "/" . uc($strains) . "/*mut";
	my @files = `ls $path 2> /dev/null`;
	print STDERR "Reading in mutations for " . $strains . "\n";
	foreach my $f (@files) {
		chomp $f;
		if(!-e $f) {
			print STDERR "Could not find $f\n";
			exit;
		}
		chomp $f;
		open FH, "<$f" or die "Can't open $f: $!";
		@filename = split('/', $f);
		@name = split("_", substr($filename[-1], 0, length($filename[-1]) - 4));
		if($name[2] > 1) { $hetero = 1; }
		if(@name == 1) {
			@name = split('\.', $filename[-1]);
			$save{$strains}{$name[0]}{$name[2]}{0} = 1;
			next;
		}
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$save{$strains}{$name[0]}{$name[2]}{$split[0]}{'o'} = $split[1];
			$save{$strains}{$name[0]}{$name[2]}{$split[0]}{'m'} = $split[2];
		}
	}
}
if($hetero == 0) {
	if($private == 0) {
		#Run all pairwise comparisons
		print STDERR "Strain comparison\t#SNPs\t#InDels\n";
		for(my $i = 0; $i < @strains - 1; $i++) {
			for(my $j = $i +1 ; $j < @strains; $j++) {
				$snps = 0;
				$indels = 0;
				#Run through all mutations in strain 1 and see if they exist in strain 2
				if(exists $save{$strains[$i]}) {
					foreach my $chr (keys %{$save{$strains[$i]}}) {
						foreach my $pos (keys %{$save{$strains[$i]}{$chr}{'1'}}) {
							#Mutation does not exist - check length and count indel/snp up
							if(!exists $save{$strains[$j]}{$chr}{'1'}{$pos}) {
								if(length($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'}) > 1 || length($save{$strains[$i]}{$chr}{'1'}{$pos}{'m'}) > 1) {
									$indels++;
								} else {
									$snps++;
								}
							#Mutation does exist - do string comparison to see if different - check length and count indel/snp up
							} else {
								if($save{$strains[$j]}{$chr}{'1'}{$pos}{'o'} ne $save{$strains[$i]}{$chr}{'1'}{$pos}{'o'}) {
									$indels++;
								} elsif($save{$strains[$j]}{$chr}{'1'}{$pos}{'m'} ne $save{$strains[$i]}{$chr}{'1'}{$pos}{'m'}) {
									if(length($save{$strains[$i]}{$chr}{'1'}{$pos}{'m'}) == 1 && length($save{$strains[$j]}{$chr}{'1'}{$pos}{'m'}) == 1) {
										$snps++;
									} else {
										$indels++;
									}
								}
							}
						}
					}
				}
				#Check if mutation exists in second strain
				if(exists $save{$strains[$j]}) {
					foreach my $chr (keys %{$save{$strains[$j]}}) {
						foreach my $pos (keys %{$save{$strains[$j]}{$chr}{'1'}}) {
							#Just count up if mutation does not exist in first strain - if it exists in both strains it was already counted up in the loop for strain 1
							if(!exists $save{$strains[$i]}{$chr}{'1'}{$pos}) {
								if(length($save{$strains[$j]}{$chr}{'1'}{$pos}{'o'}) > 1 || length($save{$strains[$j]}{$chr}{'1'}{$pos}{'m'}) > 1) {
									$indels++;
								} else {
									$snps++;
								}
							}
						}
					}
				}
				#Add comma seperators for 1000 pos
				$snps =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
				$indels =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
				print STDERR "" . $strains[$i] . " vs " . $strains[$j] . "\t" . $snps . "\t" . $indels . "\n";
			}
		}
	} else {
		my $exist = 0;
		my %private_snp;
		my %private_indel;
		for(my $i = 0; $i < @strains; $i++) {
			if(exists $save{$strains[$i]}) {
				foreach my $chr (keys %{$save{$strains[$i]}}) {
					foreach my $pos (keys %{$save{$strains[$i]}{$chr}{'1'}}) {
						$exist = 0;
						for(my $j = 0; $j < @strains; $j++) {
							if($i == $j) { next; }
							if(exists $save{$strains[$j]}{$chr}{'1'}{$pos}) {
								$exist = 1;
							}
						}
						if($exist == 0) {
							if(length($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'}) > 1 || length($save{$strains[$i]}{$chr}{'1'}{$pos}{'m'}) > 1) {
								$private_indel{$strains[$i]}++;
							} else {
								$private_snp{$strains[$i]}++;
							}
						}
					}
				}
			}
		}
		print STDERR "Strain\tprivate SNPs\tprivate InDels\n";
		foreach my $strains (keys %save) {
			if(!exists $private_snp{$strains}) { 
				$private_snp{$strains} = 0;
			} elsif(length($private_snp{$strains}) > 3) {
				$private_snp{$strains} =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
			}
			if(!exists $private_indel{$strains}) {
				$private_indel{$strains} = 0;
			} elsif(length($private_indel{$strains}) > 3) {
				$private_indel{$strains} =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
			}
			print STDERR $strains . "\t" . $private_snp{$strains} . "\t" . $private_indel{$strains} . "\n";;
		}
	}
} else {
	#If data is heterozygous summarize data
	my %summary;
	for(my $i = 0; $i < @strains; $i++) {
		if(!exists $save{$strains[$i]}) { next; }
		foreach my $chr (keys %{$save{$strains[$i]}}) {
			foreach my $pos (keys %{$save{$strains[$i]}{$chr}{'1'}}) {
				if(exists $save{$strains[$i]}{$chr}{'1'}{$pos} && exists $save{$strains[$i]}{$chr}{'2'}{$pos}) {
					$summary{$strains[$i]}{$chr}{$pos} = 2;
				} else {
					$summary{$strains[$i]}{$chr}{$pos} = 1;
				}
			}
			if($hetero == 1) {
				foreach my $pos ( keys %{$save{$strains[$i]}{$chr}{'2'}}) {
					if(exists $save{$strains[$i]}{$chr}{'1'}{$pos} && exists $save{$strains[$i]}{$chr}{'2'}{$pos}) {
						$summary{$strains[$i]}{$chr}{$pos} = 2;
					} else {
						$summary{$strains[$i]}{$chr}{$pos} = 1;
					}
				}
			}
		}
	}
	if($private == 0) {
		#Just in comparison to the reference
		print STDERR "Strain comparison\t#SNPs homo - homo\t#SNPS homo - hetero \t#SNPs hetero - homo\t#SNPs hetero - hetero\t#InDels homo - homo\t#InDels homo - hetero \t#InDels hetero - homo\t#InDels hetero -hetero\n";
		for(my $i = 0; $i < @strains - 1; $i++) {
			for(my $j = $i +1 ; $j < @strains; $j++) {
				($snps_homo_homo, $snps_homo_hetero, $snps_hetero_hetero, $snps_hetero_homo) = 0;
				($indels_homo_homo, $indels_homo_hetero, $indels_hetero_hetero, $indels_hetero_homo) = 0;
				#Run through all mutations in strain 1 and see if they exist in strain 2
				if(exists $summary{$strains[$i]}) {
					foreach my $chr (keys %{$summary{$strains[$i]}}) {
						foreach my $pos (keys %{$summary{$strains[$i]}{$chr}}) {
							if(!exists $summary{$strains[$j]}{$chr}{$pos}) {
								if((exists $save{$strains[$i]}{$chr}{'1'}{$pos} && (length($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'}) > 1 || length($save{$strains[$i]}{$chr}{'1'}{$pos}{'m'}) > 1)) || (exists $save{$strains[$i]}{$chr}{'2'}{$pos} && (length($save{$strains[$i]}{$chr}{'2'}{$pos}{'o'}) > 1 || length($save{$strains[$i]}{$chr}{'2'}{$pos}{'m'}) > 1))) {
									if($summary{$strains[$i]}{$chr}{$pos} == 1) {
										$indels_hetero_homo++;
									} else {
										$indels_homo_homo++;
									}
								} else {
									if($summary{$strains[$i]}{$chr}{$pos} == 1) {
										$snps_hetero_homo++;
									} else {
										$snps_homo_homo++;
									}
								}
							} else {
								if((exists $save{$strains[$i]}{$chr}{'1'}{$pos} && exists $save{$strains[$i]}{$chr}{'2'}{$pos} && exists $save{$strains[$j]}{$chr}{'1'}{$pos} && exists $save{$strains[$j]}{$chr}{'2'}{$pos}) && 
									((($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'} eq $save{$strains[$j]}{$chr}{'1'}{$pos}{'o'} && $save{$strains[$i]}{$chr}{'2'}{$pos}{'o'} eq $save{$strains[$j]}{$chr}{'2'}{$pos}{'o'} && $save{$strains[$i]}{$chr}{'1'}{$pos}{'m'} eq $save{$strains[$j]}{$chr}{'1'}{$pos}{'m'} && $save{$strains[$i]}{$chr}{'2'}{$pos}{'m'} eq $save{$strains[$i]}{$chr}{'2'}{$pos}{'m'})) || 
									(($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'} eq $save{$strains[$j]}{$chr}{'2'}{$pos}{'o'} && $save{$strains[$i]}{$chr}{'2'}{$pos}{'o'} eq $save{$strains[$j]}{$chr}{'1'}{$pos}{'o'} && $save{$strains[$i]}{$chr}{'1'}{$pos}{'m'} eq $save{$strains[$j]}{$chr}{'2'}{$pos}{'m'} && $save{$strains[$i]}{$chr}{'2'}{$pos}{'m'} eq $save{$strains[$i]}{$chr}{'1'}{$pos}{'m'})) || 
									(($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'} eq $save{$strains[$j]}{$chr}{'1'}{$pos}{'o'} && $save{$strains[$i]}{$chr}{'2'}{$pos}{'o'} eq $save{$strains[$j]}{$chr}{'2'}{$pos}{'o'} && $save{$strains[$i]}{$chr}{'1'}{$pos}{'m'} eq $save{$strains[$j]}{$chr}{'1'}{$pos}{'m'} && $save{$strains[$i]}{$chr}{'2'}{$pos}{'m'} eq $save{$strains[$i]}{$chr}{'2'}{$pos}{'m'})))) {
									next;
								} else {
									if((exists $save{$strains[$i]}{$chr}{'1'}{$pos}{'o'} && (length($save{$strains[$i]}{$chr}{'1'}{$pos}{'o'}) > 1 || length($save{$strains[$i]}{$chr}{'1'}{$pos}{'m'}) > 1)) || (exists $save{$strains[$i]}{$chr}{'2'}{$pos} && (length($save{$strains[$i]}{$chr}{'2'}{$pos}{'o'}) > 1 || length($save{$strains[$i]}{$chr}{'2'}{$pos}{'m'}) > 1)) || (exists $save{$strains[$j]}{$chr}{'1'}{$pos}{'o'} && (length($save{$strains[$j]}{$chr}{'1'}{$pos}{'o'}) > 1 || length($save{$strains[$j]}{$chr}{'1'}{$pos}{'m'}) > 1)) || (exists $save{$strains[$j]}{$chr}{'2'}{$pos} && (length($save{$strains[$j]}{$chr}{'2'}{$pos}{'o'}) > 1 || length($save{$strains[$j]}{$chr}{'2'}{$pos}{'m'}) > 1))) {
										if($summary{$strains[$i]}{$chr}{$pos} == 1 && $summary{$strains[$j]}{$chr}{$pos} == 1) {
											$indels_homo_homo++;
										} elsif($summary{$strains[$i]}{$chr}{$pos} == 2 && $summary{$strains[$j]}{$chr}{$pos} == 1) {
											$indels_hetero_homo++;
										} elsif($summary{$strains[$i]}{$chr}{$pos} == 1 && $summary{$strains[$j]}{$chr}{$pos} == 2) {
											$indels_homo_hetero++;
										} else {
											$indels_hetero_hetero++;
										}
									} else {
										if($summary{$strains[$i]}{$chr}{$pos} == 1 && $summary{$strains[$j]}{$chr}{$pos} == 1) {
											$snps_homo_homo++;
										} elsif($summary{$strains[$i]}{$chr}{$pos} == 2 && $summary{$strains[$j]}{$chr}{$pos} == 1) {
											$snps_hetero_homo++;
										} elsif($summary{$strains[$i]}{$chr}{$pos} == 1 && $summary{$strains[$j]}{$chr}{$pos} == 2) {
											$snps_homo_hetero++;
										} else {
											$snps_hetero_hetero++;
										}
									}
								}
							}
						}
					}
				}
				print STDERR $strains[$i] . " vs " . $strains[$j] . "\t";
				print STDERR $snps_homo_homo . "\t" . $snps_homo_hetero . "\t" . $snps_hetero_homo . "\t" . $snps_hetero_hetero . "\t";
				print STDERR $indels_homo_homo . "\t" . $indels_homo_hetero . "\t" . $indels_hetero_homo . "\t" . $indels_hetero_hetero . "\n";
			}

		}
	} else {
		print STDERR "Private mutations for heterozygous individuals\n";
		print STDERR "Function will be added shortly\n";
	}
}
