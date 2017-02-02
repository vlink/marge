#!/usr/bin/perl
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'};
use strict;
use Getopt::Long;
use config;

$_ = () for my(@split, @strains, %save, @name, @filename);
$_ = 0 for my($snps, $indels, $private);
$_ = "" for my ($data, $path);

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-ind <individuals>: Comma separated list to count mutations\n";
	print STDERR "\t-data_dir: data folder: default specified in config file\n";
	print STDERR "\t-private: Counts number of private mutations\n";
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
		@name = split("_", $filename[-1]);
		if(@name == 1) {
			@name = split('\.', $filename[-1]);
			$save{$strains}{$name[0]}{0} = 1;
			next;
		}
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$save{$strains}{$name[0]}{$split[0]} = $split[2];
		}
	}
}
if($private == 0) {
	#Run all pairwise comparisons
	print "Strain comparison\t#SNPs\t#InDels\n";
	for(my $i = 0; $i < @strains - 1; $i++) {
		for(my $j = $i +1 ; $j < @strains; $j++) {
			$snps = 0;
			$indels = 0;
			#Run through all mutations in strain 1 and see if they exist in strain 2
			if(exists $save{$strains[$i]}) {
				foreach my $chr (keys %{$save{$strains[$i]}}) {
					foreach my $pos (keys %{$save{$strains[$i]}{$chr}}) {
						#Mutation does not exist - check length and count indel/snp up
						if(!exists $save{$strains[$j]}{$chr}{$pos}) {
							if(length($save{$strains[$i]}{$chr}{$pos}) > 1) {
								$indels++;
							} else {
								$snps++;
							}
						#Mutation does exist - do string comparison to see if different - check length and count indel/snp up
						} else {
							if($save{$strains[$j]}{$chr}{$pos} ne $save{$strains[$i]}{$chr}{$pos}) {
								if(length($save{$strains[$i]}{$chr}{$pos}) == length($save{$strains[$j]}{$chr}{$pos})) {
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
					foreach my $pos (keys %{$save{$strains[$j]}{$chr}}) {
						#Just count up if mutation does not exist in first strain - if it exists in both strains it was already counted up in the loop for strain 1
						if(!exists $save{$strains[$i]}{$chr}{$pos}) {
							if(length($save{$strains[$j]}{$chr}{$pos}) > 1) {
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
			print "" . $strains[$i] . " vs " . $strains[$j] . "\t" . $snps . "\t" . $indels . "\n";
		}
	}
} else {
	my $exist = 0;
	my %private_snp;
	my %private_indel;
	for(my $i = 0; $i < @strains; $i++) {
		if(exists $save{$strains[$i]}) {
			foreach my $chr (keys %{$save{$strains[$i]}}) {
				foreach my $pos (keys %{$save{$strains[$i]}{$chr}}) {
					$exist = 0;
					for(my $j = 0; $j < @strains; $j++) {
						if($i == $j) { next; }
						if(exists $save{$strains[$j]}{$chr}{$pos}) {
							$exist = 1;
						}
					}
					if($exist == 0) {
						if(length($save{$strains[$i]}{$chr}{$pos}) > 1) {
							$private_indel{$strains[$i]}++;
						} else {
							$private_snp{$strains[$i]}++;
						}
					}
				}
			}
		}
	}
	print "Strain\tprivate SNPs\tprivate InDels\n";
	foreach my $strains (keys %save) {
		$private_snp{$strains} =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
		$private_indel{$strains} =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
		print $strains . "\t" . $private_snp{$strains} . "\t" . $private_indel{$strains} . "\n";;
	}
}

