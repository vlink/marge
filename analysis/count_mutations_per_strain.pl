#!/usr/bin/perl


use strict;
use Getopt::Long;

BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'};
use config;
#use analysis;
$_ = () for my(%peaks, @split, @strains, %save, @name, @filename);
$_ = 0 for my($snps, $indels);
$_ = "" for my ($data, $path);
sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-strains <strains>: Comma separated list to count mutations\n";
	print STDERR "\t-data: data folder: default specified in config file\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $param = config::read_config();

GetOptions(	"strains=s{,}" => \@strains, 
		"data=s" => \$data)
        or die("Error in command line options!\n");

if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}


for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

foreach my $strains (@strains) {
	if(!-e $data . "/" . uc($strains)) {
		next;
	}
	$path = $data . "/" . uc($strains) . "/*mut";
	my @files = `ls $path 2> /dev/null`;
	$snps = 0;
	$indels = 0;
	print STDERR "Reading in mutations for " . $strains . "\n";
	foreach my $f (@files) {
		if(-e $f) {
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

print "Strain comparison\t#SNPs\t#InDels\n";
for(my $i = 0; $i < @strains - 1; $i++) {
	for(my $j = $i +1 ; $j < @strains; $j++) {
		$snps = 0;
		$indels = 0;
		if(exists $save{$strains[$i]}) {
			foreach my $chr (keys %{$save{$strains[$i]}}) {
				foreach my $pos (keys %{$save{$strains[$i]}{$chr}}) {
					if(!exists $save{$strains[$j]}{$chr}{$pos}) {
						if(length($save{$strains[$i]}{$chr}{$pos}) > 1) {
							$indels++;
						} else {
							$snps++;
						}
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
		if(exists $save{$strains[$j]}) {
			foreach my $chr (keys %{$save{$strains[$j]}}) {
				foreach my $pos (keys %{$save{$strains[$j]}{$chr}}) {
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
		$snps =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
		$indels =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
		print "" . $strains[$i] . " vs " . $strains[$j] . "\t" . $snps . "\t" . $indels . "\n";
	}
}
