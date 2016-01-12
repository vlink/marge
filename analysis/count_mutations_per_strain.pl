#!/usr/bin/perl


use strict;
use Getopt::Long;

BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'};
use config;
#use analysis;
$_ = () for my(%peaks, @split, @strains);
$_ = 0 for my($pairwise);

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-strains <strains>: Comma separated list to count mutations\n";
	print STDERR "\t-pairwise: Calculates all pairwise comparisons\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $param = config::read_config();

GetOptions(	"strains=s{,}" => \@strains, 
		"pairwise" => \$pairwise)
        or die("Error in command line options!\n");

my $data = config::read_config()->{'data_folder'};

my $path;
my $snps = 0;
my $indels = 0;
my %save = ();
my @name;

if($pairwise == 0) {
	print "Strain\t#SNPs\t#InDels\n";
} else {
	print "Strain comparison\t#SNPs\t#InDels\n";
}

foreach my $strains (@strains) {
	$strains =~ s/,//g;
	$path = $data . "/" . uc($strains) . "/*mut";
	my @files = `ls $path`;
	$snps = 0;
	$indels = 0;
	if($pairwise == 0) {
		print $strains . "\t";
	} else {
		print STDERR "Reading in mutations for " . $strains . "\n";
	}
	foreach my $f (@files) {
		chomp $f;
		open FH, "<$f";
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			if($pairwise == 0) {
				if(length($split[1]) == length($split[2])) {
					$snps++;
				} else {
					$indels++;
				}
			} else {
				@name = split("_", $f);
				$save{$strains}{$name[0]}{$split[0]} = $split[2];			
			}
		}
	}
	if($pairwise == 0) {
		print $snps . "\t" . $indels . "\n";
	}
}

my %seen;

if($pairwise == 1) {
	for(my $i = 0; $i < @strains - 1; $i++) {
		for(my $j = $i + 1; $j < @strains; $j++) {
			$snps = 0;
			$indels = 0;
			foreach my $chr (keys %{$save{$strains[$i]}}) {
				%seen = ();
				foreach my $pos (keys %{$save{$strains[$i]}{$chr}}) {
					if(exists $save{$strains[$j]}{$chr}{$pos}) {
				#		print $save{$strains[$i]}{$chr}{$pos} . "\t" . $save{$strains[$j]}{$chr}{$pos} . "\n";
						if($save{$strains[$i]}{$chr}{$pos} ne $save{$strains[$j]}{$chr}{$pos}) {
							if(length($save{$strains[$i]}{$chr}{$pos}) != length($save{$strains[$j]}{$chr}{$pos}) && (length($save{$strains[$i]}{$chr}{$pos}) > 1 || length($save{$strains[$j]}{$chr}{$pos}) > 1)) {
								$indels++;
							} else {
								$snps++;
							}
						}
					} else {
					#	print $save{$strains[$i]}{$chr}{$pos} . "\tNA\n";
						if(length($save{$strains[$i]}{$chr}{$pos}) > 1) {
							$indels++;
						} else {
							$snps++;
						}
					}	
					$seen{$pos} = 1;
				#	print $snps . "\t" . $indels . "\n";
				}
				foreach my $pos (keys %{$save{$strains[$j]}{$chr}}) {
					if(exists $seen{$pos}) {
						next;
					}
				#	print "NA\t" . $save{$strains[$j]}{$chr}{$pos} . "\n";
					if(length($save{$strains[$j]}{$chr}{$pos}) > 1) {
						$indels++;
					} else {
						$snps++;
					}
				#	print $snps . "\t" . $indels . "\n";
				}
			}
			print $strains[$i] . " vs " . $strains[$j] . "\t" . $snps . "\t" . $indels . "\n";
		}
	}
}
