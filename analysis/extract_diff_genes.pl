#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
require '../general/config.pm';

sub usage() {
	print STDERR "Usage:\n";
	print STDERR "Run analyzeRepeats with -noadj and afterwards getDiffExpression.pl before you use this programm!\n";
	print STDERR "\n\n";
	print STDERR "$0 -file <diff expressed file> <options>\n";
	print STDERR "\t-FDR <value>: Use FDR additionally to filter (Only useful with replicates!)\n";
	print STDERR "\t-FDR-only: Use only FDR for filtering (default FDR: 0.01)\n";
	print STDERR "\t-p <value>: Set p-value for filtering (default: 0.001)\n";
	print STDERR "\t-both: Uses p-value and FDR to filter (default: only p-value)\n";
	print STDERR "\t-top <value>: Number of most differenlty expressed genes that are reported (default: 100)\n";
	print STDERR "\t-output <filename>: Output file\n";
	print STDERR "\t-seperate: Outputs all top candidates for each comparison in an extra file\n";
}

if(@ARGV < 1) {
	&usage();
	exit;
}

$_ = "" for my($file, $output);
my $fdr = 0.01;
my $pvalue = 0.001;
my $top = 100;
$_ = 0 for my($fdr_only, $sample, $f, $comp, $div, $rest, $count, $both, $min_pvalue, $min_name, $seperate);
$_ = () for my(%comparisons, @candidates, %num, %cand, $name);

my %mandatory = ('-file' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "file=s" => \$file,
                "output=s" => \$output,
		"FDR=s" => \$fdr,
		"FDR_only" => \$fdr_only,
		"p=s" => \$pvalue,
		"both" => \$both,
		"top=s" => \$top, 
		"seperate" => \$seperate)
        or die ("Error in command line arguments!\n");

my($filename, $dir, $suffix) = fileparse($file);
open FH, "<$file";
my @split;

if($fdr != 0.01 && $both == 0 && $fdr_only == 0) {
	print STDERR "FDR is defined but no comparions methods!\n";
	print STDERR "Using FDR and p-value for filtering!\n";
	$both = 1;
}

foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if($f == 0) {
		for(my $i = 8; $i < @split; $i++) {
			if($split[$i] =~ /vs\..*logFC/) {
				$comp = (@split - $i)/4;
				print STDERR "" . $comp . " comparisons in this file!\n";
				$rest = $comp;
				$sample = 1;
				while($rest > $sample) {
					$sample++;
					$rest = $rest/$sample;
				}
				$sample++;
				print STDERR "" . $sample . " samples were used in this file!\n";
				last;
			}
		}
		$f++;
		$count = 0;
		for(my $i = 8 + $sample; $i < @split; $i = $i + 4) {
			$name = substr($split[$i], 0, length($split[$i]) - 6);
			$name =~ s/\.//g;
			$name =~ s/\s/_/g;
			$comparisons{$count} = $name;
			$count++;
		}
		next;	
	}
	if($both == 0 && $fdr_only == 0) {
		$count = 0;
		for(my $i = 8 + $sample; $i < @split; $i = $i + 4) {
			if($split[$i + 2] <= $pvalue) {
				push(@{$candidates[$count]}, $line);		
				$cand{$count}{$split[$i+2]}{$split[0]} = $line;
			}
			$count++;	
		}
	} elsif($fdr_only == 1) {
		$count = 0;
		for(my $i = 8 + $sample; $i < @split; $i = $i + 4) {
			if($split[$i + 3] <= $fdr) {
				push(@{$candidates[$count]}, $line);
				$cand{$count}{$split[$i+2]}{$split[0]} = $line;
			}
			$count++;
		}
	} else {
		$count = 0;
		for(my $i = 8 + $sample; $i < @split; $i = $i + 4) {
			if($split[$i+2] <= $pvalue && $split[$i+3] <= $fdr) {
				push(@{$candidates[$count]}, $line);
				$cand{$count}{$split[$i+2]}{$split[0]} = $line;
			}
		}
	}
}

print STDERR "Filter out top " . $top . " candidates per comparison!\n";

my %sort;
my $print;
my %to_file;
foreach my $c (keys %cand) {
	$print = 0;
	if($seperate == 1) {
		$name = $dir . "top_" . $top . "_" . $comparisons{$c} . ".txt";
		open OUT, ">$name";
	}
	foreach my $p (sort {$a <=> $b } keys %{$cand{$c}}) {
		foreach my $tmp (keys $cand{$c}{$p}) {
			if($print < $top) {
				$to_file{$tmp} = $cand{$c}{$p}{$tmp};
				if($seperate == 1) {
					print OUT $cand{$c}{$p}{$tmp} . "\n";
				}
				$print++;
			} else { 
				if($seperate == 1) { close OUT; }
				last; 
			}
		}
	}
}
if($seperate == 0) {
	$name = $dir . "top_" . $top . "_" . $filename;
	open OUT, ">$name";
	foreach my $key (keys %to_file) {
		print OUT $to_file{$key} . "\n";
	}
	close OUT;
}
