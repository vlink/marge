#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if(@ARGV < 3) {
	print STDERR "Please insert input file, start and end\n";
	exit;
}

my $input = $ARGV[0];
my $start = $ARGV[1];
my $end = $ARGV[2];

my $matrix_start = 0;
my $matrix_end = 0;

my $header = `head -n1 $input`;
my @split = split('\t', $header);

my @name;

for(my $i = 1; $i < @split; $i++) {
	@name = split("-", $split[$i]);
	if($matrix_start == 0 && $name[2] > $start) {
		$matrix_start = $i;
	}
	if($matrix_end == 0 && $name[2] > $end) {
		$matrix_end = $i - 1;
	}	
}

open my $fh, "<", $input;
my $n = "part_" . substr($input, 0, length($input) - 4) . "_" . $start . "_" . $end . ".txt";
open OUT, ">$n";

my $count = 1;
for(my $i = $matrix_start; $i < $matrix_end; $i++) {
	@name = split("-", $split[$i]);
	print OUT "\t" . $name[2];
}
print OUT "\n";

while(my $line = <$fh>) {
	if($count > $matrix_start && $count <= $matrix_end) {
		@split = split('\t', $line);
		@name = split("-", $split[0]);
		print OUT $name[2];
		#	print OUT $split[0];
		for(my $i = $matrix_start; $i < $matrix_end; $i++) {
			print OUT "\t" . $split[$i];
		}
		print OUT "\n";
	}
	$count++;
}
close $fh;
close OUT;
