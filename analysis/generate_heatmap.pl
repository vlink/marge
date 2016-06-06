#!/usr/bin/perl
use strict;
use Getopt::Long;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
$_ = () for my(@files, @names);
$_ = 0 for my($sig, $threshold);
$_ = "" for my ($output);

sub printCMD {
	print STDERR "\t-files <comma separated list of files>\n";
	print STDERR "\t-names <comma separated list of names>\n";
	print STDERR "\t-sig: only significant motifs (threshold p-value 0.001)\n";
	print STDERR "\t-threshold: threshold for significant\n";
	print STDERR "\t-output: output name for R file and pdf (default summary_heatmap.R)\n";
	exit;
}

if(@ARGV < 1) {
        &printCMD();
}

#Check mandatory command line arguments
my %mandatory = ('-files' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "files=s{,}" => \@files,
		"names=s{,}" => \@names,
		"sig" => \$sig,
		"threshold=s" => \$threshold,
		"output=s" => \$output)
        or die("Error in command line options!\n");
#First step: Get the sequences for the peaks
for(my $i = 0; $i < @files; $i++) {
	$files[$i] =~ s/,//g;
	if(defined $names[$i]) {
		$names[$i] =~ s/,//g;
	} else {
		$names[$i] = $files[$i];
	}
}

if($output eq "") {
	$output = "summary_heatmap.R";
}

if($threshold == 0) {
	$threshold = 0.001
}

my @split;
my $strain1;
my $strain2;
my $tmnt;
my $AB;
my $current_motif;
my @file_name;

open R, ">tmp.R";
my $f;
my $length_one = 0;
my $length_two = 0;

for(my $i = 0; $i < @files; $i++) {
	$f = $files[$i];
	chomp $f;
	@file_name = split("/", $f);
	@split = split("_", $file_name[-1]);
#	if($split[-1] ne "KLA.R" && $split[-1] ne "notx.R") {
#		next;
#	}
#	$strain1 = $split[1];
#	$strain2 = $split[2];
#	$tmnt = substr($split[4], 0, length($split[4]) - 2);
#	$AB = $split[3];
	open FH, "<$f";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 4) eq "plot") {
			@split = split('\\\n', $line);
			$current_motif = substr($split[-1], 0, length($split[-1]) - 2);
		}		
		if(substr($line, 0, 3) eq "one") {
			print R $line . "\n";
			$length_one = length($line);
		}
		if(substr($line, 0,3 ) eq "two") {
			print R $line . "\n";
			$length_two = length($line);
		}
		if(substr($line, 0, 6) eq "p_both" && $length_one > 10 && $length_two > 10) {
			print R $line . "\n";
			print R "print(paste(\"" . $names[$i] . "\", \":" . $current_motif . " \", p_both, sep=\"\"))\n";  
		}
	}	
	close FH;
}

`Rscript tmp.R > tmp_output_rfile.txt`;

open FH, "<tmp_output_rfile.txt";

my @name;
my %save;
my %comp;

my $print_line;
my $sig = 0;

foreach my $line (<FH>) {
	chomp $line;
	@split = split(" ", substr($line, 5, length($line) - 6));
	@name = split(":", $split[0]);
	$save{$name[1]}{$name[0]} = $split[1];
#	$comp{$name[0] . "_" . $name[1] . "_" . $name[2]} = 1;
}
close FH;
open OUT, ">matrix_" . $output;
my $print;

my $first = 0;
foreach my $motif (keys %save) {
	if($first == 0) {
		print OUT "TF";
		for(my $i = 0; $i < @names; $i++) {
#		foreach my $comp (keys %comp) {
			print OUT "\t" . $names[$i];
		}
		print OUT "\n";
		$first++;
	}
	$print_line = $motif;
	$print = 0;
	for(my $i = 0; $i < @names; $i++) {
#	foreach my $comp (keys %comp) {
		if(exists $save{$motif}{$names[$i]}) {
			if($save{$motif}{$names[$1]} < 1.0e-40) {
				$save{$motif}{$names[$i]} = 1.0e-40;
			}
			$print_line .= "\t" . $save{$motif}{$names[$i]};
			if($save{$motif}{$names[$i]} < $threshold) {
				$print = 1;
			}
		} else {
			$print_line .= "\t1";
		}
	}
	if($sig == 1 && $print == 1) {
		print OUT $print_line . "\n";
	} elsif($sig == 0) {
		print OUT $print_line . "\n";
	}
}
close OUT;

open OUT, ">$output";
print OUT "library(\"gplots\")\n";
print OUT "matrix <- as.matrix(read.delim(\"matrix_" . $output . "\", row.names=1))\n";
print OUT "pdf(\"" . substr($output, 0, length($output) - 2) . ".pdf\", width=5, height=10)\n";
print OUT "custom.color.fun = colorRampPalette(c(\"blue\", \"white\", \"red\"), bias = 1, space = \"rgb\")\n";
print OUT "heatmap.2(matrix, main=\"Summary p-values \", trace = \"n\", cexCol = 0.4, scale = \"n\", col = custom.color.fun(100), Colv=NA, Rowv=NA, dendrogram='n', cex=0.5, margin=c(8,12), cexRow=0.4)\n";
print OUT "dev.off()\n"; 
close OUT;
`Rscript $output`;
