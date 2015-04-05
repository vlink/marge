#!/usr/bin/perl -w

use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';

if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "\t-organism: specify organsim\n";
	print STDERR "\t-strains <list of strains>: Creates output for every strain\n";
	print STDERR "\n\nFormats for input file\n";
	print STDERR "\t-homer: homer peak or RNA-Seq file (default)\n";
        print STDERR "\t-format: checks formats.txt to find the user defined format called annotation_format\n";
	exit;
}


$_ = "" for my($organism, $dbh, $sth, $split_sign, $chr_pos, $shift_order, $total_chr, $file_name, $current_chr);
$_ = 0 for my($homer, $format, $count, $min);
$_ = () for my(@strains, @split, @l_ref, @l_mut);
my %order;

GetOptions(     "organism=s" => \$organism,
                "strains=s{,}" => \@strains,
                "homer" => \$homer,
                "format" => \$format)
        or die ("Error in command line arguments!\n");

if($organism eq "") {
	print STDERR "Please specify organsim!\n";
	exit;
}

if($homer == 0 && $format == 0) {
	$homer = 1;
}

if(@strains == 0) {
        print STDERR "Please specify strains (-strains)!\n";
        exit;
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
}

$current_chr = 0;
my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
for(my $s = 0; $s < @strains; $s++) {
	print $strains[$s] . "\n";
	$file_name = "ucsc_mutation_file_$strains[$s].txt";
	open OUT, ">$file_name";
	print STDERR "Collect mutations!\n";
	$sth = $dbh->prepare("SELECT * FROM " . $organism . "_mutations_" . $strains[$s] . " ORDER BY chr, position");
	$sth->execute();
	print STDERR "Start creating file\n";
	print OUT "track name=\"all-mutations-$strains[$s]\"\n";
	while(my $m = $sth->fetchrow_hashref()) {
		print OUT "chr" . $m->{'chr'} . "\t" . ($m->{'position'} - 1) . "\t" . ($m->{'position'} - 1 + length($m->{'reference'})) . "\t" . $m->{'reference'} . "->" . $m->{'strain'} . "\t1\t+\n";
		if($m->{'chr'} ne $current_chr) {
			print STDERR "chromosome $m->{'chr'}\n";
			$current_chr = $m->{'chr'};
		}
	}	
	$sth->finish();
	close OUT;
	if(-e $file_name . ".gz") {
		`rm $file_name.gz`;
	}
	`gzip $file_name`;
}
