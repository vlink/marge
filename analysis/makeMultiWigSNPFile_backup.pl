#!/usr/bin/perl -w

use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';

if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "\t-organism: specify organsim\n";
	print STDERR "\t-strains <list of strains>: Creates output for every strain\n";
	print STDERR "\t-genome <genome>\n";
	print STDERR "\n\nFormats for input file\n";
	print STDERR "\t-homer: homer peak or RNA-Seq file (default)\n";
        print STDERR "\t-format: checks formats.txt to find the user defined format called annotation_format\n";
	print STDERR "\t-hub <hub name>: Name of the hub";
	print STDERR "\t-force: Overwrittes hub if it already exists\n";
	exit;
}

my $bigwig = `which bedGraphToBigWig`;
if($bigwig eq "") {
	print STDERR "\nCould not detect bedGraphToBigWig programm!\n";
	print STDERR "Please install!\n";
	exit;
}

my $config = config::read_config();

$_ = "" for my($organism, $dbh, $sth, $split_sign, $chr_pos, $shift_order, $total_chr, $filename_ref, $filename_strain, $current_chr, $genome, $command, $hub_path, $hub_url, $hub, $wd);
$_ = 0 for my($homer, $format, $count, $min, $force);
$_ = () for my(@strains, @split, @l_ref, @l_mut);
my %order;

GetOptions(     "organism=s" => \$organism,
                "strains=s{,}" => \@strains,
                "homer" => \$homer,
		"genome=s" => \$genome,
                "format" => \$format, 
		"hub=s" => \$hub,
		"force" => \$force)
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

if($genome eq "") {
	print STDERR "Please specify genome!\n";
	exit;
}

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
}

if($hub eq "") {
	print STDERR "Please specify hub name!\n";
	exit;
}

my $color = &generate_colors();

if($config->{'homer_path'} eq "") {
	print STDERR "Here is still some work to do!\n";
	exit;
}

#Time to get all the paths
open HOMER_CONFIG, "<$config->{'homer_path'}/config.txt";
foreach my $line (<HOMER_CONFIG>) {
	chomp $line;
	@split = split('=', $line);
	if($split[0] eq "hubsDir") {
		$hub_path = $split[1];
	}
	if($split[0] eq "hubsUrl") {
		$hub_url = $split[1];
	}
}
close HOMER_CONFIG;
$wd = $hub_path . "/" . $hub;
if(-e $wd) {
	if($force == 0) {
		print STDERR "Hub exists!\n";
		print STDERR "Use -force to overwrite!\n";
		exit;
	} else {
		`rm -rf "$wd"`;
	}
}

`mkdir -p $wd`;
open GENOMES, ">$wd/genomes.txt";
print GENOMES "genome " . $genome . "\n";
print GENOMES "trackDb " . $genome . "/trackDb.txt";
close GENOMES;

my $username = getpwuid( $< );

open HUB, ">$wd/hub.txt";
print HUB "hub $hub\n";
print HUB "shortLabel $hub\n";
print HUB "longLabel $hub\n";
print HUB "genomesFile genomes.txt\n";
print HUB "email " . $username . "\@ucsd.edu\n";
close HUB;

$wd .= "/" . $genome;
`mkdir $wd`;

open TRACK, ">$wd/trackDb.txt";
print TRACK "track $hub\n";
print TRACK "container multiWig\n";
print TRACK "noInherit on\n";
print TRACK "shortLabel $hub\n";
print TRACK "longLabel $hub\n";
print TRACK "type bigWig\n";
print TRACK "configurable on\n";
print TRACK "visibility full\n";
print TRACK "aggregate transparentOverlay\n";
print TRACK "showSubtrackColorOnUi on\n";
print TRACK "autoScale on\n";
print TRACK "windowingFunction average\n";
print TRACK "priority 1.4\n";
print TRACK "alwaysZero on\n";
print TRACK "yLineMark 0\n";
print TRACK "yLineOnOff on\n";		
print TRACK "maxHeightPixels 100:75:11\n\n";


$current_chr = 0;
my $end_strain_old = 0;
my $start_strain_old = 0;

#my $con = config::database_connection();
#$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});


$dbh = DBI->connect("DBI:Pg:dbname=$config->{'database'};host=$config->{'host'}",  $config->{'database_user'}, $config->{'database_pw'}, {'RaiseError' => 1});
 
my $run = 0;

for(my $s = 0; $s < @strains; $s++) {
	print $strains[$s] . "\n";
	$filename_ref = "ref_$strains[$s]";
	open REF, ">$wd/$filename_ref.wig";
	$filename_strain = "strain_$strains[$s]";
	open STRAIN, ">$wd/$filename_strain.wig";
	print STDERR "Collect mutations!\n";
	$sth = $dbh->prepare("SELECT * FROM " . $organism . "_mutations_" . $strains[$s] . " ORDER BY chr, position");
	$sth->execute();
	print STDERR "Start creating file\n";
#	print OUT "track name=\"all-mutations-$strains[$s]\"\n";
	while(my $m = $sth->fetchrow_hashref()) {
		if($m->{'chr'} ne $current_chr) {
			print STDERR "chromosome " . $m->{'chr'} . "\n";
			$current_chr = $m->{'chr'};
			$start_strain_old = 0;
			$end_strain_old = 0;
		}
		print REF "chr" . $m->{'chr'} . "\t" . $m->{'position'} . "\t" . ($m->{'position'} + length($m->{'reference'})) . "\t1\n";
		if($m->{'position'} < $end_strain_old) {
			next;
		}

		if($end_strain_old > $m->{'position'} && ($end_strain_old + 1) < ($m->{'position'} + length($m->{'strain'}))) {
			print STRAIN "chr" . $m->{'chr'} . "\t" . ($end_strain_old + 1) . "\t" . ($m->{'position'} + length($m->{'strain'})) . "\t-1\n";
		} elsif($end_strain_old < $m->{'position'}) {
			print STRAIN "chr" . $m->{'chr'} . "\t" . $m->{'position'} . "\t" . ($m->{'position'} + length($m->{'strain'})) . "\t1\n";
		}
		$end_strain_old = ($m->{'position'} + length($m->{'strain'})); 
		$start_strain_old = $m->{'position'};
	}	
	$sth->finish();
	close REF;
	close STRAIN;
	print STDERR "Convert to bigwig\n";
	$command = "bedGraphToBigWig " . $wd . "/" . $filename_ref . ".wig " . $config->{'homer_path'} . "/data/genomes/" . $genome . "/chrom.sizes " . $wd . "/" . $filename_ref . ".bigWig";
	`$command`;
	$command = "bedGraphToBigWig " . $wd . "/" . $filename_strain . ".wig " . $config->{'homer_path'} . "/data/genomes/" . $genome . "/chrom.sizes " . $wd . "/" . $filename_strain . ".bigWig";
	`$command`;
	print TRACK "track $filename_ref\n";
	print TRACK "bigDataUrl $filename_ref.bigWig\n";
	print TRACK "shortLabel $filename_ref\n";
	print TRACK "longLabel $filename_ref\n";
	print TRACK "type bigWig\n";
	print TRACK "parent $hub\n";
	print TRACK "color $color->[$run]\n\n"; 
	$run++;

	print TRACK "track $filename_strain\n";
	print TRACK "bigDataUrl $filename_strain.bigWig\n";
	print TRACK "shortLabel $filename_strain\n";
	print TRACK "longLabel $filename_strain\n";
	print TRACK "type bigWig\n";
	print TRACK "parent $hub\n";
	print TRACK "color $color->[$run]\n\n"; 
	$run++;
}

`rm $wd/*wig`;
close TRACK;

print STDERR "\tLoad hub: " . $hub_url . "/" . $hub . "/hub.txt\n\n";
sub generate_colors{
	my @color;
	my $run = 0;
	my $add = 255;
	my $steps = 0;
	my $temp = @strains;
	while($temp > 0) {
		$temp = $temp - 3;
		$steps++;
	}
	$add = int($add/$steps);
	for(my $i = 0; $i < @strains; $i++) {
		if($i%3 == 0) {
			#Change first one
			$color[$run] = (int(($i/3) + 1) * $add) . ",0,0";
			$run++;
			$color[$run] = (int(($i/3) + 1) * $add) . ",204,204";
			$run++;
		} elsif($i%3 == 1) {
			$color[$run] = "0," . (int(($i/3) + 1) * $add) . ",0";
			$run++;
			$color[$run] = "204," . (int(($i/3) + 1) * $add) . ",204";
			$run++;
			#Change second one
		} else {
			$color[$run] = "0,0," . (int(($i/3) + 1) * $add);
			$run++;
			$color[$run] = "204,204," . (int(($i/3) + 1) * $add);
			$run++;
		}
	}
	return \@color
}
