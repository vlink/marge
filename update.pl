#/usr/bin/env perl
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use strict;
use warnings;
use Getopt::Long;
use config;
use Data::Dumper;

$_ = () for my(@split, @MARGE_folders, @motif_bg, @motif_analysis_folders, %save_output, @tmp, $list, @STAR_index, @bowtie_index, @complete);
$_ = 0 for my($MARGE_folders_list, $motif_window, $motif_analysis_list, $bowtie_index_list, $STAR_index_list, $motif_bg_list);
$_ = "" for my($tmp_output, $tmp_error, $command, $MARGE_genome, $file);

my $url = "http://homer.ucsd.edu/MARGE/";
 
sub printCMD{
	print STDERR "Script to download already pre-processed data for MARGE\n";
	print STDERR "The script offers downloads for MARGE folders and individualzed genomes, as well as background scans for motifs and pre-processed files for motif analysis per strain (200bp)\n";
	print STDERR "\n\n";
	print STDERR "Usage:\n";
	print STDERR "\t-genome <genome>: The genome for which the motifs were scanned for (e.g. mm10/hg19)\n";
	print STDERR "\t-MARGE_folders <list of folders>: Downloads the list of all folders specified in the list.\n\t\tIf set to all: downloads all available folders\n\t\tIf empty: lists all available folders\n";
	print STDERR "\t-motif_bg <list of motifs>: Downloads all the motif background scans for the motifs listed here.\n\t\tIf set to all: downloads all availabile motifs\n\t\tIf empty: lists all motif backgrounds for genome specified in motif_genome\n";
	print STDERR "\t-motif_analysis_folders <list of folders> : Downloads all the preparsed motif analysis backgrounds specified in the list.\n\t\tIf set to all: downloads all availabe preparsed motif backgrounds\n\t\tIf empty: lists all available preparsed motif backgrounds\n";
	print STDERR "\t-motif_window <bp>: Motif window for which files should be downloaded (default in HOMER motif analysis - 200bp)\n";
	print STDERR "\t-STAR_index <list of individuals/strains>: Downloads the STAR indicies for all individuals specified.\n\t\tIf set to all: downloads all available STAR indicies\n\t\tIf empty: lists all available STAR indicies\n";
	print STDERR "\t-bowtie_index <ist of individuals/strains>: Downloads the bowtie indices for all individuals specified.\n\t\tIf set to all: downloads all availabe bowtie indicies\n\t\tIf empty: lists all available bowtie indicies\n";
	print STDERR "\t-complete <list of individulas/strains>: Downloads all data per strain (MARGE folder, motif background, STAR index and bowtie index>\n\t\tIf set to all: downloads all available folders (NOT RECOMMENDED!)\n\t\tIf empty: Lists all available MARGE folders\n";
	print STDERR "\nList of all content to download:\n";
	print STDERR "\t-MARGE_folders_list: lists all available MARGE folders\n";
	print STDERR "\t-motif_analysis_list: lists all available motif analysis backgrounds\n";
	print STDERR "\t-bowtie_index_list: lists all available bowtie indices\n";
	print STDERR "\t-STAR_index_list: lists all available STAR indicies\n";
	print STDERR "\t-motif_bg_list: lists all available backgrounds for for distance plots\n";
	print STDERR "\n\n\n";
	exit;
}
#CHECK FOR UPDATES

if(@ARGV < 1) {
        &printCMD();
}

GetOptions(	"MARGE_folders:s{,}" => \@MARGE_folders,
		"bowtie_index:s{,}" => \@bowtie_index,
		"STAR_index:s{,}" => \@STAR_index,
		"genome=s" => \$MARGE_genome,
		"motif_bg:s{,}" => \@motif_bg,
		"complete:s{,}" => \@complete,
		"motif_analysis_folders:s{,}" => \@motif_analysis_folders,
		"motif_window=s" => \$motif_window,
		"MARGE_folders_list" => \$MARGE_folders_list,
		"motif_analysis_list" => \$motif_analysis_list,
		"bowtie_index_list" => \$bowtie_index_list,
		"motif_bg_list" => \$motif_bg_list,
		"STAR_index_list" => \$STAR_index_list)
	or die(&printCMD());

my $config = config::read_config();

if(@complete == 1 && $complete[0] eq "") {
	&list_MARGE_folders(1);
	exit;
}

if((@MARGE_folders == 1 && $MARGE_folders[0] eq "") || $MARGE_folders_list == 1) {
#First check for all the lists
	&list_MARGE_folders(1);
	exit;
}

if((@motif_bg == 1 && $motif_bg[0] eq "") || $motif_bg_list == 1) {
	&list_motif_bg(1);
	exit;
}

if((@bowtie_index == 1 && $bowtie_index[0] eq "") || $bowtie_index_list == 1) {
	&list_bowtie_index(1);
	exit;
}

if((@STAR_index == 1 && $STAR_index[0] eq "") || $STAR_index_list == 1) {
	&list_STAR_index(1);
	exit;
}

if((@motif_analysis_folders && $motif_analysis_folders[0] eq "") || $motif_analysis_list == 1) {
	&list_motif_analysis(1);
	exit;
}

if($MARGE_genome eq "") {
	print STDERR "Please specify genome!\n\n";
	exit;
}

if(@complete == 1) {
	@tmp = split(",", $complete[0]);
	if(@tmp > @complete) {
		@complete = @tmp;
	}
}

for(my $i = 0; $i < @complete; $i++) {
	$complete[$i] =~ s/,//g;
	$complete[$i] = uc($complete[$i]);
	$MARGE_folders[$i] = $complete[$i];
	$STAR_index[$i] = $complete[$i];
	$bowtie_index[$i] = lc($complete[$i]);
}

if(@complete >= 1) {
	my $list_marge_folders = &list_MARGE_folders(0);
	my $list_motif_bg = &list_motif_analysis(0);
	my $list_star = &list_STAR_index(0);
	my $list_bowtie = &list_bowtie_index(0);
	if(@complete == 1 && uc($complete[0]) eq "ALL") {
		print STDERR "This is a lot of data to download.\nThis option is not recommended!\nPlease download data one strain/individual at a time\n\n";
		print STDERR "The download will take some time.\n";
		print STDERR "Continuting in 10 seconds\n";
		print STDERR "To abort press Ctrl + C\n";
		print STDERR "Waiting for 10 seconds:\n";
		for(my $i = 0; $i < 10; $i++) {
			print STDERR ".";
			sleep(1);
		}
		print STDERR "\n";
		@MARGE_folders = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @MARGE_folders, $s;
		}
		@motif_analysis_folders = ();
		foreach my $s (keys %{$list->{$MARGE_genome}->{$motif_window}}) {
			push @motif_analysis_folders, $s;
		}
		@STAR_index = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @STAR_index, $s;
		}
		@bowtie_index = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @bowtie_index, $s;
		}
	}
	for(my $i = 0; $i < @MARGE_folders; $i++) {
		if(!exists $list_marge_folders->{$MARGE_genome}->{$MARGE_folders[$i]}) {
			print STDERR "Could not find " . $MARGE_folders[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifying any indicies or with -MARGE_folders_list to see a list of all available indicies\n";
		} else {
			&download_file_MARGE($MARGE_genome, $MARGE_folders[$i]);
		}
	}
	for(my $i = 0; $i < @motif_analysis_folders; $i++) {
		if(!exists $list_motif_bg->{$MARGE_genome}->{$motif_window}->{$motif_analysis_folders[$i]}) {
			print STDERR "Could not find " . $motif_analysis_folders[$i] . " for " . $MARGE_genome . " with window " . $motif_window . "\n";
			print STDERR "Run this script without specifying any motif analysis folders or with parameter -motif_analysis_list to see what is available\n";
		} else {
			&download_motif_analysis($MARGE_genome, $motif_analysis_folders[$i], $motif_window);
		}
	}
	for(my $i = 0; $i < @STAR_index; $i++) {
		if(!exists $list_star->{$MARGE_genome}->{$STAR_index[$i]}) {
			print STDERR "Could not find " . $STAR_index[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifying any indicies or with -STAR_index_list to see a list of all available indicies\n";
		} else {
			&download_STAR_index($MARGE_genome, $STAR_index[$i]);
		}
	}
	for(my $i = 0; $i < @bowtie_index; $i++) {
		if(!exists $list_bowtie->{$MARGE_genome}->{$bowtie_index[$i]}) {
			print STDERR "Could not find " . $bowtie_index[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifying any indicies or with -list_bowtie_index to see a list of all available indicies\n";
		} else {
			&download_bowtie_index($MARGE_genome, $bowtie_index[$i]);
		}
	}
}

if(@MARGE_folders == 1) {
	@tmp = split(",", $MARGE_folders[0]);
	if(@tmp > @MARGE_folders) {
		@MARGE_folders = @tmp;
	}
}

for(my $i = 0; $i < @MARGE_folders; $i++) {
	$MARGE_folders[$i] =~ s/,//g;
	$MARGE_folders[$i] = uc($MARGE_folders[$i]);
}

if(@MARGE_folders >= 1) {
	$list = &list_MARGE_folders(0);
	if(@MARGE_folders == 1 && uc($MARGE_folders[0]) eq "ALL") {
		@MARGE_folders = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @MARGE_folders, $s;
		}
	}
	for(my $i = 0; $i < @MARGE_folders; $i++) {
		if(!exists $list->{$MARGE_genome}->{$MARGE_folders[$i]}) {
			print STDERR "Could not find " . $MARGE_folders[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifying any indicies or with -MARGE_folders_list to see a list of all available indicies\n";
		} else {
			&download_file_MARGE($MARGE_genome, $MARGE_folders[$i]);
		}
	}
}

if(@motif_analysis_folders == 1) {
	@tmp = split(",", $motif_analysis_folders[0]);
	if(@tmp > @motif_analysis_folders) {
		@motif_analysis_folders = @tmp;
	}
}

for(my $i = 0; $i < @motif_analysis_folders; $i++) {
	$motif_analysis_folders[$i] =~ s/,//g;
	$motif_analysis_folders[$i] = uc($motif_analysis_folders[$i]);
}

if(@motif_analysis_folders >= 1) {
	if($motif_window == 0) {
		print STDERR "Motif window was not specified!\n";
		print STDERR "Set it to default - 200bp\n";
		$motif_window = 200;
	} 
	$motif_window =~ s/bp//g;
	$motif_window =~ s/BP//g;
	$motif_window = $motif_window . "bp";

	$list = &list_motif_analysis(0);
	if(@motif_analysis_folders == 1 && uc($motif_analysis_folders[0]) eq "ALL") {
		@motif_analysis_folders = ();
		foreach my $s (keys %{$list->{$MARGE_genome}->{$motif_window}}) {
			push @motif_analysis_folders, $s;
		}
	}
	for(my $i = 0; $i < @motif_analysis_folders; $i++) {
		if(!exists $list->{$MARGE_genome}->{$motif_window}->{$motif_analysis_folders[$i]}) {
			print STDERR "Could not find " . $motif_analysis_folders[$i] . " for " . $MARGE_genome . " with window " . $motif_window . "\n";
			print STDERR "Run this script without specifying any motif analysis folders or with parameter -motif_analysis_list to see what is available\n";  
		} else {
			&download_motif_analysis($MARGE_genome, $motif_analysis_folders[$i], $motif_window);
		}
	}
}

if(@STAR_index == 1) {
	@tmp = split(",", $STAR_index[0]);
	if(@tmp > @STAR_index) {
		@STAR_index = @tmp;
	}
}

for(my $i = 0; $i < @STAR_index; $i++) {
	$STAR_index[$i] =~ s/,//g;
	$STAR_index[$i] = lc($STAR_index[$i]);
}

if(@STAR_index >= 1) {
	$list = &list_STAR_index(0);
	if(@STAR_index == 1 && $STAR_index[0] eq "ALL") {
		@STAR_index = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @STAR_index, $s;
		}
	}
	for(my $i = 0; $i < @STAR_index; $i++) {
		if(!exists $list->{$MARGE_genome}->{$STAR_index[$i]}) {
			print STDERR "Could not find " . $STAR_index[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifying any indicies or with -STAR_index_list to see a list of all available indicies\n";
		} else {
			&download_STAR_index($MARGE_genome, $STAR_index[$i]);
		}
	}	
}

if(@bowtie_index == 1) {
	@tmp = split(",", $bowtie_index[0]);
	if(@tmp > @bowtie_index) {
		@bowtie_index = @tmp;
	}
}

for(my $i = 0; $i < @bowtie_index; $i++) {
	$bowtie_index[$i] =~ s/,//g;
	$bowtie_index[$i] = lc($bowtie_index[$i]);
}

if(@bowtie_index >= 1) {
	$list = &list_bowtie_index(0);
	if(@bowtie_index == 1 && $bowtie_index[0] eq "ALL") {
		@bowtie_index = ();
		foreach my $s (keys %{$list->{$MARGE_genome}}) {
			push @bowtie_index, $s;
		}
	}
	for(my $i = 0; $i < @bowtie_index; $i++) {
		if(!exists $list->{$MARGE_genome}->{$bowtie_index[$i]}) {
			print STDERR "Could not find " . $bowtie_index[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifying any indicies or with -list_bowtie_index to see a list of all available indicies\n";
		} else {
			&download_bowtie_index($MARGE_genome, $bowtie_index[$i]);
		}
	}	
}



if(@motif_bg == 1) {
	@tmp = split(",", $motif_bg[0]);
	if(@tmp > @motif_bg) {
		@motif_bg = @tmp;
	}
}

for(my $i = 0; $i < @motif_bg; $i++) {
	$motif_bg[$i] =~ s/,//g;
	$motif_bg[$i] = $motif_bg[$i];
}

if(@motif_bg >= 1) {
	$list = &list_motif_bg(0);
	if(@motif_bg == 1 && $motif_bg[0] eq "ALL") {
		@motif_bg = ();
		foreach my $m (keys %{$list->{$MARGE_genome}}) {
			push @motif_bg, $m;
		}
	}
	for(my $i = 0; $i < @motif_bg; $i++) {
		if(!exists $list->{$MARGE_genome}->{$motif_bg[$i]}) {
			print STDERR "Could not find " . $motif_bg[$i] . " for " . $MARGE_genome . "\n";
			print STDERR "Run this script without specifiying any motif backgrounds or with -list_motif_bg to see a list of all available motif backgorunds\n";
		} else {
			&download_motif_bg($MARGE_genome, $motif_bg[$i]);
		}
	}
}

sub download_file_MARGE{
	my $g = $_[0];
	my $f = $_[1];
	print STDERR "Downloading " . $f . "\n";
	$command = "wget " . $url . "/MARGE_folders/" . $g . "_" . $f . ".tar.gz";
	`$command`;
	#Unzip file
	print STDERR "Extracting and moving " . $f . "\n";
	$command = "tar xfvz " . $g . "_" . $f . ".tar.gz";
	`$command`;
	$command = "mv " . uc($f) . " " . $config->{'data_folder'};
	#Move file
	`$command`;
	$command = "rm " . $g . "_" . $f . ".tar.gz";
	if(-e $g . "_" . $f . ".tar.gz") {
		`$command`;
	}
}

sub download_motif_bg{
	my $g = $_[0];
	my $f = $_[1];
	print STDERR "Downloading " . $f . "\n";
	$command = "wget " . $url . "/motif_bg/" . $g . "_" . $f . ".txt.tar.gz";
	`$command`;
	#Unzip file
	print STDERR "Extracting and moving " . $f . "\n";
	$command = "tar xfvz " . $g . "_" . $f . ".txt.tar.gz";
	`$command`;
	$command = "mv " . $g . "_" . $f . ".txt " . $config->{'motif_bg'};
	`$command`;
	$command = "rm " . $g . "_" . $f . ".txt.tar.gz";
	if(-e $g . "_" . $f . ".txt.tar.gz") {
		`$command`;	
	}
}

sub download_STAR_index{
	my $g = $_[0];
	my $f = $_[1];
	print STDERR "Downloading " . $f . "\n";
	$command = "wget " . $url . "/STAR_index/" . $g . "_" . lc($f) . ".tar.gz";
	`$command`;
	#Unzip file
	print STDERR "Extracting file " . $f . "\n";
	$command = "tar xvfz " . $g . "_" . $f . ".tar.gz";
	`$command`;
	$command = "rm " . $g . "_" . $f . ".tar.gz";
	if(-e $g . "_" . $f . ".tar.gz") {
		`$command`;
	}
	print STDERR "Index " . $f . " successfully downloaded and extracted - please move to your STAR index directory\n";
}


sub download_bowtie_index{
	my $g = $_[0];
	my $f = $_[1];
	print STDERR "Downloading " . $f . "\n";
	$command = "wget " . $url . "/bowtie_index/" . $g . "_" . lc($f) . "_index.tar.gz";
	`$command`;
	#Unzip file
	print STDERR "Extracting file " . $f . "\n";
	$command = "tar xvfz " . $g . "_" . $f . "_index.tar.gz";
	`$command`;
	$command = "rm " . $g . "_" . $f . "_index.tar.gz";
	if(-e $g . "_" . $f . "_index.tar.gz") {
		`$command`;
	}
	print STDERR "Index " . $f . " successfully downloaded and extracted - please move to your bowtie2 index directory\n";
}

sub download_motif_analysis{
	my $g = $_[0];
	my $f = $_[1];
	my $w = $_[2];
	if(!-e $config->{'data_folder'} . "/" . $f) {
		print STDERR "There is no data folder for $f!\n";
		print STDERR "Please download first the MARGE folder for this individual genome!\n";
		exit;
	}

	print STDERR "Downloading background for findMotifs for " . $f . " for " . $w . "\n";
	$command = "wget " . $url . "/motif_analysis/" . $g . "_preparsed_" . $f . "_" . $w . ".tar.gz";
	print $command . "\n";
	`$command`; 
	#Extract
	print STDERR "Extracting and moving " . $f . "\n";
	$command = "tar xfvz " . $g . "_preparsed_" . $f . "_" . $w . ".tar.gz";
	print $command . "\n";
	`$command`;
	if(!-e $config->{'data_folder'} . "/" . $f . "/preparsed") {
		print STDERR "Generating folder for preparsed data!\n";
		$command = "mkdir " . $config->{'data_folder'} . "/" . $f . "/preparsed";
		`$command`;
	}
	$command = "cp " . $g . "_preparsed_" . $f . "_" . $w . "/* " . $config->{'data_folder'} . "/" . $f . "/preparsed/";	
	print $command . "\n";
	`$command`;
	$command = "rm -rf " . $g . "_preparsed_" . $f . "_" . $w;
	if(-e $g . "_preparsed_" . $f . "_" . $w) {
		print $command . "\n";
		`$command`;
	}
	if(-e $g . "_preparsed_" . $f . "_" . $w . ".tar.gz") {
		$command = "rm " . $g . "_preparsed_" . $f . "_" . $w . ".tar.gz";
		print $command . "\n";
		`$command`;
	}
}

sub list_MARGE_folders {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/MARGE_folders/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", substr($line, 0, length($line) - 7));
		$save_output{$split[0]}{$split[1]} = 1;
	}
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print STDERR "Folders for genome " . $genome . "\n";
			foreach my $strain (sort {$a cmp $b} keys %{$save_output{$genome}}) {
				print STDERR "\t" . $strain . "\n";
			}
		}
	}
	return \%save_output;
}

sub list_motif_bg {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/motif_bg/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", substr($line, 0, length($line) - 11));
		$save_output{$split[0]}{$split[1]} = 1;
	}
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print STDERR "Motifs that were scanned for genome " . $genome . "\n";
			foreach my $strain (sort {$a cmp $b} keys %{$save_output{$genome}}) {
				print STDERR "\t" . $strain . "\n";
			}
		}
	}
	return \%save_output;
}

sub list_motif_analysis {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/motif_analysis/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", substr($line, 0, length($line) - 7));
		$save_output{$split[0]}{$split[3]}{$split[2]} = 1;
	}
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print STDERR "Motifs that were scanned for genome " . $genome . "\n";
			foreach my $window (keys %{$save_output{$genome}}) {
				print STDERR "Window that was used for scan: " . $window . "\n";
				foreach my $strain (sort {$a cmp $b} keys %{$save_output{$genome}{$window}}) {
					print STDERR "\t\t" . $strain . "\n";
				}
			}
		}	
	}
	return \%save_output;
}

sub list_bowtie_index {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	`wget -O $tmp_output $url/bowtie_index/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", $line);
		$save_output{$split[0]}{$split[1]} = 1;
	}
	close FH;
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print STDERR "Bowtie indicies for individuals based on " . $genome . " genome\n";
			foreach my $s (sort {$a cmp $b} keys %{$save_output{$genome}}) {
				print STDERR "\t\t" . $s . "\n";
			}
		}
	}
	return \%save_output;
}

sub list_STAR_index {
	my $print = $_[0];
	$tmp_output = "output_" . rand(5) . ".txt";
	$tmp_error = "error_" . rand(5) . ".txt";
	my @a;
	`wget -O $tmp_output $url/STAR_index/list.txt 2> $tmp_error`;
	open FH, "<$tmp_output";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split("_", $line);
		@a = split('\.', $split[1]);
		$save_output{$split[0]}{$a[0]} = 1;
	}
	close FH;
	`rm $tmp_output`;
	`rm $tmp_error`;
	if($print == 1) {
		foreach my $genome (keys %save_output) {
			print STDERR "STAR indicies for individuals based on " . $genome . " genome\n";
			foreach my $s (sort {$a cmp $b} keys %{$save_output{$genome}}) {
				print STDERR "\t\t" . $s . "\n";
			}
		}
	}
	return \%save_output;
}

