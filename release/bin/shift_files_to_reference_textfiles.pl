#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Storable;
use config;
use general;
use Set::IntervalTree;
use Data::Dumper;

$_ = 0 for my ($sam, $peak, $tag, $help, $print_header, $allele, $bed, $bismark);
$_ = "" for my ($dir, $chr, $last, $data_dir, $out_name, $last_strain);
$_ = () for my (@files, @strains, @split, %last, %lookup, %tree, @tmp, $fetch, $pos_shifted);

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
	print STDERR "\t-dir <directory with files to shift>: has to be the same strain\n";
	print STDERR "\t-files <list with files>: comma separated list\n";
	print STDERR "\t-ind: one or several individuals - comma separated\n";
	print STDERR "\tIf only one strain is specified all files are shifted from this strain\n";
	print STDERR "\tIf several individuals are specified it has to match the number of files specified!\n";
	print STDERR "\t-data_dir <directory>: directory with data: default specified in config\n";
	print STDERR "\n\nFormat to shift (default = sam):\n";
	print STDERR "\t-sam: shifts sam file (result file from mapping with bowtie or STAR\n";
	print STDERR "\t-homer: shifts HOMER peak and RNA files\n";
	print STDERR "\t-tag: shifts tag directory\n";
	print STDERR "\t-bed: shifts bed file\n";
	print STDERR "\t-bismark: shifts cytosine report output form BISMARK\n";
	print STDERR "\n\n";
	print STDERR "\t-h | --help: prints help\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my %mandatory = ('-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "dir=s" => \$dir,
                "files=s{,}" => \@files,
                "ind=s{,}" => \@strains,
                "data_dir=s" => \$data_dir,
                "sam" => \$sam,
                "homer" => \$peak,
                "tag" => \$tag,
		"bed" => \$bed,
		"bismark" => \$bismark,
                "h" => \$help,
                "help" => \$help)
        or die (&printCMD());

#Define variables
if($help == 1) {
        &printCMD();
}
if($sam == 0 && $peak == 0 && $tag == 0 && $bed == 0 && $bismark == 0) {
	$sam = 1;
}
if(@files == 0 && $dir eq "") {
	print STDERR "Files for shifting are missing!\n";
	exit;
}
#Read in files from directory
if(@files == 0) {
	if($sam == 1) {
		@files = `ls $dir/*sam`;
	}
	if($tag == 1) {
		@files = `ls $dir/* -d`;
	}
	if($peak == 1) {
		@files = `ls $dir/*`;
	}
	if($bed == 1) {
		@files = `ls $dir/*bed`;
	}
	if($bismark == 1) {
		print STDERR "Please specify files!\n";
		exit;
	}
}
#Check if files and strains have the same length
if(@strains > 1 && @strains != @files) {
	print STDERR "Error: Length of file list and strains list is different\n";
	exit;
}
#Remove commas
@tmp = split('\s', $files[0]);
if(@tmp > @files) {
	@files = @tmp;
}
for(my $i = 0; $i < @files; $i++) {
	chomp $files[$i];
	$files[$i] =~ s/,//g;
}
for(my $i = 0; $i < @strains; $i++) {
	chomp $strains[$i];
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}
if(@strains == 1 && @files > 1) {
	@strains = (@strains, ($strains[0]) x (@files - 1));
}
if($data_dir eq "") {
	$data_dir = config::read_config()->{'data_folder'};
}
#Save all vectors befor shifting!
#There is only one strain specified - all files are shifted with the same vector
print STDERR "Save shifting vector\n";
for(my $i = 0; $i < @strains; $i++) {
	if($last_strain ne $strains[$i]) {
		my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data_dir, "strain_to_ref");
		%tree = %{$tree_ref};
		%last = %{$last};
		%lookup = %{$lookup};
		print STDERR "\tsuccessful for $strains[$i]\n";
	}
	if($sam == 1) {
		$out_name = substr($files[$i], 0, length($files[$i]) - 4) . "_shifted_from_" . $strains[$i] . ".sam";
		&shift_sam_file($files[$i], $strains[$i], $out_name);
	} elsif($peak == 1) {
		$out_name = $files[$i] . "_shifted_from_" . $strains[$i] . ".txt";
		&shift_peak_file($files[$i], $strains[$i], $out_name);
	} elsif($bed == 1) {
		$out_name = substr($files[$i], 0, length($files[$i]) - 4) . "_shifted_from_" . $strains[$i] . ".txt";
		&shift_bed_file($files[$i], $strains[$i], $out_name);
	} elsif($bismark == 1) {
		$out_name = substr($files[$i], 0, length($files[$i]) - 4) . "_shifted_from_" . $strains[$i] . ".txt";
		&shift_bismark($files[$i], $strains[$i], $out_name);
	} else {
		$out_name = $files[$i];
		while(substr($out_name, length($out_name) - 1) eq "/") {
			$out_name = substr($out_name, 0, length($out_name) - 1);
		}
		$out_name .= "_shifted_from_" . $strains[$i];
		&shift_tag_directory($files[$i], $out_name);
	}
	$last_strain = $strains[$i];
}

#Shifts position
sub shift{
	my $chr = $_[0];
	my $chr_num = $chr;
	if(exists $lookup{$chr}) {
		$chr_num = $lookup{$chr};
	}
	my $allele = $_[1];
	my $pos_to_shift = $_[2];
	if($pos_to_shift > $last{$chr}{$allele}{'pos'}) {
		$pos_shifted = $pos_to_shift + $last{$chr}{$allele}{'shift'};
	} else {
		$fetch = $tree{$chr_num}->{$allele}->fetch($pos_to_shift, $pos_to_shift + 1);
		if(scalar(@$fetch) == 0) {
			$pos_shifted = $pos_to_shift;
		} else {
			$pos_shifted = $pos_to_shift + $fetch->[0]->{'shift'};
		}
	}
	return $pos_shifted;
}



sub shift_sam_file{
	#Shift files when all or one chromosomes are saved in memory and sam file is not
	my $file = $_[0];
	my $strain = $_[1];
	my $out_file = $_[2];
#	open OUT, ">$out_file" or die "Can't open $out_file: $!\n";
	open my $out, ">", $out_file or die "Can't open $out_file: $!\n";
#	open FH, "<$file";
	open my $fh, "<", $file or die "Can't open $file: $!\n";
	print STDERR "Reading in " . $file . "\n";
	while (my $line = <$fh>) {
#	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@tmp = split("_", $split[2]);
		#print sam file header
		if(substr($line, 0, 1) eq "@" && $print_header == 0) {
		#	print OUT $line . "\n";
			print $out $line . "\n";
			next;
		}
		$print_header++;
		#print all sequences that could not be mapped
		if(length($split[2]) < 4) {
		#	print OUT $line. "\n";
			print $out $line. "\n";
			next;
		}
		$chr = substr($tmp[0], 3);
		if(@tmp > 1) {
			$allele = $tmp[2];
		} else {
			$allele = 1;
		}
		#chromosomes without mutations can not be shifted		
		if(!exists $last{$chr}{$allele}) {
		#	print OUT $line . "\n";
			print $out $line . "\n";
			next;
		}
	#	print OUT $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . &shift($chr, $allele, $split[3]);
		print $out $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . &shift($chr, $allele, $split[3]);
		for(my $j = 4; $j < @split; $j++) {
		#	print OUT "\t" . $split[$j];
			print $out "\t" . $split[$j];
		}
	#	print OUT "\n";
		print $out "\n";
	}
#	close OUT;
	close $out;
	close $fh;
}

sub shift_peak_file{
	my $file = $_[0];
	my $strain = $_[1];
	my $out_file = $_[2];
	my $start;
	my $end;
	open OUT, ">$out_file" or die "Can't open $out_file: $!\n";
	open FH, "<$file";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@tmp = split("_", $split[1]);
		$chr = substr($tmp[0], 3);
		if(@tmp > 1) {
			$allele = $tmp[2];
		} else {
			$allele = 1;
		}
		#chromosomes without mutations can not be shifted
		if(!exists $last{$chr}{$allele}) {
			print OUT $line . "\n";
			next;
		}
		$start = &shift($chr, $allele, $split[2]);
		$end = &shift($chr, $allele, $split[3]);
		if($end < $start) { $end = $start + 1; }
		print OUT $split[0] . "\t" . $split[1] . "\t" . $start . "\t" . $end;
		for(my $j = 4; $j < @split; $j++) {
			print OUT "\t" . $split[$j];
		}
		print OUT "\n";
	}
	close OUT;
}

sub shift_tag_directory{
	my $directory = $_[0];
	print STDERR "shifting " . $directory . "\n";
	my @files = `ls $directory/*tags.tsv`;
	my $output_dir = $_[1];
	#Shifted tag directory exists already
	if(-e $output_dir) {
		print STDERR "" . $output_dir . " already exists!\n";
		print STDERR "Overwriting directory in 10 seconds\n";
		print STDERR "Press Ctrl + C to stop\n";
		print STDERR "Waiting: ";
		for(my $i = 0; $i < 10; $i++) {
			print STDERR ".";
			sleep(1)
		}
		print STDERR "\n";
		`rm -rf $output_dir/*`;
	} else {
		my $error = "tmp" . rand(15);
		`mkdir $output_dir 2> $error`;
		my $lll = `wc -l $error`;
		if((split('\s+', $lll))[0] > 0) {
			`cat $error`;
			`rm $error`;
			exit;
		}
		`rm $error`;
	}
	#Copy txt files
	`cp $directory/*txt $output_dir`;
	foreach my $file (@files) {
		chomp $file;
		@split = split("/", $file);
		$file = $split[-1];
		open OUT, ">$output_dir/$file";
		open FH, "<$directory/$file";
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			@tmp = split("_", $split[1]);
			$chr = substr($tmp[0], 3);
			if(@tmp > 1) {
				$allele = $tmp[2];
			} else {
				$allele = 1;
			}
			print OUT $split[0] . "\t" . $split[1] . "\t" . &shift($chr, $allele, $split[2]) . "\t" . $split[3] . "\t" . $split[4] . "\t" . $split[5] . "\n";
		}
		close OUT;
		close FH;
	}
}

sub shift_bed_file{
	my $file = $_[0];
	my $strain = $_[1];
	my $out_file = $_[2];
	print STDERR "shifting " . $file . "\n";
	open OUT, ">$out_file" or die "Can't open file!\n";
	open FH, "<$file";
	my $start;
	my $end;
	my $allele = 1;
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@tmp = split("_", $split[0]);
		$chr = substr($tmp[0], 3);
		if(@tmp > 1) {
			$allele = $tmp[2];
		} else {
			$allele = 1;
		}
		if(!exists $last{$chr}{$allele}) {
			print OUT $line . "\n";
		} else {
			$start = &shift($chr, $allele, $split[1]);
			$end = &shift($chr, $allele, $split[2]);
			if($end < $start) { $end = $start + 1; }
			print OUT $tmp[0] . "\t" . $start . "\t" . $end . "\t" . $split[3] . "\n";
		}
	}
	close FH;
	close OUT;
}

sub shift_bismark{
	my $file = $_[0];
	my $strain = $_[1];
	my $out_file = $_[2];
	print STDERR "shifting " . $file . "\n";
	open OUT, ">$out_file" or die "Can't open file!\n";
	open FH, "<$file";
	my $start;
	my $end;
	my $allele = 1;
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		@tmp = split("_", $split[0]);
		$chr = substr($tmp[0], 3);
		if(@tmp > 1) {
			$allele = $tmp[2];
		} else {
			$allele = 1;
		}
		if(!exists $last{$chr}{$allele}) {
			print OUT $line . "\n";
		} else {
			$start = &shift($chr, $allele, $split[1]);
			print OUT $tmp[0] . "\t" . $start;
			for(my $i = 2; $i < @split; $i++) {
				print OUT "\t" . $split[$i];
			}
			print OUT "\n";
		}
	}
	close FH;
	close OUT;
}
