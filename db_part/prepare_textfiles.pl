#!/usr/bin/perl -w
use strict;
use Getopt::Long;
require '../general/config.pm';
require '../general/system_interaction.pm';
use threads;
use Data::Dumper;

$_ = "" for my($genome, $ref_genome, $annotation, $genes, $name, $current_line, $last_line, $snp, $indel, @chr, $current_chr);
$_ = 0 for my($filter, $same, $core_user, $multithreading, $resume, $help, $check_muts, $check_ref, $check_annotation, $check_genes, $change, $homo, $chr_name_length, $length_mut, $length_ref, $force, $lines_f1, $lines_f2, $lines_all, $line_count);
my($num, $lines, @merge_line, @end, %hash, $dbh, $sth, $header, %check_offset, @split, $snps, $indels, @s, @ref, $ref, @i, $table_name, $f1, $f2, $command, %check_strains, %strains_to_use, @mut_files, $total_chr, %seen, @header);
my $data = config::read_config()->{'data_folder'};


sub printCMD{
	print STDERR "Usage:\n";
	print STDERR "\n\nInput files:\n";
	print STDERR "\t-files <file>: Files with mutations - comma spearated list (max 2 - one file for snps, one file for indels)\n";
	print STDERR "\t-homo - Assumes phenotype is homozygouse - if there are two different variations reported just the first one is taken into acocunt\n";
	print STDERR "\n\n";
	print STDERR "Filter vcf file.\n";
	print STDERR "\t-filter - Filters out all mutations that did not pass all filters\n"; 	
	print STDERR "\t-same - Filters out all mutations that are homozyous\n";
	print STDERR "\tif position is homozygous and -same is not defined the first allele is taken without any further evaluation\n";
	print STDERR "\t-force: Overwrites existing tables\n";
	print STDERR "\n\n";
	print STDERR "-h | --help: shows help\n\n\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

#TODO add resume then clean up code and then done
GetOptions(	"files=s{,}" => \@mut_files,
		"filter" => \$filter,
		"same" => \$same,
		"homo" => \$homo,
		"h" => \$help,
		"force" => \$force,
		"help" => \$help)
or &printCMD(); 

if(@mut_files == 0) {
	print STDERR "No mutation files specified. Just adding genome to database!\n";
}

if($help == 1) {
	&printCMD();
}

#define header
if(@mut_files > 0) {
	$header = `grep -m 1 \'#CHROM\' $mut_files[0]`;
	chomp $header;
	if($header eq "") {
		print STDERR "Mutation file does not have a header - naming of the strains not possible!\n";
		exit;
	}
	@header = split /[\t\s]+/, uc($header);
}

#Merge SNP and indel file if both exists
if(@mut_files > 0) {
	print STDERR "Reading in file $mut_files[0]\n";
	#open files and filter out header
	$snp = $mut_files[0];
	$lines_f1 = `wc -l $mut_files[0]`;	
	open($f1, $snp);
	$snps = read_file_line($f1);
	while(substr($snps, 0, 1) eq "#") {
		$line_count++;
		$snps = read_file_line($f1);
		next;
	}
	@s = split('\t', $snps);

	if(@mut_files == 2) {
		print STDERR "Reading in file $mut_files[1]\n";
		$indel = $mut_files[1];
		open($f2, $indel);
		$lines_f2 = `wc -l $mut_files[1]`;
		$indels = read_file_line($f2);
		while(substr($indels, 0, 1) eq "#") {
			$line_count++;
			$indels = read_file_line($f2);
			next;
		}
		@i = split('\t', $indels);
	}
	print STDERR "Merge files and mutations!\n";
	$lines_all = (split('\s+', $lines_f1))[0] + (split('\s+', $lines_f2))[0];
	my $check_snp;
	my $check_indel;
	my $s_c;
	my $i_c;
	print STDERR "Reading in file and caching!\n";

	while($snps and $indels) {
		@s = split('\t', $snps);
		@i = split('\t', $indels);
		if($current_chr eq "") { $current_chr = $s[0]; }
		if($s[0] ne $current_chr && $i[0] ne $current_chr) {
			$current_chr = $s[0];
		}
		if($s[0] ne $current_chr) {
			&merge($indels);
			$indels = read_file_line($f2);		
		}
		if($i[0] ne $current_chr) {
			&merge($snps);
			$snps = read_file_line($f1);
		}
		if($s[0] eq $i[0]) {
			if($s[1] < $i[1]) {
				&merge($snps);
				$snps = read_file_line($f1);
			} else {
				&merge($indels);
				$indels = read_file_line($f2);
			}
		}
	}
	while($snps) {
		&merge($snps);
		$snps = read_file_line($f1);
	}
	while($indels) {
		&merge($indels);
		$indels = read_file_line($f2);
	}
}


#print Dumper $lines;
my $header_exists = 0;
my $outfile_open = 0;
my @fileHandlesMutation;
my @fileHandlesMutation2;
my @fileHandlesShift;
my @fileHandlesShift2;
my @ref_pos;
my @strain_pos;
my @ref_pos2;
my @strain_pos2;
my $filename;

if($homo == 1) {
	$a = 1;
} else {
	$a = 2;
}
print "\n";
my @pos_ref;
my @pos_strain;
for(my $i = 0; $i < @header - 9; $i++) {
	$pos_ref[$i] = 0;
	$pos_strain[$i] = 0;
	$chr[$i] = "";
	$header[$i+9] = uc($header[$i+9]);
}
my $out;


for(my $i = 0; $i < $a; $i++) {
	for(my $h = 0; $h < @header - 9; $h++) {
		print STDERR "Processing " . uc($header[$h+9]) . "\n";
		foreach my $l (@{$lines->[$i]->[$h]}) {
			@split = split('\t', $l);
			$out = $split[1] . "\t" . $split[2] . "\t" . $split[3];
			if($chr[$h] eq "") {
				&check_file_existance($header, $h);
			}
			if($chr[$h] ne $split[0]) {
				$chr[$h] = $split[0];
				if($outfile_open == 1) {
					&close_filehandles($h - 1);
				}
				$pos_ref[$h] = $pos_strain[$h] = $split[1];
				&open_filehandles($h);				
			}
			$fileHandlesMutation[0]->print($out . "\n");
			if(length($split[2]) != length($split[3])) {
				($pos_ref[$h], $pos_strain[$h]) = &shift_vector($split[2], $split[3], $fileHandlesShift[0], $pos_ref[$h], $pos_strain[$h], $split[1]);
			}
		}
	}
}
&close_filehandles();

print STDERR "Data is stored in $data\n";
print STDERR "Processing data successfully finished!\n";

#SUBFUNCTIONS
sub shift_vector{
	my $ref = $_[0];
	my $strain = $_[1];
	my $fileHandle = $_[2];
	my $pos_ref = $_[3];
	my $pos_strain = $_[4];
	my $pos_mut = $_[5];
	my $diff = $pos_ref - $pos_strain;
	$pos_ref = $pos_mut;
	$pos_strain = $pos_ref - $diff;
        if(length($ref) > length($strain)) {
                if(length($strain) > 1) {
                        #we are already at the firs tposition of this deletion, so just add length(strain) - 1
                        for(my $i = 1; $i < length($strain); $i++) {
                                $pos_ref++;
                                $pos_strain++;
                        }
		}
		#Next we jump in the reference sequence to the end of the insertion in the ref (deletion in strain), therefore we have to add 1 (because we look at the next position) and the length of the difference
		$pos_ref++;
		$pos_strain++;
		$pos_ref = $pos_ref + (length($ref) - length($strain));
		$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
        #Reference is shorter than mutated strain
        } else {
                if(length($ref) > 1) {
                        #Add the overlap. We are already at the first position of the insertion, so we just add -1 less
                        for(my $i = 0; $i < length($ref); $i++) {
                                $pos_strain++;
                                $pos_ref++;
                        }
		} else {
		#We are already at the current position of the mutation
			$pos_strain++;
			$pos_ref++;
		}
		for(my $i = 0; $i < length($strain) - length($ref); $i++) {
			$pos_strain++;
		}
		$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
        }
        return ($pos_ref, $pos_strain);
}

sub check_file_existance{
        my $line = $_[0];
	my $header_number = $_[1] + 9;
        #This is the header
        #Save that a header exists
        $header_exists = 1;
        #Now go through all strains and see if the folder already exists
#	for(my $i = 10; $i < @header; $i++) {
	if(-e $data . "/" . $header[$header_number]) {
		print STDERR "Folder for " . $header[$header_number] . " already exists!\n";
		if($force == 1) {
			print STDERR "Deleting folder " . $header[$header_number] . "!\n";
			print STDERR "Waiting for 10 seconds\n";
			print STDERR "Press Ctrl + C to interrupt\n";
			for(my $j = 0; $j < 10; $j++) {
				print STDERR ".";
			       sleep(1);
			}
			print STDERR "\n";
			`rm -rf $data/$header[$header_number]/*`;
		} else {
			print STDERR "If you want to overwrite this folder use -force!\n";
			exit;
		}
	} else {
		`mkdir $data/$header[$header_number]`;
	}
#        }
}

sub open_filehandles{
	my $header_number = $_[0];
	$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr[$header_number] . "_allele_1.mut";
	open my $fh_mut, ">", "$filename" or die "Can't open $filename: $!\n";
	$fileHandlesMutation[0] = $fh_mut;
#	push @fileHandlesMutation, $fh_mut;
	$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr[$header_number] . "_allele_1.shift";
	open my $fh_shift, ">", "$filename" or die "Can't open $filename: $!\n";
	$fileHandlesShift[0] = $fh_shift;	
#	push @fileHandlesShift, $fh_shift;
	$ref_pos[$header_number] = 0;
	$strain_pos[$header_number] = 0;
	if($homo == 0) {
		$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr[$header_number] . "_allele_2.mut";
		open my $fh_mut2, ">", "$filename" or die "Can't open $filename: $!\n";
	#	push @fileHandlesMutation2, $fh_mut2;
		$fileHandlesMutation2[0] = $fh_mut2;
		$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr[$header_number] . "_allele_2.shift";
		open my $fh_shift2, ">", "$filename" or die "Can't open $filename: $!\n";
	#	push @fileHandlesShift2, $fh_shift2;
		$fileHandlesShift2[0] = $fh_shift2;
		$ref_pos2[$header_number] = 0;
		$strain_pos2[$header_number] = 0;
	}
        $outfile_open = 1;
}

sub close_filehandles{
	my $header_number = 0;
#	for(my $i = 0; $i < @split - 10; $i++) {
	close $fileHandlesMutation[$header_number];
	close $fileHandlesShift[$header_number];
	if($homo == 0) {
		close $fileHandlesMutation2[$header_number];
		close $fileHandlesShift2[$header_number];
	}
#	}
        undef $fileHandlesMutation[$header_number];
        undef $fileHandlesShift[$header_number];
        if($homo == 0) {
                undef $fileHandlesMutation2[$header_number];
                undef $fileHandlesShift2[$header_number];
        }
}

my $merge_line = "";
my @len;
my $end = 0;
my $overlap = 0;
my @current;
my @points;
my @last = (0);
my $min_start;
my $number;
my $real_length_last;
my $real_length_current;
my $final_length;
my $n;

sub merge{
	$line_count++;
	print STDERR "Status: " . int(($line_count/$lines_all)*100) . "% Completed\r";
	my $line = $_[0];
	if($line eq "") {
		return;
	}
	chomp $line;
	@current = split('\t', $line);
	#Save reference and position
	my $pos = $current[1];
	my $ref = $current[3];
	my $merge_line;
	my @last;
	my @allele;
	my @all;
	$current[4] =~ s/\.//g;
	my @variants = split(',', $current[4]);

	unshift @variants, $current[3];
	if($filter == 1 && $current[6] ne "PASS") {
		return 1;
	}
	my $var;
	for(my $i = 9; $i < @current; $i++) {
		@all = split('/', (split(":", $current[$i]))[0]);
		if($all[0] eq ".") {
			$all[0] = 0;
		}
		if(@all > 1) { 
			if($all[1] eq ".") {
				$all[1] = 0;
			}
		}
		if($same == 1) {
			if(@all > 1 && $all[0] ne $all[1]) {
				next;
			}
		}
		if($homo == 1) {
			shift @all;
		}
		@merge_line = ();
		for(my $h = 0; $h < @all; $h++) {
			$var = $variants[$all[$h]];
			$ref = $current[3];
			$pos = $current[1];
			#Check if reference and mutations are different
			if($ref ne $var) {
				if(length($ref) > 1 && length($var) > 1) {
					#Remove the last bases if they are the same
					while(substr($ref, -1) eq substr($var, -1) && length($ref) > 1 && length($var) > 1) {
						chop $ref;
						chop $var;
					}
					#Now start from the beginning
					my $max = (length($ref) > length($var)) ? length($var) : length($ref);	
					for(my $l = 0; $l < $max - 1; $l++) {
						if(length($ref) > 1 && length($var) > 1 && substr($ref, 0, 1) eq substr($var, 0, 1)) {
							$pos++;
							$ref = unpack "xA*", $ref;
							$var = unpack "xA*", $var;
						}
					}
				}
				if(length($ref) > $length_mut) { $length_mut = length($ref); }
				if(length($var) > $length_mut) { $length_mut = length($var); }
				#Check if there was already a mutation before this one
				if(defined $lines->[$h]->[$i-9]) {
					@last = split('\t', $lines->[$h]->[$i-9]->[-1]);
				} else {
					$last[0] = 0;
					$last[1] = 0;
					$last[2] = "";
				}
				if($current[0] eq $last[0] && $current[1] < $last[1] + length($last[2])) {
					next;
				}
				$merge_line = $current[0] . "\t" . $pos . "\t" . $ref . "\t" . $var;
				push(@{ $lines->[$h]->[$i-9]}, $merge_line);
			}
		}
	}
}

#First step filter mutations if just one file is given
#If SNP and Indel file are given - merge the files
sub read_file_line {
	my $fh = shift;
	my @split;

	if ($fh and my $line = <$fh>) {
		chomp $line;
#		return [ split('\t', $line) ];
		return $line;
	}
	return;
}
