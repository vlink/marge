#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Storable;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
use Set::IntervalTree;
use threads;
use Data::Dumper;

$_ = "" for my($snp, $indel, @chr, $current_chr, $data, $filename, $out, $genome, $merge_line, $genome_dir, $ref_name, $sort_check);
$_ = 0 for my($filter, $same, $help, $hetero, $add, $force, $lines_f1, $lines_f2, $lines_all, $line_count, $num_strains, $outfile_open, $last_h, $none_number_chromosome, $core, $no_genome, $wait);
$_ = () for my($lines, @merge_line, $header, @split, $snps, $indels, $f1, $f2, %strains_to_use, @mut_files, @header, @strains_to_use, @s, @i, %lookup_no_number, %lookup_number, @last, @allele, @all, @test_spaces);

$none_number_chromosome = 1000;

sub printCMD{
	print STDERR "Usage:\n";
	print STDERR "\nInput files:\n";
	print STDERR "\t-files <file>: Files with mutations in VCF format\n";
	print STDERR "\t\teither one file containing SNPs and InDels or two files (one with SNPs, one with InDels)\n";
	print STDERR "\t\twhen using two files the order is SNP_file, InDel_file (either comma separated or with white space)\n";
	print STDERR "\nAdditional parameters\n";
	print STDERR "\t-ind <list of individuals>: comma separated list of strains to include - when empty every strain in VCF file is considered\n";
	print STDERR "\t-hetero: Assumes phenotype is heterozygous - default: homozygous\n";
	print STDERR "\t-genome: path to fasta files per chromosome: input directory of the reference genome per chromosome in fasta file format\n";
	print STDERR "\t-no-genome: Does not create strains specific genomes (default: off)\n";
	print STDERR "\t-ref <name>: Software generates a reference genome folder with all necessary files for further analysis - Folder is called REFERENCE, if nothing is specified here\n";
	print STDERR "\t-force: Overwrites existing folder\n";
	print STDERR "\t-add: Adds data to existing folder - if file exists it is overwritten\n";
	print STDERR "\t-core <number of cores>: default 1\n";
	print STDERR "\nConfig parameters - these parameters are defined in the config, but can be changed\n";
	print STDERR "\tCAUTION: If these parameters are changed, the config file either needs to be adjusted or they need to be defined for every script\n";
	print STDERR "\t-dir: output directory for mutation files for strains folders - default: folder specified in config file\n";
	print STDERR "\t-genome_dir: output directory for fasta genome file per strain - default: folder specified in config file\n";
	print STDERR "\nFilter VCF files:\n";
	print STDERR "\tFor a more sophisticated filtering - use vcftools\n";
	print STDERR "\t-filter - Filters out all mutations that did not pass all filters\n"; 	
	print STDERR "\t-same - Filters out all mutations that are homozygous\n";
	print STDERR "\t\twhen position is homozygous and -same is not defined the first allele is taken without any further evaluation\n";
	print STDERR "\n\n";
	print STDERR "\t-h | --help: shows help\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-files' => 1);
#my %mandatory = ('-files' => 1, '-ind' => 1, '-genome' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

#Readin command line parameters
GetOptions(	"files=s{,}" => \@mut_files,
		"dir=s" => \$data,
		"genome_dir=s" => \$genome_dir,
		"no-genome" => \$no_genome,
		"ind=s{,}" => \@strains_to_use,
		"ref=s" => \$ref_name,
		"filter" => \$filter,
		"same" => \$same,
		"hetero" => \$hetero,
		"h" => \$help,
		"force" => \$force,
		"add" => \$add,
		"genome=s" => \$genome,
		"core=s" => \$core,
		"help" => \$help)
or die(&printCMD()); 

#Set up the default parameters
if($core == 0) {
	$core = 1;
} else {
	print STDERR "Using multithreading for IO writing\n";
}
	
$num_strains = @strains_to_use;
for(my $i = 0; $i < @strains_to_use; $i++) {
	$strains_to_use[$i] =~ s/,//g;
	$strains_to_use{uc($strains_to_use[$i])} = 1;
}

if($core > $num_strains) {
	$core= $num_strains;
}
if($ref_name eq "") {
	$ref_name = "REFERENCE";
}

if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}
if($genome_dir ne "" && $genome eq "") {
	print STDERR "Specify input directory for the reference genome!\n";
	exit;
}

if($genome_dir eq "") {
	$genome_dir = $data;
}
if($help == 1) {
	&printCMD();
}


#define header
if(@mut_files > 0) {
	@test_spaces = split('\t', $mut_files[0]);
	if(@test_spaces > @mut_files) {
		@mut_files = @test_spaces;
	}
	$mut_files[0] =~ s/,//g;
	$header = `grep -m 1 \'#CHROM\' $mut_files[0]`;
	chomp $header;
	if($header eq "") {
		print STDERR "Mutation file does not have a header - naming of the strains not possible!\n";
		exit;
	}
	@header = split /[\t\s]+/, uc($header);
	for(my $i = 9; $i < @header; $i++) {
		if(@strains_to_use > 0 && !exists $strains_to_use{$header[$i]}) {
			next;
		}
		&check_file_existance($header[$i]);
	}

}
if($hetero == 1) {
	$a = 2;
} else {
	$a = 1;
}
#Define chromsome, allele and header
for(my $i = 0; $i < @header - 9; $i++) {
	$chr[$i] = 0;
	$allele[$i] = 0;
	$header[$i+9] = uc($header[$i+9]);
}

my $min_one = 0;
for(my $i = 9; $i < @header; $i++) {
	if(exists $strains_to_use{$header[$i]}) {
		$min_one++;
	}
}
if($min_one == 0) {
	print STDERR "None of the individuals you specified is in the VCF files!\n";
	print STDERR "Here are all individuals specified in the header of your file:\n";
	for(my $i = 9; $i < @header; $i++) {
		print STDERR "\t" . $header[$i] . "\n";
	}
	print STDERR "Processing all individuals in 5 seconds\n";
	print STDERR "Print Ctrl + C to abort\n";
	for(my $i = 0; $i < 5; $i++) {
		print STDERR ".";
		sleep(1);
	}
	print STDERR "\n";
	print STDERR "Starting to process all individuals\n";
	for(my $i = 9; $i < @header; $i++) {
		$strains_to_use{$header[$i]} = 1;
	}
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
	print STDERR "Reading in file and caching!\n";
	#Merge INDEL and SNP files
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
	print STDERR "\n\n";
}

#Set up routine for threading
if($core > 1) {
	$_ = () for my (%data_for_threading, @running, @Threads, $current_thread, $current_thread_level1, $current_thread_level2);
	#Add threads to queue
	for(my $i = 0; $i < $a; $i++) {
		for(my $h = 0; $h < @header - 9; $h++) {
			if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
				next;
			}
			$data_for_threading{$i}{$h} = 1;
		}
		
	}
	#Run through queue and write strains data
	while (scalar @Threads < $num_strains) {
		@running = threads->list(threads::running);
		if(scalar @running < $core) {
			foreach my $keys (keys %data_for_threading) {
				foreach my $keys2 (keys %{$data_for_threading{$keys}}) {
					$current_thread = $keys . "_" . $keys2;
					$current_thread_level1 = $keys;
					$current_thread_level2 = $keys2;
				}
			}
			delete $data_for_threading{$current_thread_level1}{$current_thread_level2};
			#Start new thread and execute sub-function thread-routine
			my $thread = threads->new( sub { thread_routine($current_thread_level1, $a, $current_thread_level2) });
			push (@Threads, $thread);
			my $tid = $thread->tid;
			@running = threads->list(threads::running);
		}
		@running = threads->list(threads::running);
		foreach my $thr (@Threads) {
			if($thr->is_running()){ 
				my $tid = $thr->tid;
			} elsif($thr->is_joinable()){
				my $tid = $thr->tid;
				$thr->join;
			}
		}
		@running = threads->list(threads::running);
	}	
	#Wait for the remaining threads to finish
	while (scalar @running > 0) {
		foreach my $thr (@Threads) {
			$thr->join if ($thr->is_joinable());
		}
		@running = threads->list(threads::running);
	}
} else {
	#No threading
	#Execute sub-function thread routine iteratively for each strain and allele
	for(my $i = 0; $i < $a; $i++) {
		for(my $h = 0; $h < @header - 9; $h++) {
			if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
				next;
			}
			&thread_routine($i, $a, $h);
		}
	}
}

print STDERR "Data is stored in $data\n";
print STDERR "Processing data successfully finished!\n";

#Create a folder for the reference, so in downstream analysis there is a folder the program can use
&check_file_existance(uc($ref_name));
#If strains do not have mutations still create files necessary for downstream analysis
&touch_last_shift(\%strains_to_use, uc($ref_name));

#Generate genomes if varibale is set
if($no_genome == 0) {
	print STDERR "Generating genomes per strain\n";
	print STDERR $genome . "\n";
	my @genome_files = `ls $genome/*fa`;
	$_ = () for my(@tmp, $chr, $mut_file, $command);
	foreach my $g_file (@genome_files) {
		chomp $g_file;
		@tmp = split('/', $g_file);
		$chr = substr($tmp[-1], 3, -3);
		if(exists $lookup_number{$chr}) {
			$chr = $lookup_number{$chr};
		}
		for(my $i = 0; $i < $a; $i++) {
			#Copy fastq file to reference genome for downstream analysis
			$command = "cp " . $g_file . " " . $genome_dir . "/" . uc($ref_name) . "/chr" . $chr . "_allele_" . ($i + 1) . ".fa";
			`$command`;
		}
		#Generate genome
		processing::create_genome($chr, \%strains_to_use, $data, $genome_dir, $genome . "/chr" . $chr . ".fa", $a, $hetero);
	}
}


#SUBFUNCTIONS
#Routine to start threading
sub thread_routine {
	#Define variables new - method is called several times in parallel - varibales need to be local
	#I is the current allele we are looking at
	my $i = $_[0];
	my $allele = $i + 1;
	my $a = $_[1];
	my $h = $_[2];
	$_ = 0 for my($last_h, $outfile_open);	
	$_ = () for my (%thread_last_shift_ref, $thread_ref_ref, %thread_last_shift_strain, $thread_ref_strain, %thread_lookup, @split, $out, $fh_local);
	if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
		next;
	}
	print STDERR "Processing " . uc($header[$h+9]) . "\n";
	#Run through all merged lines for strain $h
	foreach my $l (@{$lines->[$i]->[$h]}) {
		@split = split('\t', $l);
		$out = $split[1] . "\t" . $split[2] . "\t" . $split[3];
		if($split[0] !~ /\d+/) {
			#Time to convert letters into numbers
			if(!exists $lookup_no_number{$split[0]}) {
				$lookup_no_number{$split[0]} = $none_number_chromosome;
				$thread_lookup{$none_number_chromosome} = $split[0];
				$none_number_chromosome++;
			}
		}
		#The new line is a different chromsome, so we have to open a new file
		if($chr[$h] ne $split[0] || $allele[$h] ne $allele) {
			if($outfile_open == 1) {
				&close_filehandles($fh_local);
				#Create offset vectors to save them for downstream analysis
				$thread_ref_ref = processing::create_offset_ref_to_strain($chr[$h], $allele, $header[$h+9], \%thread_last_shift_ref, $data, \%thread_lookup);
				%thread_last_shift_ref = %$thread_ref_ref;
				$thread_ref_strain = processing::create_offset_strain_to_ref($chr[$h], $allele, $header[$h+9], \%thread_last_shift_strain, $data, \%thread_lookup);
				%thread_last_shift_strain = %$thread_ref_strain;
			}
			$chr[$h] = $split[0];
			$allele[$h] = $allele;
			#Open new file handles for the current chromosome and allele
			$fh_local = &open_filehandles($h, $allele, $fh_local);
			$outfile_open = 1;
		}
		$fh_local->print($out . "\n");
	}
	$last_h = $h;
	if($outfile_open == 1) {
		&close_filehandles($fh_local);
	}
	$thread_ref_ref = processing::create_offset_ref_to_strain($chr[$last_h], $allele, $header[$last_h+9], \%thread_last_shift_ref, $data);
	%thread_last_shift_ref = %$thread_ref_ref;
	$thread_ref_strain = processing::create_offset_strain_to_ref($chr[$last_h], $allele, $header[$last_h+9], \%thread_last_shift_strain, $data);
	%thread_last_shift_strain = %$thread_ref_strain;
	#Save the last position for shifting in a seperate file
	&write_last_shift($thread_ref_strain, $thread_ref_ref);
	$outfile_open = 0;
}

#Generate the minimum set of files required for downstream analysis in case this strain did not have any mutations
sub touch_last_shift {
	my %strains = %{$_[0]};
	my $ref = $_[1];
	$strains{$ref} = 1;
	#Create empty hash that can be saved
	my %empty = ();
	foreach my $s (keys %strains) {
		print STDERR $s . "\n";
		next;
		if(!-e $data . "/" . $s . "/last_shift_strain.txt") {
			store \%empty, "$data/$s/last_shift_strain.txt";
			store \%empty, "$data/$s/last_shift_ref.txt";
			store \%lookup_no_number, "$data/$s/lookup_table_chr.txt";
		}
	} 
}

#Output the last shift files for all strains that were processed
sub write_last_shift{
	my $last_shift_pos_strain = $_[0];
	my $last_shift_pos_ref = $_[1];
	foreach my $strain (keys %{$last_shift_pos_strain}) {
		$out = $data . "/" . $strain . "/last_shift_strain.txt";
		store \%{$last_shift_pos_strain->{$strain}}, "$out";
		$out = $data . "/" . $strain . "/last_shift_ref.txt";
		store \%{$last_shift_pos_ref->{$strain}}, "$out";
		$out = $data . "/" . $strain . "/lookup_table_chr.txt";
		store \%lookup_no_number, "$out";
	}
}

#Check existance of folder to not overwrite by accident
sub check_file_existance{
	$header = $_[0];
        #Save that a header exists
        #Now go through all strains and see if the folder already exists
	if(-e $data . "/" . $header) {
		print STDERR "Folder for " . $header . " already exists!\n";
		if($force == 1) {
			$wait = 1;
		} elsif($add == 1) {
			$wait = 1;
		} else {
			print STDERR "If you want to overwrite this folder use -force!\n";
			exit;
		}
	} else {
		`mkdir -p $data/$header`;
	}
}

if($wait == 1) {
	if($force == 1) {
		print STDERR "Deleting folder " . $header . "!\n";
		print STDERR "Waiting for 10 seconds\n";
		print STDERR "Press Ctrl + C to interrupt\n";
		for(my $j = 0; $j < 10; $j++) {
			print STDERR ".";
		       sleep(1);
		}
		print STDERR "\n";
		`rm -rf $data/$header/*`;
	} elsif($add == 1) {
		print STDERR "If data already exists in this folder it will be overwritten!\n";
		print STDERR "\nWaiting for 3 seconds\n";
		print STDERR "Press Ctrl + C to interrupt\n";
		for(my $j = 0; $j < 3; $j++) {
			print STDERR ".";
			sleep(1);
		}
		print STDERR "\n";
	}
}
#Open filehandles and save them in a array
sub open_filehandles{
	my $header_number = $_[0];
	my $allele_number = $_[1];
	my $fh = $_[2];
	my $chr = $chr[$header_number];
	if(exists $lookup_number{$chr}) {
		$chr = $lookup_number{$chr};
	}
	print STDERR "\t\tchromosome " . $chr . " allele " . $allele_number . "\n";
	$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr . "_allele_" . $allele_number . ".mut";
	open $fh, ">", "$filename" or die "Can't open $filename: $!\n";
	$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr . "_allele_" . $allele_number . ".shift";
        $outfile_open = 1;
	return $fh;
}

sub close_filehandles{
	my $fh = $_[0];
	close $fh;
	undef $fh;
}

#Merge mutations from SNP and InDel file in one large file and shorten variances to the shortest possible solution per strain
sub merge{
	$line_count++;
	print STDERR "Status: " . int(($line_count/$lines_all)*100) . "% Completed\r";
	my $line = $_[0];
	if($line eq "") {
		return;
	}
	chomp $line;
	#Split current line of the file and split into different parts
	my @current = split('\t', $line);
	#Save reference and position
	my $pos = $current[1];
	my $ref = $current[3];
	$_ = () for my($var, $max, $copy, $pos_save);
	my $seq_before = "";
	$current[4] =~ s/\.//g;
	#The human VCF files contain more sophisticated annotations <CN>, <INS>
	#Just keep CN at the moment, filter out the other annotations
	if($current[4] =~ m/[\:\<\>]/) {
		return;
	}
	my @variants = split(',', $current[4]);
	#If CN is 0 add one basepair to the left, shift position one down, so the mutation that would be empty gets one basepair
	for(my $i = 0; $i < @variants; $i++) {
		if(substr($variants[$i], 0, 3) eq "<CN") {
			$copy = substr($variants[$i], 3, -1);
			if($copy == 0) {
				$seq_before = processing::get_reference_position($genome . "/chr" . $current[0] . ".fa", $current[0], $pos - 1);
			}
			$variants[$i] = ($ref) x $copy;
		} elsif(substr($variants[$i], 0, 1) eq "<") {
			last;
		}
	}
	#Add the base before the actual mutation, because CN is 0
	if($seq_before ne "") {
		$pos--;
		$ref = $seq_before . $ref;
		for(my $i = 0; $i < @variants; $i++) {
			$variants[$i] = $seq_before . $variants[$i];
		}
	}
	$pos_save = $pos;
	#Add reference variance at the beginning of the variants array (if strain has reference it is annotated as 0)
	unshift @variants, $ref;
	if($filter == 1 && $current[6] ne "PASS") {
		return 1;
	}
	#Run through the strains
	for(my $i = 9; $i < @current; $i++) {
		if($num_strains > 0 && !exists $strains_to_use{$header[$i]}) {
			next;
		}
		@all = split /[\/,\|]+/, ((split(":", $current[$i]))[0]);
		#No mutation can also be annotated as . - reference allele is in variants[0] s0 change , to 0
		if($all[0] eq ".") {
			$all[0] = 0;
			$all[1] = 0;
		}
		#Sometimes a variance for the second allele is annotated - if that is . change it to 0
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
		#organism is homozygous - one varaint can be ignored
		if($hetero == 0) {
			shift @all;
		}
		@merge_line = ();
		$merge_line = "";
		#Run though variance annotation for the number of alleles we are looking at 
		for(my $h = 0; $h < @all; $h++) {
			$var = $variants[$all[$h]];
			$ref = $variants[0];
			$pos = $pos_save;
			#Check if reference and mutations are different
			if($ref ne $var) {
				if(length($ref) > 1 && length($var) > 1) {
					#Remove the last bases if they are the same
					while(substr($ref, -1) eq substr($var, -1) && length($ref) > 1 && length($var) > 1) {
						chop $ref;
						chop $var;
					}
					#Now start from the beginning
					$max = (length($ref) > length($var)) ? length($var) : length($ref);	
					for(my $l = 0; $l < $max - 1; $l++) {
						if(length($ref) > 1 && length($var) > 1 && substr($ref, 0, 1) eq substr($var, 0, 1)) {
							$pos++;
							$ref = unpack "xA*", $ref;
							$var = unpack "xA*", $var;
						}
					}
				}
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

#Read in file line by line
sub read_file_line {
	my $fh = shift;
	if ($fh and my $line = <$fh>) {
		chomp $line;
		return $line;
	}
	return;
}


