#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Storable;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
use system_interaction;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
use Set::IntervalTree;

use threads;
use Data::Dumper;

$_ = "" for my($snp, $indel, @chr, $current_chr, $data, $filename, $out, $genome);
$_ = 0 for my($filter, $same, $help, $homo, $length_mut, $length_ref, $force, $lines_f1, $lines_f2, $lines_all, $line_count, $num_strains, $header_exists, $outfile_open, $last_h, $run_through, $none_number_chromosome, $core);
$_ = () for my($lines, @merge_line, $header, @split, $snps, $indels, $f1, $f2, %strains_to_use, @mut_files, @header, %last_shift_pos_ref, %last_shift_pos_strain, @strains_to_use, @s, @i, $check_snp, $check_indel, $s_c, $i_c, @fileHandlesMutation, @fileHandlesMutation2, @fileHandlesShift, @fileHandlesShift2, @ref_pos, @strain_pos, @ref_pos2, @strain_pos2, @pos_ref, @pos_strain, %lookup_no_number, %lookup_number, %multi_strains);
#Config with default data folder

$none_number_chromosome = 1000;

sub printCMD{
	print STDERR "Usage:\n";
	print STDERR "\n\nInput files:\n";
	print STDERR "\t-files <file>: Files with mutations - comma separated list (max 2 - one file for snps, one file for indels)\n";
	print STDERR "\t-homo: Assumes phenotype is homozygouse - if there are two different variations reported just the first one is taken into acocunt\n";
	print STDERR "\t-dir: output directory for genome folders - default: folder specified in config file\n"; 
	print STDERR "\t-strains <list of strains>: comma separated list of strains to include - when empty every strain in vcf file is considered\n";
	print STDERR "\n\n";
	print STDERR "Filter vcf file:\n";
	print STDERR "\t-filter - Filters out all mutations that did not pass all filters\n"; 	
	print STDERR "\t-same - Filters out all mutations that are homozyous\n";
	print STDERR "\tif position is homozygous and -same is not defined the first allele is taken without any further evaluation\n";
	print STDERR "\t-force: Overwrites existing folder\n";
	print STDERR "\n\n";
	print STDERR "\t-genome <path to reference fastq files>: generates strain specific genome\n\n\n";
	print STDERR "-h | --help: shows help\n\n\n";
	print STDERR "-core <number of cores>: default 1\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

print STDERR "Add generate genome as external method\n";

#Readin command line parameters
GetOptions(	"files=s{,}" => \@mut_files,
		"dir=s" => \$data,
		"strains=s{,}" => \@strains_to_use,
		"filter" => \$filter,
		"same" => \$same,
		"homo" => \$homo,
		"h" => \$help,
		"force" => \$force,
		"genome=s" => \$genome,
		"core=s" => \$core,
		"help" => \$help)
or &printCMD(); 

if($core == 0) {
	$core = 1;
} else {
	print STDERR "Using multithreading for IO writing\n";
}	
$num_strains = @strains_to_use;
for(my $i = 0; $i < @strains_to_use; $i++) {
	$strains_to_use[$i] =~ s/,//g;
	$multi_strains{uc($strains_to_use[$i])} = 1;
	$strains_to_use{uc($strains_to_use[$i])} = 1;
}

if($data eq "") {
	$data = config::read_config()->{'data_folder'};
}

if(@mut_files == 0) {
	print STDERR "No mutation files specified.\n";
}

if($help == 1) {
	&printCMD();
}

#define header
if(@mut_files > 0) {
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
		$multi_strains{$header[$i]} = 1;
		&check_file_existance($i);
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

if($homo == 1) {
	$a = 1;
} else {
	$a = 2;
}

for(my $i = 0; $i < @header - 9; $i++) {
	$pos_ref[$i] = 0;
	$pos_strain[$i] = 0;
	$chr[$i] = 0;
	$header[$i+9] = uc($header[$i+9]);
}

if($core > 1) {
	my %data_for_threading;
	for(my $i = 0; $i < $a; $i++) {
		for(my $h = 0; $h < @header - 9; $h++) {
			if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
				next;
			}
			$data_for_threading{$i}{$h} = 1;
		}
		
	}
	my @running = ();
	my @Threads;
	my $current_thread;
	my $current_thread_level1;
	my $current_thread_level2;
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
	while (scalar @running > 0) {
		foreach my $thr (@Threads) {
			$thr->join if ($thr->is_joinable());
		}
		@running = threads->list(threads::running);
	}
} else {
	for(my $i = 0; $i < $a; $i++) {
		for(my $h = 0; $h < @header - 9; $h++) {
			if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
				next;
			}
			&thread_routine($i, $a, $h);
	#		foreach my $l (@{$lines->[$i]->[$h]}) {
	#			@split = split('\t', $l);
	#			$out = $split[1] . "\t" . $split[2] . "\t" . $split[3];
	#			if($split[0] !~ /\d+/) {
	#				#Time to convert letters into numbers
	#				if(!exists $lookup_no_number{$split[0]}) {
	#					$lookup_no_number{$split[0]} = $none_number_chromosome;
	#					$lookup_number{$none_number_chromosome} = $split[0];
	#					$none_number_chromosome++;
	#				}
	#				$split[0] = $lookup_no_number{$split[0]};
	#			}
	#			if($chr[$h] ne $split[0]) {
	#				if($outfile_open == 1) {
	#					&close_filehandles($h);
	#					&create_offset_ref_to_strain($chr[$h], ($i + 1), $header[$h+9]);
	#					&create_offset_strain_to_ref($chr[$h], ($i + 1), $header[$h+9]);
	#				}
	#				$chr[$h] = $split[0];
	#				$pos_ref[$h] = $pos_strain[$h] = $split[1];
	#				$chr[$h] = $split[0];
	#				&open_filehandles($h);	
	#			}
	#			$fileHandlesMutation[0]->print($out . "\n");
	#			if($split[2] ne "" && $split[3] ne "" && length($split[2]) != length($split[3])) {
	#				($pos_ref[$h], $pos_strain[$h]) = &shift_vector($split[2], $split[3], $fileHandlesShift[0], $pos_ref[$h], $pos_strain[$h], $split[1]);
	#			}
	#		}
	#		$last_h = $h;
	#		$run_through++;
	#		&close_filehandles($last_h);
	#		&create_offset_ref_to_strain($chr[$last_h], $a, $header[$last_h+9]);
	#		&create_offset_strain_to_ref($chr[$last_h], $a, $header[$last_h+9]);
	#		&write_last_shift();
	#		$outfile_open = 0;
		}
	}
}
print STDERR "Data is stored in $data\n";
print STDERR "Processing data successfully finished!\n";


if($genome ne "") {
	print STDERR "Generating genomes per strain\n";
	my @genome_files = `ls $genome`;
	my $chr;
	my $mut_file;
	foreach my $g_file (@genome_files) {
		$chr = substr($g_file, 3, length($g_file) - 7);
		if(exists $lookup_number{$chr}) {
			$chr = $lookup_number{$chr};
		}
		for(my $h = 0; $h < @header - 9; $h++) {
			for(my $i = 0; $i < $a; $i++) {
				$filename = $data . "/" . $header[$h+9] . "/chr" . $chr . "_allele_" . ($i + 1) . ".fa";
				$mut_file = $data . "/" . $header[$h+9] . "/chr" . $chr . "_allele_" . ($i + 1) . ".mut";	
				if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
					next;
				}
				processing::create_genome($chr, $filename, $genome . "/" . $g_file, $mut_file);
			}

		}
	}
}


#SUBFUNCTIONS
#Routine to start threading
sub thread_routine {
	my $i = $_[0];
	my $a = $_[1];
	my $h = $_[2];
	my $outfile_open = 0;	
	my %thread_last_shift_ref;
	my $thread_ref_ref;
	my %thread_last_shift_strain;
	my $thread_ref_strain;
	my %thread_lookup;
	my @split;
	my $out;
	if($num_strains > 0 && !exists $strains_to_use{$header[$h+9]}) {
		next;
	}
	print STDERR "Processing " . uc($header[$h+9]) . "\n";
	foreach my $l (@{$lines->[$i]->[$h]}) {
		@split = split('\t', $l);
		$out = $split[1] . "\t" . $split[2] . "\t" . $split[3];
		if($split[0] !~ /\d+/) {
			#Time to convert letters into numbers
			if(!exists $lookup_no_number{$split[0]}) {
				$lookup_no_number{$split[0]} = $none_number_chromosome;
				$lookup_number{$none_number_chromosome} = $split[0];
				$none_number_chromosome++;
			}
			$split[0] = $lookup_no_number{$split[0]};
		}
		if($chr[$h] ne $split[0]) {
			if($outfile_open == 1) {
				&close_filehandles($h);
				$thread_ref_ref = &create_offset_ref_to_strain($chr[$h], ($i + 1), $header[$h+9], \%thread_last_shift_ref);
				%thread_last_shift_ref = %$thread_ref_ref;
				$thread_ref_strain = &create_offset_strain_to_ref($chr[$h], ($i + 1), $header[$h+9], \%thread_last_shift_strain);
				%thread_last_shift_strain = %$thread_ref_strain;
			}
			$chr[$h] = $split[0];
			$pos_ref[$h] = $pos_strain[$h] = $split[1];
			$chr[$h] = $split[0];
			&open_filehandles($h);
			$outfile_open = 1;
		}
		$fileHandlesMutation[0]->print($out . "\n");
		if($split[2] ne "" && $split[3] ne "" && length($split[2]) != length($split[3])) {
			($pos_ref[$h], $pos_strain[$h]) = &shift_vector($split[2], $split[3], $fileHandlesShift[0], $pos_ref[$h], $pos_strain[$h], $split[1]);
		}
	}
	$last_h = $h;
	$run_through++;
	&close_filehandles($last_h);
	$thread_ref_ref = &create_offset_ref_to_strain($chr[$last_h], $a, $header[$last_h+9], \%thread_last_shift_ref);
	%thread_last_shift_ref = %$thread_ref_ref;
	$thread_ref_strain = &create_offset_strain_to_ref($chr[$last_h], $a, $header[$last_h+9], \%thread_last_shift_strain);
	%thread_last_shift_strain = %$thread_ref_strain;
	&write_last_shift(\%thread_last_shift_strain, \%thread_last_shift_ref, $header[$last_h+9]);
	$outfile_open = 0;
}

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
                        #we are already at the first position of this deletion, so just add length(strain) - 1
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

sub create_offset_ref_to_strain {
        my $chr = $_[0];
	my $a = $_[1];
	my $strain = $_[2];
	my %last_shift_pos_ref = %{$_[3]};
	if(exists $lookup_number{$chr}) {
		$chr = $lookup_number{$chr};
	}
        $filename = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".mut";
	my $out = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".ref_to_strain.vector";
	$_ = () for my(@array_shift, @tmp_split);
        if(-e $filename) {
                open FH, "<$filename";
		open OUT, ">$out";
                my $run = 0;
		my $last_shift = 0;
                my @a;
                foreach my $line (<FH>) {
                        chomp $line;
                        @tmp_split = split('\t', $line);
			if(length($tmp_split[1]) == length($tmp_split[2])) { next; }
		#	if($run == $tmp_split[0]) {
		#		next;
		#	} 
			print OUT $last_shift . "\t" . $run . "\t" . $tmp_split[0]. "\n";
			$last_shift = $last_shift + (length($tmp_split[2]) - length($tmp_split[1]));
			$run = $tmp_split[0] + 1;
			$last_shift_pos_ref{$strain}{$chr}{$a}{'pos'} = $tmp_split[0] + 1;
			$last_shift_pos_ref{$strain}{$chr}{$a}{'shift'} = $last_shift;
                }
		close OUT;
        }
	return \%last_shift_pos_ref;
}

sub create_offset_strain_to_ref {
	my $chr = $_[0];
	my $a = $_[1];
	my $strain = $_[2];
	my %last_shift_pos_strain = %{$_[3]};
	if(exists $lookup_number{$chr}) {
		$chr = $lookup_number{$chr};
	}
        my $last_shift = 0;
	$_ = () for my (@array_shift, @tmp_split);
	my $start = 0;
	$filename = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".mut";
	my $out = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".strain_to_ref.vector";
	if(-e $filename) {
		open FH, "<$filename";
		open OUT, ">$out";
		foreach my $line (<FH>) {
			chomp $line;
			@tmp_split = split('\t', $line);
		#	while($start < $tmp_split[1] - 1) {
		#		$array_shift[$start] = $last_shift;
		#		$start++;
		#	}
			#tmp_split[0] is reference position - tmp_split[1] is strains position - calculate shifting vector out of it
			if($start == $tmp_split[0] - 1) {
				next;
			}	
			if(length($tmp_split[1]) == length($tmp_split[2])) { next; }
			print OUT $last_shift . "\t" . $start . "\t" . ($tmp_split[0] + length($tmp_split[2])) . "\n";
			$start = $tmp_split[0] + length($tmp_split[2]) + 1;
			$last_shift = $last_shift + (length($tmp_split[1]) - length($tmp_split[2]));
		#	$array_shift[$start] = $last_shift;
		#	$start++;
			$last_shift_pos_strain{$strain}{$chr}{$a}{'pos'} = $start;
			$last_shift_pos_strain{$strain}{$chr}{$a}{'shift'} = $last_shift;
		}
		close FH;
		close OUT;
	#	store \@array_shift, "$out";
	}
	return \%last_shift_pos_strain;
}

sub write_last_shift{
	my $out = "";
	my %last_shift_pos_strain = %{$_[0]};
	my %last_shift_pos_ref = %{$_[1]};
	foreach my $strain (keys %last_shift_pos_strain) {
		$out = $data . "/" . $strain . "/last_shift_strain.txt";
		store \%{$last_shift_pos_strain{$strain}}, "$out";
		$out = $data . "/" . $strain . "/last_shift_ref.txt";
		store \%{$last_shift_pos_ref{$strain}}, "$out";
		$out = $data . "/" . $strain . "/lookup_table_chr.txt";
		store \%lookup_no_number, "$out";
	}
}

sub check_file_existance{
	my $header_number = $_[0];
	        #This is the header
        #Save that a header exists
        $header_exists = 1;
        #Now go through all strains and see if the folder already exists
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
}

sub open_filehandles{
	my $header_number = $_[0];
	my $chr = $chr[$header_number];
	if(exists $lookup_number{$chr}) {
		$chr = $lookup_number{$chr};
	}
	print STDERR "\t\tchromosome " . $chr . "\n";
	$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr . "_allele_1.mut";
	open my $fh_mut, ">", "$filename" or die "Can't open $filename: $!\n";
	$fileHandlesMutation[0] = $fh_mut;
	$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr . "_allele_1.shift";
	open my $fh_shift, ">", "$filename" or die "Can't open $filename: $!\n";
	$fileHandlesShift[0] = $fh_shift;	
	$ref_pos[$header_number] = 0;
	$strain_pos[$header_number] = 0;
	if($homo == 0) {
		$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr . "_allele_2.mut";
		open my $fh_mut2, ">", "$filename" or die "Can't open $filename: $!\n";
		$fileHandlesMutation2[0] = $fh_mut2;
		$filename = $data . "/" . $header[$header_number+9] . "/chr" . $chr . "_allele_2.shift";
		open my $fh_shift2, ">", "$filename" or die "Can't open $filename: $!\n";
		$fileHandlesShift2[0] = $fh_shift2;
		$ref_pos2[$header_number] = 0;
		$strain_pos2[$header_number] = 0;
	}
        $outfile_open = 1;
}

sub close_filehandles{
	my $header_number = 0;
	close $fileHandlesMutation[$header_number];
	close $fileHandlesShift[$header_number];
	if($homo == 0) {
		close $fileHandlesMutation2[$header_number];
		close $fileHandlesShift2[$header_number];
	}
        undef $fileHandlesMutation[$header_number];
        undef $fileHandlesShift[$header_number];
        if($homo == 0) {
                undef $fileHandlesMutation2[$header_number];
                undef $fileHandlesShift2[$header_number];
        }
}

sub merge{
	$line_count++;
	print STDERR "Status: " . int(($line_count/$lines_all)*100) . "% Completed\r";
	my $line = $_[0];
	if($line eq "") {
		return;
	}
	chomp $line;
	my @current = split('\t', $line);
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
		if($num_strains > 0 && !exists $strains_to_use{$header[$i]}) {
			next;
		}
		@all = split('/', (split(":", $current[$i]))[0]);
		if($all[0] eq ".") {
			$all[0] = 0;
			$all[1] = 0;
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
	if ($fh and my $line = <$fh>) {
		chomp $line;
#		return [ split('\t', $line) ];
		return $line;
	}
	return;
}
