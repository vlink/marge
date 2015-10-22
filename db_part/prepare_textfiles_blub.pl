#!/usr/bin/perl -w
use strict;
use Getopt::Long;
require '../general/config.pm';
require '../general/system_interaction.pm';
use threads;
use Data::Dumper;

$_ = "" for my($genome, $ref_genome, $annotation, $genes, $name, $current_line, $last_line, $snp, $indel, $chr, $filename);
$_ = 0 for my($filter, $same, $core_user, $multithreading, $resume, $help, $check_muts, $check_ref, $check_annotation, $check_genes, $change, $homo, $chr_name_length, $length_mut, $length_ref, $force);
my($num, $lines, @merge_line, @end, %hash, $dbh, $sth, $header, %check_offset, @split, $snps, $indels, @s, @ref, @i, $table_name, $f1, $f2, $command, @h, %check_strains, %strains_to_use, @mut_files, $total_chr, %seen, @muts, @alleles);


sub printCMD{
	print STDERR "Usage:\n";
	print STDERR "General commands:\n";
	print STDERR "\n\nInput files:\n";
	print STDERR "\t-files <file>: Files with mutations - comma spearated list (max 2 - one file for snps, one file for indels)\n";
	print STDERR "\t-homo - Assumes phenotype is homozygouse - if there are two different variations reported just the first one is taken into acocunt\n";
	print STDERR "\t-force - overwrites strains folders when they already exist!\n";
	print STDERR "\n\n";
	print STDERR "Filter vcf file.\n";
	print STDERR "\t-filter - Filters out all mutations that did not pass all filters\n"; 	
	print STDERR "\t-same - Filters out all mutations that are homozyous\n";
	print STDERR "\tif position is homozygous and -same is not defined the first allele is taken without any further evaluation\n";
	print STDERR "\n\nMultithreading:\n";
	print STDERR "\t-no-multithreading: does not use multithreading to create vectors for database\n";
	print STDERR "\t-core <num>: Number of cores to use for multithreading - if not specified and multithreading is not turned off and possible it uses 1/4 of all cores\n";
	print STDERR "\n\n";
	print STDERR "-h | --help: shows help\n\n\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-files' => 1);
my %commandline = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%commandline);
#TODO add resume then clean up code and then done
GetOptions(	"files=s{,}" => \@mut_files,
		"same" => \$same,
		"core=s" => \$core_user,
		"no-multithreading" => \$multithreading,
		"homo" => \$homo,
		"h" => \$help,
		"help" => \$help, 
		"force" => \$force)
or &printCMD(); 

#Only one big vcf file
my $data = config::read_config()->{'data_folder'};


#Read in file

my $header_exists = 0;
my $outfile_open = 0;
my @header;
my @fileHandlesMutation;
my @fileHandlesMutation2;
my @fileHandlesShift;
my @fileHandlesShift2;
my @ref_pos;
my @strain_pos;
my @ref_pos2;
my @strain_pos2;

$snp = $mut_files[0];
open($f1, $snp);
$snps = read_file_file($f1);

if(@mut_files > 1) {
	$indel = $mut_files[1];
	open($f2, $indel);
	$indels = read_file_line($f2);
	while(substr($indels, 0, 1) eq "#") {
		$indels = read_file_line($f2);
		next;
	}
}


while($snps and $indels) {
	@snps = split('\t', $snps);
	@indels = split('\t', $indels);
	for(my $i = 0; $i < @snps - 10; $i++) {
		@alleles_snps = split('/', $snps[$i+10]);
		@alleles_indels = split('/', $indels[$i+10]);
		@muts_snps = split(",", $snps[4]);
		@muts_indel = split(",", $indles[4]);

		#Fix chromosomes here!
		if($snps[0] ne $indels[0]) {
			#Chromosome is different!
			if($snps[0] eq $chr) {
				$fileHandlesMutation[$i]->print($snps[1] . "\t" . $snps[3] . "\t" . $mut_snps[$alleles_snps[0] - 1] . "\n");
				$snps = read_file_line($f2);
			}
			if($indels[0] eq $chr) {
				$fileHandlesMutation[$i]->print($indels[1] . "\t" . $indels[3] . "\t" . $mut_indels[$alleles_indels[0] - 1] . "\n");
				$indels = read_file_line($f1);
			}
		} else {
			if($chr ne $snps[0]) {
				print STDERR "Preparing chromosome " . $snps[0] . "\n";
				$chr = $snps[0];
				if($outfile_open == 1) {
					&close_filehandles();
				}
				#open all files and set positions
				&open_filehandles();
			}

			if($snps[1] > $indels[1]) {
				if($alleles_indels[0] ne "0") {
					$fileHandlesMutation[$i]->print($indels[1] . "\t" . $indels[3] . "\t" . $mut_indel[$alleles_indels[0] - 1] . "\n");
				}
				if(length($indels[3]) != length($mut_indel[$alleles_indels[0] - 1])) {
					($ref_pos[$i], $strain_pos[$i]) = &shift_pos($indels[3], $mut_indel[$alleles_indels[0] - 1], $ref_pos[$i], $strain_pos[$i], $fileHandlesShift[$i], $indels[1]);
				}
				$indels = read_file_line($f2);
			} 
			if($snps[1] < $indels[1]) {
				if($alleles_snps[0] ne "0") {
					$fileHandlesMutation[$i]->print($snps[1] . "\t" . $snps[3] . "\t" . $mut_snps[$alleles_snps[0] - 1] . "\n");
				}
				if(length($snps[3]) != length($mut_snps[$alleles_snps[0] - 1])) {
					($ref_pos[$i], $strain_pos[$i]) = &shift_pos($snps[3], $mut_snps[$alleles_snps[0] - 1], $ref_pos[$i], $strain_pos[$i], $fileHandlesShift[$i], $snps[1]);
				}
				$snps = read_file_line($f1);
			}
		}
	}
}
#			&merge($snps);




while($snps) {
	$snps = read_file_line($f1);
	chomp $snps;
	if(substr($snps, 0, 1) eq "#") {
		&check_file_existance($snps);
	} else {
		if($header_exists == 0) {
			print STDERR "This is not a vcf file with header! Can not proceed with file parsing!\n";
			exit;
		} else {
			@split = split('\t', $snps);
			if($chr ne $split[0]) {
				print STDERR "Preparing chromosome " . $split[0] . "\n";
				$chr = $split[0];
				if($outfile_open == 1) {
					&close_filehandles();
				}
				#open all files and set positions
				&open_filehandles();
			}
			for(my $i = 0; $i < @split - 10; $i++) {
				#Look at the strains and split according to alleles
				@alleles = split('/', $split[$i+10]);		
								if($homo == 0) {
					if($alleles[1] ne "0") {
						$fileHandlesMutation2[$i]->print($split[1] . "\t" . $split[3] . "\t" . $muts[$alleles[1] - 1] . "\n");
						if(length($split[3]) != length($muts[$alleles[1] - 1])) {
							while($ref_pos2[$i] < $split[1]) {
								$ref_pos2[$i]++;
								$strain_pos2[$i]++;
					#			$fileHandlesShift2[$i]->print($ref_pos2[$i] . "\t" . $strain_pos2[$i] . "\n");
							}
							($ref_pos2[$i], $strain_pos2[$i]) = &shift_pos($split[3], $muts[$alleles[1] - 1], $ref_pos2[$i], $strain_pos2[$i], $fileHandlesShift2[$i]);	
						}
					}
				}
			}
		}
	}
}

while($indels) {
	&merge($indels);
	$indels = read_file_line($f2);
}



} else {
	$snp = $mut_files[0];
	open($f1, $snp);
	$indel = $mut_files[1];
	open($f2, $indel);
	$snps = read_file_file($f1);
	$indels = read_file_line($f2);
	while(substr($indels, 0, 1) eq "#") {
		$indels = read_file_line($f2);
		next;
	}
        while($snps and $indels) {
                @s = split('\t', $snps);
                $s_c = &convert_chr_to_num($s[0]);
                @i = split('\t', $indels);
                $i_c = &convert_chr_to_num($i[0]);
                if($s_c > $i_c) {
                        if($change == 0) {
                                print STDERR "\tchromosome " . $s[0] . "\n";
                        }
                        &merge($indels);
                        $indels = read_file_line($f2);
                        $change = 1;
                }
                if($i_c > $s_c) {
                        if($change == 0) {
                                print STDERR "\tchromosome " . $i[0] . "\n";
                        }
                        &merge($snps);
                        $snps = read_file_line($f1);
                        $change = 1;
                }

                if($s_c == $i_c) {
                        $change = 0;
                        if($s[1] < $i[1]) {
                                &merge($snps);
                                $snps = read_file_line($f1);
                        } else {
                                &merge($indels);
                                $indels = read_file_line($f2);
                        }
                }
        }
                &merge($snps);
                $snps = read_file_line($f1);
        }
        while($indels) {
                &merge($indels);
                $indels = read_file_line($f2);
        }

	@i = split('\t', $indels);
}

sub open_filehandles{
	for(my $i = 0; $i < @split - 10; $i++) {
		$filename = $data . "/" . $header[$i+10] . "/chr" . $chr . "_allele_1.mut";
		open my $fh_mut, ">", "$filename" or die "Can't open $filename: $!\n";
		push @fileHandlesMutation, $fh_mut;
		$filename = $data . "/" . $header[$i+10] . "/chr" . $chr . "_allele_1.shift";
		open my $fh_shift, ">", "$filename" or die "Can't open $filename: $!\n";
		push @fileHandlesShift, $fh_shift;
		$ref_pos[$i] = 0;
		$strain_pos[$i] = 0;
		if($homo == 0) {
			$filename = $data . "/" . $header[$i+10] . "/chr" . $chr . "_allele_2.mut";
			open my $fh_mut2, ">", "$filename" or die "Can't open $filename: $!\n";
			push @fileHandlesMutation2, $fh_mut2;
			$filename = $data . "/" . $header[$i+10] . "/chr" . $chr . "_allele_2.shift";
			open my $fh_shift2, ">", "$filename" or die "Can't open $filename: $!\n";
			push @fileHandlesShift2, $fh_shift2;
			$ref_pos2[$i] = 0;
			$strain_pos2[$i] = 0;
		}
	}	
	$outfile_open = 1;		
}

sub close_filehandles{
	for(my $i = 0; $i < @split - 10; $i++) {
		close $fileHandlesMutation[$i];
		close $fileHandlesShift[$i];
		if($homo == 0) {
			close $fileHandlesMutation2[$i];
			close $fileHandlesShift2[$i];
		}
	}
	undef @fileHandlesMutation;
	undef @fileHandlesShift;
	if($homo == 0) {
		undef @fileHandlesMutation2;	
		undef @fileHandlesShift2;
	}
}

sub check_file_existance{
	my $line = $_[0];
	#This is the header
	#Save that a header exists
	$header_exists = 1;
	@header = split('\t', $line);
	#Now go through all strains and see if the folder already exists
	for(my $i = 10; $i < @header; $i++) {
		print $header[$i] . "\n";
		if(-e $data . "/" . $header[$i]) {
			print STDERR "Folder for " . $header[$i] . " already exists!\n";
			if($force == 1) {
				print STDERR "Deleting folder " . $header[$i] . "!\n";
				print STDERR "Waiting for 10 seconds\n";
				print STDERR "Press Ctrl + C to interrupt\n";
				for(my $j = 0; $j < 10; $j++) {
					print STDERR ".";
				#	sleep(1);
				}
				print STDERR "\n";
				print STDERR "ADD sleep in again!\n";
				`rm -rf $data/$header[$i]/*`;
			} else {
				print STDERR "If you want to overwrite this folder use -force!\n";
				exit;
			}
		} else {
			`mkdir $data/$header[$i]`;
		}
	} 
}

sub shift_pos{
	my $ref = $_[0];
	my $strain = $_[1];
	my $pos_ref = $_[2];
	my $pos_strain = $_[3];
	my $fileHandle = $_[4];
	my $mut_pos = $_[5];
	while($ref_pos[$i] < $mut_pos) {
		$ref_pos[$i]++;
		$strain_pos[$i]++;
	}	

	if(length($ref) > length($strain)) {
		if(length($strain) > 1) {
			#we are already at the firs tposition of this deletion, so just add length(strain) - 1
			for(my $i = 1; $i < length($strain); $i++) {
				$pos_ref++;
				$pos_strain++;
				$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
			}
			#Next we jump in the reference sequence to the end of the insertion in the ref (deletion in strain), therefore we have to add 1 (because we look at the next position) and the length of the difference
			$pos_ref++;
			$pos_strain++;
			$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
			$pos_ref = $pos_ref + length($ref) - length($strain);
		} else {
			$pos_ref++;
			$pos_strain++;
			$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
			$pos_ref = $pos_ref + length($ref) - length($strain);
		}
	#Reference is shorter than mutated strain
	} else {
		if(length($ref) > 1) {
			#Add the overlap. We are already at the first position of the insertion, so we just add -1 less
			for(my $i = 0; $i < length($ref); $i++) {
				$pos_strain++;
				$pos_ref++;
				$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
			}
			#Add the positions of the query and do not change the position of the reference
			for(my $i = 0; $i < length($strain) - length($ref); $i++) {
				$pos_strain++;
				$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n"); 
			}
		} else {
			#We are already at the current position of the mutation
			$pos_strain++;
			$pos_ref++;
			$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
			for(my $i = 0; $i < length($strain) - length($ref); $i++) {
				$pos_strain++;
				$fileHandle->print($pos_ref . "\t" . $pos_strain . "\n");
			}
		}
	}
	return ($pos_ref, $pos_strain);
}


sub merge{
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
#       $current[4] =~ s/^\./$current[3]/g;
        $current[4] =~ s/\.//g;
        my @variants = split(',', $current[4]);

        unshift @variants, $current[3];

        if($filter == 1 && $current[6] ne "PASS") {
                return 1;
        }

        my $var;
        for(my $i = 9; $i < @current; $i++) {
                @all = split('/', ((split(':', $current[$i]))[0]));
                if($all[0] eq ".") {
                        $all[0] = 0;
                }
                if(@all == 1) {
                        $all[1] = 0;
                } else {
                        if($all[1] eq ".") {
                                $all[1] = 0;
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
                                        #Remove the first bases if they are the same
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
                                @last = split(',', $lines->[$h]->[$i-9]->[-1]) if @{$lines->[$h]->[$i-9]} > 0;
                                if(@{$lines->[$h]->[$i-9]} < 1 || $last[0] ne $current[0]) {
                                        $merge_line = $current[0] . "," . $pos . "," . $ref . "," . $var;
                                        push(@{ $lines->[$h]->[$i-9] }, $merge_line);
                                        next;
                                }
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
                                if($last[1] == $pos) { next; }
                                if($current[1] < $last[1] + length($last[2])) {
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
                                        next;
                                }
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
                                $merge_line = $current[0] . "," . $pos . "," . $ref . "," . $var;
                                push(@{ $lines->[$h]->[$i-9]}, $merge_line);
                        }
                }
        }
}


sub read_file_line {
        my $fh = shift;
        my @split;

        if ($fh and my $line = <$fh>) {
                chomp $line;
#               return [ split('\t', $line) ];
                return $line;
        }
        return;
}
