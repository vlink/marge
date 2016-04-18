#!/usr/bin/perl

package textfiles_module;
use strict;
use Getopt::Long;
use Storable;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
use Set::IntervalTree;

use threads;
use Data::Dumper;

$_ = "" for my($snp, $indel, @chr, $current_chr, $data, $filename, $out, $genome);
$_ = 0 for my($filter, $same, $help, $homo, $length_mut, $length_ref, $force, $lines_f1, $lines_f2, $lines_all, $line_count, $num_strains, $header_exists, $outfile_open, $last_h, $run_through, $none_number_chromosome, $core);
$_ = () for my($lines, @merge_line, $header, @split, $snps, $indels, $f1, $f2, %strains_to_use, @mut_files, @header, %last_shift_pos_ref, %last_shift_pos_strain, @strains_to_use, @s, @i, $check_snp, $check_indel, $s_c, $i_c, @fileHandlesMutation, @fileHandlesMutation2, @fileHandlesShift, @fileHandlesShift2, @ref_pos, @strain_pos, @ref_pos2, @strain_pos2, @pos_ref, @pos_strain, %lookup_no_number, %lookup_number, %multi_strains);
#Config with default data folder

#$none_number_chromosome = 1000;

#SUBFUNCTIONS
#Routine to start threading

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
		print $pos_ref . "\t" . $pos_strain . "\n";
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
	my $data = $_[4];
	if(exists $lookup_number{$chr}) {
		$chr = $lookup_number{$chr};
	}
        my $last_shift = 0;
	$_ = () for my (@array_shift, @tmp_split);
	my $start = 0;
	$filename = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".mut";
	print $filename . "\n";
	my $out = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".strain_to_ref.vector";
	my $diff = 0;
	if(-e $filename) {
		open FH, "<$filename";
		open OUT, ">$out";
		print $filename . "\n";
		my $strain_pos = 0;
		my $ref_pos = 0;
		my $last_strain = 0;
		my $last_ref = 0;
		my $shift = 0;
		my $last_shift = 0;
		foreach my $line (<FH>) {
			chomp $line;
			@tmp_split = split('\t', $line);
		#	print $line . "\n";
			while($ref_pos < $tmp_split[0]) {
				$strain_pos++;
				$ref_pos++;
			}
		#	#tmp_split[0] is reference position - tmp_split[1] is strains position - calculate shifting vector out of it
		#	if($start == $tmp_split[0] - 1) {
		#		next;
		#	}	
			if(length($tmp_split[1]) == length($tmp_split[2])) { next; }
			#Insertion in strain
			#Strain length > Ref length
			if(length($tmp_split[2]) > length($tmp_split[1])) {
			#	print $line . "\n";
			#	print "split1: " . length($tmp_split[1]) . "\n";
			#	print "split2: " . length($tmp_split[2]) . "\n";
			#	print "strain: " . $strain_pos . "\t" . "ref: " . $ref_pos . "\n";
			#	print "TEST\n";
			#	print "difference: " . (length($tmp_split[1]) - length($tmp_split[2])) . "\n";
				$diff = (length($tmp_split[1]) - length($tmp_split[2]));
				$shift = $shift + $diff;
			#	print "shift: " . $shift . "\n";
				$strain_pos = $strain_pos + abs($diff);
			#	print "new strain pos: " . $strain_pos . "\n";
				print $last_shift . "\t" . $last_strain . "\t" . $strain_pos . "\n";
				$last_shift = $shift;
				$strain_pos++;
				$ref_pos++;
				$last_strain = $strain_pos;
			} else {
			#	print $line . "\n";
			#	print "del split1: " . length($tmp_split[1]) . "\n";
			#	print "del split2: " . length($tmp_split[2]) . "\n";
			#	print "" . (length($tmp_split[1]) - length($tmp_split[2])) . "\n";
				$ref_pos = $ref_pos + length($tmp_split[1]) - 1;	
			#	print "strain: " . $strain_pos . "\tref: " . $ref_pos . "\n";
				$diff = (length($tmp_split[1]) - length($tmp_split[2]));
				$shift = $shift + $diff;
			#	print "shift: " . $shift . "\n";
			#	print $last_shift . "\t" . $last_strain . "\t" . $strain_pos . "\n";
				$last_shift = $shift;
				$strain_pos++;
				$ref_pos++;
				$last_strain = $strain_pos;
			}
		#	print OUT $last_shift . "\t" . $start . "\t" . ($tmp_split[0] + length($tmp_split[2])) . "\n";
			print $last_shift . "\t" . $start . "\t" . ($tmp_split[0] + length($tmp_split[2]) - 1) . "\n";
			$start = $tmp_split[0] + length($tmp_split[2]);
		#	$last_shift = $last_shift + (length($tmp_split[1]) - length($tmp_split[2]));
		#	$array_shift[$start] = $last_shift;
		#	$start++;
		#	$last_shift_pos_strain{$strain}{$chr}{$a}{'pos'} = $start;
		#	$last_shift_pos_strain{$strain}{$chr}{$a}{'shift'} = $last_shift;
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
1;
