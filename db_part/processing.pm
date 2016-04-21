#!/usr/bin/perl -w

package processing;

sub create_genome {
	my $chr = $_[0];
	my %strains = %{$_[1]};
	my $input = $_[2];
	my $output_dir = $_[3];
	my $reference_chr = $_[4];
	my $allele = $_[5];
	print STDERR "Generating chromosome " . $chr . " for allele " . $allele . "\n";
	open FH, "<$reference_chr" or die "Can't open: " . $reference_chr . "\n";
	#First save the genome in an array
	$_ = () for my (@cache, @save_cache, @split, @split_mut);
	my $strain_seq;
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			next;
		}
		@split = split('', $line);
		foreach my $s (@split) {
			push(@cache, $s);
		}
	}
	#Save reference genome array, so when there are several strains we can reset the genome for every strain
	@save_cache = @cache;
	foreach my $strains (keys %strains) {
		$mutation_file = $input . "/" . $strains . "/chr" . $chr . "_allele_" . $allele . ".mut";
		$output = $output_dir . "/" . $strains . "/chr" . $chr . "_allele_" . $allele . ".fa";
		#Now open the mutation file and change the mutated sequences
		if(!(-e $mutation_file)) {
			print STDERR $mutation_file . " does not exist - go to next strain\n";
			next;
		}
		if(!(-e $output_dir . "/" . $strains)) {
			`mkdir -p $output_dir/$strains`;
		}
		open FH, "<$mutation_file" or die "Can't open $mutation_file\n";
		#Add SNPs and InDels into reference genome array
		foreach my $line (<FH>) {
			chomp $line;
			@split_mut = split('\t', $line);
			if($split_mut[0] > @cache) { last; }
			if(length($split_mut[1]) > 1) {
				for(my $i = 1; $i < length($split_mut[1]); $i++) {
					$cache[$split_mut[0] + $i - 1] = "";
				}
			} else {
				$cache[$split_mut[0] - 1] = $split_mut[2];
			}
		}
		close FH;

		#Write strain genome output - use 50bg per line (default HOMER uses and needed for grabbing sequenes with byte offsets later on)
		print STDERR "\t\t" . $strains . "\n";
		$strain_seq = join("", @cache);
		open OUT, ">$output" or die "Can't open $output\n";
		print OUT ">chr" . $chr . "\n";
		while(length($strain_seq) > 50) {
			print OUT substr($strain_seq, 0, 50) . "\n";
			$strain_seq =~ s/^.{50}//s; 
		}
		print OUT $strain_seq . "\n";
		close OUT;
		@cache = @save_cache;
	}
}

sub create_offset_ref_to_strain {
        my $chr = $_[0];
	my $a = $_[1];
	my $strain = $_[2];
	my %last_shift_pos_ref = %{$_[3]};
	my $data = $_[4];
	if(exists $lookup_number{$chr}) {
		$chr = $lookup_number{$chr};
	}
        $filename = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".mut";
	my $out = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".ref_to_strain.vector";
	$_ = () for my(@array_shift, @tmp_split, @a);
	$_ = 0 for my ($run, $last_shift);
        if(-e $filename) {
                open FH, "<$filename";
		open OUT, ">$out";
                foreach my $line (<FH>) {
                        chomp $line;
                        @tmp_split = split('\t', $line);
			if(length($tmp_split[1]) == length($tmp_split[2])) { next; }
			#Print last shift vector for position where run variable is currently and where the next mutation starts
			print OUT $last_shift . "\t" . $run . "\t" . $tmp_split[0]. "\n";
			#Change last shift to the current shift 
			$last_shift = $last_shift + (length($tmp_split[2]) - length($tmp_split[1]));
			$run = $tmp_split[0] + 1;
			#Save last shift 
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
	$_ = 0 for my ($start, $diff);
	$filename = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".mut";
	my $out = $data . "/" . $strain . "/chr" . $chr . "_allele_" . $a . ".strain_to_ref.vector";
	if(-e $filename) {
		open FH, "<$filename";
		open OUT, ">$out";
		my $strain_pos = 0;
		my $ref_pos = 0;
		my $last_strain = 0;
		my $last_ref = 0;
		my $shift = 0;
		my $last_shift = 0;
		#Now we need to keep track of the current position in the strain - mutations are annotated in reference - so we need two variables to keep track of the current places in strain and reference
		foreach my $line (<FH>) {
			chomp $line;
			@tmp_split = split('\t', $line);
			while($ref_pos < $tmp_split[0]) {
				$strain_pos++;
				$ref_pos++;
			}
		#	#tmp_split[0] is reference position - tmp_split[1] is strains position - calculate shifting vector out of it
			if(length($tmp_split[1]) == length($tmp_split[2])) { next; }
			#Insertion in strain
			#Strain length > Ref length
			if(length($tmp_split[2]) > length($tmp_split[1])) {
				$diff = (length($tmp_split[1]) - length($tmp_split[2]));
				$shift = $shift + $diff;
				$strain_pos = $strain_pos + abs($diff);
			#Deletion in strain
			} else {
				$ref_pos = $ref_pos + length($tmp_split[1]) - 1;	
				$diff = (length($tmp_split[1]) - length($tmp_split[2]));
				$shift = $shift + $diff;
			}
			print OUT $last_shift . "\t" . $last_strain . "\t" . $strain_pos . "\n";
			$last_shift = $shift;
			$strain_pos++;
			$ref_pos++;
			$last_strain = $strain_pos;
		}
		close FH;
		close OUT;
	}
	return \%last_shift_pos_strain;
}

sub get_refseq_for_gene{
	my $gene_file = $_[0];
	my $gene = $_[1];
        open FH, "<$gene_file";
	$_ = "" for my($transcript, $refseq_NM);
        my $gene_found = 0;
        foreach my $line (<FH>) {
                chomp $line;
                @split = split('\t', $line);
                for(my $i = 0; $i < @split; $i++) {
                        if(uc($split[$i]) eq uc($gene)) {
                                $gene_found = 1;
                        }
                        if(substr($split[$i], 0, 2) eq "NM") {
                                $refseq_NM = $split[$i];
                        }
                }
                if($gene_found == 1) {
                        $transcript = $refseq_NM;
                        last;
                }
        }
        close FH;
	return $refseq_NM;
}


sub save_transcript {
	my $refseq_file = $_[0];
	my $transcript = $_[1];
        $_ = "" for my ($stop, $chr, $strand);
	$_ = () for my (@introns, @parts, %peaks, %exons);
	$_ = 0 for my ($start, $exon);
        open FH, "<$refseq_file";
        foreach my $line (<FH>) {
                chomp $line;
                @split = split('\t', $line);
                if($split[0] eq $transcript) {
                        $chr = $split[1];
                        $strand = $split[4];
                        @introns = split(",", $split[5]);
                        if($split[4] eq "-") {
                                for(my $i = @introns - 1; $i > 0; $i--) {
                                        @parts = split(":", $introns[$i]);
                                        if(length($parts[0]) > 4) {
                                                next;
                                        }
                                        if(substr($parts[0], 0, 1) eq "E") {
                                                $start = $parts[1] - 1;
                                                $i--;
                                                @parts = split(":", $introns[$i]);
                                                if(length($parts[0]) < 6 && substr($parts[0], 0, 1) eq "I") {
                                                        print STDERR "Weird annotation!\n";
                                                }
                                                $stop = $parts[1] - 1;
                                                $peaks{substr($split[1], 3)}{$start} = $stop;
                                                $exons{$start . "_" . $stop} = $exon;
                                                $exon++;
                                        }
                                }
                        } else {
                                for(my $i = 0; $i < @introns; $i++) {
                                        @parts = split(':', $introns[$i]);
                                        if(length($parts[0]) > 4) {
                                                next;
                                        }
                                        if(substr($parts[0], 0, 1) eq "E") {
                                                $start = $parts[1] - 1;
                                                $i++;
                                                if(@introns > $i) {
                                                        @parts = split(':', $introns[$i]);
                                                        if(length($parts[0]) < 4 && substr($parts[0], 0, 1) ne "I") {
                                                                print STDERR "Weird pos strand\n";
                                                        }
                                                        $stop = $parts[1] - 1;
                                                } else {
                                                        $stop = $split[3];
                                                }
                                                $peaks{substr($split[1], 3)}{$start} = $stop;
                                                $exons{$start . "_" . $stop} = $exon;
                                                $exon++;
                                        }
                                }
                        }
                }
        }
        return ($strand, $chr, \%peaks, \%exons);
}

sub check_homer_files{
	my $refseq_file = $_[0];
	my $gene_file = $_[1];
	my $config = $_[2];
	my $gene_number = $_[3];
	my $genome = $_[4];
	my $path;
	my $species;
        if($refseq_file eq "" && exists $config->{'homer_path'} && $config->{'homer_path'} ne "") {
                #Homer exists now check if the genome and the annotation files exist
                $path = $config->{'homer_path'} . "/data/genomes/" . $genome . "/" . $genome . ".rna";
                if(!(-e $path)) {
                        print STDERR "File $path is missing!\n";
                        exit;
                }
                $refseq_file = $path;
                if($gene_number > 0) {
                        #First find out which species this genome belongs to - homer config
                        $path = $config->{'homer_path'} . "/config.txt";
                        open FH, "<$path" or die "Can not find $path!\n";
                        foreach my $line (<FH>) {
                                chomp $line;
                                @split = split('\t', $line);
                                if($split[0] eq $genome) {
                                        my @a = split(",", $split[-1]);
                                        $species = $a[0];
                                }
                        }
                        if($species eq "") {
                                print STDERR "Could not figure out which species this genome belongs to - please specify it!\n";
                                exit;
                        }
                        $path = $config->{'homer_path'} . "/data/accession/" . lc($species) . "2gene.tsv";
                        if(!(-e $path)) {
                                print STDERR "File $path is missing!\n";
                                exit;
                        }
                        $gene_file = $path;
                }
        }
        if($refseq_file eq "") {
                print STDERR "No RefSeq file specified and RefSeq File could not be found!\n";
                exit;
        }
        if($gene_number > 0 && $gene_file eq "") {
                print STDERR "No gene file specified and Gene file cound not be found!\n";
                exit;
        }
	return($refseq_file, $gene_file);
}
1;
