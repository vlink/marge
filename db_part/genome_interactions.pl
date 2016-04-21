#!/usr/bin/perl -w

use strict;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use Getopt::Long;
use config;
use general;
use Data::Dumper;
my $config = config::read_config();
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
use processing;
use analysis_tree;
use Set::IntervalTree;


sub printCMD {
	print STDERR "Usage\n";
        print STDERR "\t-genome: Genome\n";
	print STDERR "\t-species: species of genome e.g. mm10 is mouse\n";
	print STDERR "\t-refseq_file: File with RefSeq IDs (check HOMER format)\n";
	print STDERR "\t-gene_file: File with Gene IDs (check HOMER format)\n";
        print STDERR "\t-method <gene|protein|both|genomic|align_gene|align_protein|align_both|align_genomic_gene|align_genomic_protein|align_genomic_both|genomic_protein|genomic_both>\n";
        print STDERR "\t-gene <gene name> (Comma separated list)\n";
        print STDERR "\t-transcript <RefSeq transcript name> (Comma separated list)\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - if not defined only the reference is used - if all is specified all strains are used\n";
        print STDERR "\t-start <pos>: Start position of genomic region (Comma separated list - same length as end, chr, and strand)\n";
        print STDERR "\t-end <pos>: Stop position of genomic region (Comma separated list - same length as start, chr, and strand)\n";
        print STDERR "\t-chr <chromosome>:Chromosome for genomic region (Comma separated list - same length as start, end, and strand)\n";
        print STDERR "\t-strand <+|->: Default: + (Comma separated list - same length as start, end, and chr)\n";
	print STDERR "\t-data: Path to data directory for strains - default specified in config\n";
        print STDERR "\t-hetero: Data is heterozygous\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

$_ = "" for my ($genome, $path, $species, $refseq_file, $gene_file, $data, $refseq_NM, $method, $strand, $chr, $gene, $transcript);
$_ = 0 for my ($hetero, $html, $exon, $gene_found, $allele, $shift, $longest_seq, $filter_no_mut, $start);
$_ = () for my (%strains, @split, @strains, @strains2, $seqs, @parts, @introns, %exons, %peaks, %gene, %codons, %tree, %last_strain, %lookup_strain, $f1, $mut_line, $f1_file, @fileHandles, @exons, %align_nt, %align_nt_seq, @gene, @transcript, @list, @list_id, @start, @end, @chr, @strand);


GetOptions(     "genome=s" => \$genome,
		"species=s" => \$species,
		"refseq_file=s" => \$refseq_file,
		"gene_file=s" => \$gene_file,
                "method=s" => \$method,
                "gene=s{,}" => \@gene,
                "transcript=s{,}" => \@transcript,
                "strains=s{,}" => \@strains,
                "start=s{,}" => \@start,
                "end=s{,}" => \@end,
                "chr=s{,}"=> \@chr,
                "hetero" => \$hetero,
		"data=s" => \$data,
                "html" => \$html)
        or die("Error in command line options!\n");

if($hetero == 1) {
	$allele = 2;
} else {
	$allele = 1;
}

for(my $i = 0; $i < @gene; $i++) {
	$gene[$i] =~ s/,//g;
	push(@list, $gene[$i]);
	push(@list_id, "gene");
}
for(my $i = 0; $i < @transcript; $i++) {
	$transcript[$i] =~ s/,//g;
	push(@list, $transcript[$i]);
	push(@list_id, "refseq");
}
if(@start != @end || @start != @chr) {
	print STDERR "Length of start, end, chr, and strand are different\n";
	exit;
}

if(@strand < @start) {
	for(my $i = 0; $i < @start - @strand; $i++) {
		$strand[$i] = "+";
	}
}

for(my $i = 0; $i < @start; $i++) {
	$start[$i] =~ s/,//g;
	$end[$i] =~ s/,//g;
	$chr[$i] =~ s/,//g;
	$strand[$i] =~ s/,//g;
	if($i == 0) {
		$shift = @list;
	}
	push(@list, $start[$i]);
	push(@list_id, "other");
}

#Need to check if HOMER is installed so we can access the gene file or alternativel if a file is specified
if(@gene > 0 || @transcript > 0) {
	&check_homer_files();
}
&write_codons();

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

if($data eq "") {
	$data = $config->{'data_folder'};
}

print STDERR "Loading shift vectors\n";
for(my $i = 0; $i < @strains; $i++) {
        my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data, $allele, "ref_to_strain");
        $tree{$strains[$i]} = $tree_ref;
        $lookup_strain{$strains[$i]} = $lookup;
        $last_strain{$strains[$i]} = $last;
}

for(my $entry_id = 0; $entry_id < @list; $entry_id++) {
	%peaks = ();
	%exons = ();
	#Create hash with sequences we want to get
	if($list_id[$entry_id] eq "gene") {
		$gene = $list[$entry_id];
		print $gene . "\n";
		&get_refseq_for_gene();
		if($transcript eq "") {
			print STDERR "Could not find Ref-Seq entry for " . $gene . "\n";
			print STDERR "Skip\n";
			next;
		}
		($strand, $chr) = &save_transcript();
	} elsif($list_id[$entry_id] eq "refseq") {
		$transcript = $list[$entry_id];
		print $transcript . "\n";
		($strand, $chr) = &save_transcript();
	} else {
		if(substr($chr[$entry_id - $shift], 0, 3) ne "chr") {
			$chr[$entry_id - $shift] = "chr" . $chr[$entry_id - $shift];
		}
		$strand = $strand[$entry_id - $shift];
		$chr = $chr[$entry_id - $shift];
		print $chr . "\t" . $start[$entry_id - $shift] . "\t" . $end[$entry_id - $shift] . "\t" . $strand . "\n";
		$peaks{substr($chr[$entry_id - $shift], 3)}{$start[$entry_id - $shift]} = $end[$entry_id - $shift];
		$exons{$start[$entry_id - $shift] . "_" . $end[$entry_id - $shift]} = 1;
	}
	$seqs = ();
	if(-e "tmp") { `rm tmp`; }
	($seqs, $longest_seq, $filter_no_mut) = analysis::get_seq_for_peaks("tmp", \%peaks, \@strains, $data, $allele, $exon, 0, 0, \%tree, \%lookup_strain, \%last_strain);
	if($method eq "gene" || $method eq "genomic") {
		&assemble_gene(1);
	} elsif($method eq "protein" || $method eq "genomic_protein") {
		&assemble_gene(0);
		foreach my $s (keys %gene) {
			print $s . "\n";
			if($strand eq "-") {
				$gene{$s} = analysis::rev_comp($gene{$s});
			}
			&translate_seq($gene{$s}, 1);
		}
	} elsif($method eq "both" || $method eq "genomic_both") {
		&assemble_gene(1);
		foreach my $s (keys %gene) {
			print $s . "\n";
			if($strand eq "-") {
				$gene{$s} = analysis::rev_comp($gene{$s});
			}
			&translate_seq($gene{$s}, 1);
		}
	} elsif($method eq "align_gene" || $method eq "align_genomic_gene") {
		&assemble_gene(0);
		&assemble_alignment(1);
	} elsif($method eq "align_protein" || $method eq "align_genomic_protein") {
		&assemble_gene(0);
		&assemble_alignment(0);
		&align_protein();
	} elsif($method eq "align_both" || $method eq "align_genomic_both") {
		&assemble_gene(0);
		&assemble_alignment(1);
		&align_protein();
	}
}
sub align_protein {
	my $local_seq = "";
	my $local_seq_align = "";
	my $protein_seq = "";
	my $align_prot = "";
	my $num_gaps = 0;
	my $num_gaps_tmp = 0;
	my $removed_gaps = 0;
	my $prot = "";
	my $x = "-";
	my $prot_pos = 0;
	my %prot_seq;
	my %align_seq;
	my %align_prot_seq;
	for(my $i = 0; $i < @strains; $i++) {
		$local_seq = $align_nt_seq{$i};
		$local_seq =~ s/-//g;
		$prot_seq{$i} = &translate_seq($local_seq, 0);
	}
	for(my $i = 0; $i < @strains; $i++) {
		$prot = "";
		$num_gaps = 0;
		$num_gaps_tmp = 0;
		$align_prot = "";
		$prot_pos = 0;
		$local_seq = $align_nt_seq{$i};
		$local_seq_align = $align_nt{$i};	
	#	print $local_seq . "\n";
		$local_seq =~ s/-//g;
		while(length($local_seq_align) >= 3) {
			$num_gaps_tmp =  () = substr($local_seq_align, 0, 3) =~ /$x/g; 
			$num_gaps += $num_gaps_tmp;
			$removed_gaps = 0;
	#		print substr($local_seq_align, 0, 3) . "\n";
			if(index(substr($local_seq_align, 0, 3), "-") != -1) {
	#			print "num or gaps: " . $num_gaps . "\n";
				if($num_gaps >= 3) {
					$prot .= "-";
					$align_prot .= "-";
					$num_gaps = $num_gaps - 3;
					$removed_gaps = 1;
				} else {
					$prot .= substr($prot_seq{$i}, $prot_pos, 1);
					$align_prot .= "|";
					$prot_pos++;

				}
			} elsif(index(substr($local_seq_align, 0, 3), '.') != -1) {
				$prot .= substr($prot_seq{$i}, $prot_pos, 1);
				$align_prot .= ".";
				$prot_pos++;
			} else {
	#			print "everything is fine\n";
				$prot .= substr($prot_seq{$i}, $prot_pos, 1);
				$align_prot .= "|";
				$prot_pos++;
			}
			if(length($local_seq_align) <= 3) {
				$local_seq_align = "";
				$local_seq = ""; 
			} else {
				if($removed_gaps == 0) {
					$local_seq = substr($local_seq, 3);
				}
				$local_seq_align = substr($local_seq_align, 3);
			}
		}
		$align_prot_seq{$i} = $prot;
		$align_seq{$i} = $align_prot;
	}
	#Check for synoymous mutations
	my @alignment_seq;
	my $diff_as = 0;
	my $current_pos = 0;
	for(my $i = 0; $i < @strains; $i++) {
		@alignment_seq = split('', $align_seq{$i});
		$current_pos = 0;
		foreach my $c (@alignment_seq) {
			#Check aa of all strains
			if($c eq ".") {
				$diff_as = 0;
				for(my $k = 0; $k < @strains; $k++) {
					if(substr($align_prot_seq{$i}, $current_pos, 1) ne substr($align_prot_seq{$k}, $current_pos, 1)) {
						$diff_as = 1;
					}
				}
				if($diff_as == 0) {
					$align_seq{$i} = substr($align_seq{$i}, 0, $current_pos) . "*" . substr($align_seq{$i}, $current_pos + 1);
				}
			}
			$current_pos++;
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		print $strains[$i] . "\n";
		print $align_prot_seq{$i} . "\n";
		print $align_seq{$i} . "\n";
	}
}

sub translate_seq{
	my $local_seq = $_[0];
	my $print = $_[1];
	my $protein_seq = "";
	while(length($local_seq) >= 3) {
		$protein_seq .= $codons{substr($local_seq, 0, 3)};
		$local_seq = substr($local_seq, 3);
	}
	if($print == 1) {
		print $protein_seq . "\n";
	}
	return $protein_seq;
}

sub assemble_alignment{
	&open_filehandles();	
	my $print = $_[0];
	my @mutation_array;
	my %mutation_hash;
	my $run_index = 0;
	my $complete_index;
	my $last_index = 0;
	my @align;
	my @align_seq;
	my @orig_seq;
	my @mut_shift;
	for(my $i = 0; $i < @strains; $i++) {
		@{$orig_seq[$i]} = split('', $gene{uc($strains[$i])});
		$mut_shift[$i] = 0;
	}
	my $current_mut;
	my $strain_mut;
	my $mut_exists = 0;
	for(my $i = 0; $i < @strains; $i++) {
		$run_index = 0;
		$complete_index = 0;
		if(!defined $fileHandles[$i]) { next; }
		$mut_line = read_file_line($fileHandles[$i]);
		@split = split('\t', $mut_line);
		foreach my $e (sort {$a cmp $b} keys %exons) {
			$mut_exists = 0;
			@exons = split("_", $e);
			while($split[0] < $exons[0]) {
				$mut_line = read_file_line($fileHandles[$i]);
				@split = split('\t', $mut_line);
			}
			while($split[0] < $exons[1]) {
				$run_index = $complete_index + ($split[0] - $exons[0]) - 1;
				$mutation_hash{$i}{$run_index} = length($split[1]) - length($split[2]);
				$mut_line = read_file_line($fileHandles[$i]);
				@split = split('\t', $mut_line);
			}
			$complete_index = $complete_index + $exons[1] - $exons[0];
		}
	}

	#Check if there are mutations in every strain that are always the same
	my $num_mut = 0;
	my $diff = 0;
	foreach my $gene (keys %gene) {
		if(length($gene{$gene}) > $last_index) {
			$last_index = length($gene{$gene});
		}
	}
	foreach my $strain (keys %mutation_hash) {
		foreach my $pos (keys %{$mutation_hash{$strain}}) {
			$num_mut = 0;
			$diff = 0;
			for(my $i = 0; $i < @strains; $i++) {
				if(exists $mutation_hash{$i}{$pos}) {
					$num_mut++;
					if($mutation_hash{$i}{$pos} ne $mutation_hash{$strain}{$pos}) {
						$diff = 1;
					}
				}
			}
			if($num_mut == @strains && $diff == 0) {
				for(my $i = 0; $i < @strains; $i++) {
					delete $mutation_hash{$i}{$pos};
				}
			}
		}
	}
	for(my $run = 0; $run < $last_index; $run++) {
		for(my $i = 0; $i < @strains; $i++) {
			if(!defined $orig_seq[$i][$run - $mut_shift[$i]]) { next; }
			if(exists $mutation_hash{$i}{$run}) {
				if($mutation_hash{$i}{$run} == 0) {
					$align[$i][$run] = ".";
					$align_seq[$i][$run] = $orig_seq[$i][$run - $mut_shift[$i]];
				} elsif($mutation_hash{$i}{$run} > 0) {
					if($strand eq "-") {
						$align[$i][$run] =  ( '-' x $mutation_hash{$i}{$run} ) . "|"; 
						$align_seq[$i][$run] =  ( '-' x $mutation_hash{$i}{$run} ) . $orig_seq[$i][$run - $mut_shift[$i]];
					} else {
						$align[$i][$run] = "|" . ( '-' x $mutation_hash{$i}{$run} ) . "|"; 
						$align_seq[$i][$run] = $orig_seq[$i][$run - $mut_shift[$i]] . ( '-' x $mutation_hash{$i}{$run} );
					}
					$mut_shift[$i] += $mutation_hash{$i}{$run};
				} else {
					#First check if there are several insertions and get the longest
					my $current = $mutation_hash{$i}{$run};
					for(my $k = 0; $k < @strains; $k++) {
						if(exists $mutation_hash{$k}{$run} && $mutation_hash{$k}{$run} < $mutation_hash{$i}{$run}) {
							$i = $k;
						}
					}
					for(my $k = 0; $k < @strains; $k++) {
						if($strand eq "-") {
							if(!exists $mutation_hash{$k}{$run}) {
								$align[$k][$run] = "|" . ( '-' x abs($mutation_hash{$i}{$run}) );
								$align_seq[$k][$run] = $orig_seq[$k][$run - $mut_shift[$k]] . ( '-' x abs($mutation_hash{$i}{$run} ));
							} elsif(exists $mutation_hash{$k}{$run} && $mutation_hash{$k}{$run} != $mutation_hash{$i}{$run}) {
								if($mutation_hash{$k}{$run} < 0) {
									$align[$k][$run] = ( '-' x (abs($mutation_hash{$i}{$run}) - abs($mutation_hash{$k}{$run})));	
									$align_seq[$k][$run] = ( '-' x (abs($mutation_hash{$i}{$run}) - abs($mutation_hash{$k}{$run})));
								} else {
									$align[$k][$run] =  ( '-' x (abs($mutation_hash{$i}{$run}) + $mutation_hash{$k}{$run}) ) . "|"; 
									$align_seq[$k][$run] =  ( '-' x (abs($mutation_hash{$i}{$run}) + $mutation_hash{$k}{$run}) ) . $orig_seq[$k][$run - $mut_shift[$k]];
									$mut_shift[$k] += $mutation_hash{$k}{$run};
								}	
							} else {
								$align[$k][$run] = "|";
								$align_seq[$k][$run] = $orig_seq[$i][$run - $mut_shift[$k]];
							}
						} else {
							if(!exists $mutation_hash{$k}{$run}) {
								$align[$k][$run] = ( '-' x abs($mutation_hash{$i}{$run}) ) . "|";
								$align_seq[$k][$run] = ( '-' x abs($mutation_hash{$i}{$run} )) . $orig_seq[$k][$run - $mut_shift[$k]];
							} elsif(exists $mutation_hash{$k}{$run} && $mutation_hash{$k}{$run} != $mutation_hash{$i}{$run}) {
								if($mutation_hash{$k}{$run} < 0) {
									$align[$k][$run] = ( '-' x (abs($mutation_hash{$i}{$run}) - abs($mutation_hash{$k}{$run})));	
									$align_seq[$k][$run] = ( '-' x (abs($mutation_hash{$i}{$run}) - abs($mutation_hash{$k}{$run})));
								} else {
									$align[$k][$run] =  "|" . ( '-' x (abs($mutation_hash{$i}{$run}) + $mutation_hash{$k}{$run}) ); 
									$align_seq[$k][$run] = $orig_seq[$k][$run - $mut_shift[$k]] . ( '-' x (abs($mutation_hash{$i}{$run}) + $mutation_hash{$k}{$run}) );
									$mut_shift[$k] += $mutation_hash{$k}{$run};
								}	
							} else {
								$align[$k][$run] = "|";
								$align_seq[$k][$run] = $orig_seq[$i][$run - $mut_shift[$k]];
							}

						}
					}
					#We are through all mutations
					for(my $k = 0; $k < @strains; $k++) {
						delete $mutation_hash{$k}{$run};
					}
					$i = @strains + 1;
				}
			} else {
				if(!defined $orig_seq[$i][$run]) {
					next;
				}
			#	$k .= $strains[$i] . "\t" . $orig_seq[$i][$run - $mut_shift[$i]] . "(" . ($run) . ")\t\t";
				$align_seq[$i][$run] = $orig_seq[$i][$run];
				$align[$i][$run] = "|";
			}
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		$align_nt_seq{$i} = "" . join("", @{$align_seq[$i]});
		$align_nt{$i} =  "" . join("", @{$align[$i]});
		if($strand eq "-") {
			$align_nt_seq{$i} = analysis::rev_comp($align_nt_seq{$i});
			$align_nt{$i} = reverse($align_nt{$i});
		}
		if($print == 1) {
			print $strains[$i] . "\n";
			print $align_nt_seq{$i} . "\n";
			print $align_nt{$i} . "\n";
		}
	}
}

sub open_filehandles{
	for(my $i = 0; $i < @strains; $i++) {
		my $filename = $data . "/" . uc($strains[$i]) . "/chr" . substr($chr, 3) . "_allele_1.mut";
		if(-e $filename) {
			open my $fh, "<", $filename or die "Can't open $filename: $!\n";
			$fileHandles[$i] = $fh;
		}
	}
}

sub assemble_gene{
	my $print = $_[0];
	my $assembly = "";
	for(my $i = 0; $i < @strains; $i++) {
		$assembly = "";
		foreach my $e (sort {$a cmp $b} keys %exons) {
			$assembly .= $seqs->{$chr . "_" . $e . "_" . $strains[$i]};
		}
		$gene{$strains[$i]} = $assembly;
		if($print == 1) {
			if($strand eq "-") {
				$assembly = analysis::rev_comp($assembly);
			}
			print $strains[$i] . "\n";
			print $assembly;
			print "\n";
		}
	}

}

sub get_refseq_for_gene{
	open FH, "<$gene_file";
	$transcript = "";
	$gene_found = 0;
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
			close FH;
			last;
		}
	}
	close FH;
}

sub save_transcript {
	my $id = $_[0];
	my $stop;
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
	return ($strand, $chr);
}

sub check_homer_files{
	if($refseq_file eq "" && exists $config->{'homer_path'} && $config->{'homer_path'} ne "") {
		#Homer exists now check if the genome and the annotation files exist
		$path = $config->{'homer_path'} . "/data/genomes/" . $genome . "/" . $genome . ".rna";
		if(!(-e $path)) {
			print STDERR "File $path is missing!\n";
			exit;
		}
		$refseq_file = $path;
		if(@gene > 0) {
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
	if(@gene > 0 && $gene_file eq "") {
		print STDERR "No gene file specified and Gene file cound not be found!\n";
		exit;
	}
}

sub read_file_line {
        my $fh = shift;
        if ($fh and my $line = <$fh>) {
                chomp $line;
                return $line;
        }
        return;
}

sub write_codons{
 (%codons) = ('TCA'=>'S', #Serine
                'TCC'=>'S', #Serine
                'TCG'=>'S', #Serine
                'TCT'=>'S', #Serine
                'TTC'=>'F', #Phenylalanine
                'TTT'=>'F', #Phenylalanine
                'TTA'=>'L', #Leucine
                'TTG'=>'L', #Leucine
                'TAC'=>'Y', #Tyrosine
                'TAT'=>'Y', #Tyrosine
                'TAA'=>'_', #Stop
                'TAG'=>'_', #Stop
                'TGC'=>'C', #Cysteine
                'TGT'=>'C', #Cysteine
                'TGA'=>'_', #Stop
                'TGG'=>'W', #Tryptophan
                'CTA'=>'L', #Leucine
                'CTC'=>'L', #Leucine
                'CTG'=>'L', #Leucine
                'CTT'=>'L', #Leucine
                'CCA'=>'P', #Proline
                'CAT'=>'H', #Histidine
                'CAA'=>'Q', #Glutamine
                'CAG'=>'Q', #Glutamine
                'CGA'=>'R', #Arginine
                'CGC'=>'R', #Arginine
                'CGG'=>'R', #Arginine
                'CGT'=>'R', #Arginine
                'ATA'=>'I', #Isoleucine
                'ATC'=>'I', #Isoleucine
                'ATT'=>'I', #Isoleucine
                'ATG'=>'M', #Methionine
                'ACA'=>'T', #Threonine
                'ACC'=>'T', #Threonine
                'ACG'=>'T', #Threonine
                'ACT'=>'T', #Threonine
                'AAC'=>'N', #Asparagine
                'AAT'=>'N', #Asparagine
                'AAA'=>'K', #Lysine
                'AAG'=>'K', #Lysine
                'AGC'=>'S', #Serine#Valine
                'AGT'=>'S', #Serine
                'AGA'=>'R', #Arginine
                'AGG'=>'R', #Arginine
                'CCC'=>'P', #Proline
                'CCG'=>'P', #Proline
                'CCT'=>'P', #Proline
                'CAC'=>'H', #Histidine
                'GTA'=>'V', #Valine
                'GTC'=>'V', #Valine
                'GTG'=>'V', #Valine
                'GTT'=>'V', #Valine
                'GCA'=>'A', #Alanine
                'GCC'=>'A', #Alanine
                'GCG'=>'A', #Alanine
                'GCT'=>'A', #Alanine
                'GAC'=>'D', #Aspartic Acid
                'GAT'=>'D', #Aspartic Acid
                'GAA'=>'E', #Glutamic Acid
                'GAG'=>'E', #Glutamic Acid
                'GGA'=>'G', #Glycine
                'GGC'=>'G', #Glycine
                'GGG'=>'G', #Glycine
                'GGT'=>'G', #Glycine
                );
}
