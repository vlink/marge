#!/usr/bin/env perl
use strict;
use warnings;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use Getopt::Long;
use config;
use general;
my $config = config::read_config();
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/analysis'};
use processing;
use analysis_tree;
use Set::IntervalTree;
use Data::Dumper;

sub printCMD {
	print STDERR "\n\nUsage\n\n";
        print STDERR "\t-method <gene|protein|both|genomic|genomic_protein|genomic_both|align_gene|align_protein|align_both|align_genomic_gene|align_genomic_protein|align_genomic_both>\n";
	print STDERR "\t\tgene: outputs the gene nucleotid sequence based on RefSeq of Gene ID annotation\n";
	print STDERR "\t\tprotein: outputs the amino acid sequence of a gene based on RefSeq or Gene ID annotation\n";
	print STDERR "\t\tboth: outputs the nucleotid and amino acid sequence of a gene based on RefSeq or GENE ID annotation\n";
	print STDERR "\t\t\tFile with RefSeq or Gene ID coordinates have to be specified with -refseq_file or _gene_file\n";
	print STDERR "\t\t\tIDs for genes have to be specified with -gene or -transcript\n";
	print STDERR "\t\tgenomic: outputs the nucleotid sequence for any specified locus\n";
	print STDERR "\t\t\tgenomic locus has to be specified with -chr, -start, -end and if desired -strand\n";
	print STDERR "\t\tgenomic_protein: outputs the amino acid sequence of the specified locus\n";
	print STDERR "\t\tgenomic_both: outputs the nucleotid and amino acid sequence of specified locus\n";
	print STDERR "\t\talign_gene: outputs the alignment of nucleotid gene sequence based on mutation annotations in VCF file\n";
	print STDERR "\t\talign_protein: outputs the alignment of protein sequence. For this the nucleotid sequences are aligned based on mutations from the VCF file and then translated into amino acid sequences\n";
	print STDERR "\t\talign_both: outputs the alignment of the nucleotid and amino acid sequence for a gene based on the mutations in VCF file\n";
	print STDERR "\t\talig_genomic_gene: outputs the alignment the nucleotid sequence of a locus specified with -chr, -start, -end (and -strand)\n";
	print STDERR "\t\talign_genomic_protein: outputs the alignment of the amino acid sequence of a locus sepcified with -chr, -start, -end (and -strand). For this the nucleotid sequences are aligned based on mutations from the VCF file and then translated into amino acid sequences\n";
	print STDERR "\t\talign_genomic_both: outputs the alignment of the nucleotid and amino acid sequence for a locus specified with -chr, -start, -end (and -strand)\n\n";
        print STDERR "\t-ind <individuals>: Comma-separated list of individuals\n";
	print STDERR "\nParameters for genes\n";
        print STDERR "\t-genome: Genome\n";
	print STDERR "\t-refseq_file: File with RefSeq IDs (check HOMER format)\n";
	print STDERR "\t-gene_file: File with Gene IDs (check HOMER format)\n";
        print STDERR "\t-gene <gene name> (Comma separated list)\n";
        print STDERR "\t-transcript <RefSeq transcript name> (Comma separated list)\n";
	print STDERR "\nParameters for genomic locations\n";
        print STDERR "\t-start <pos>: Start position of genomic region (Comma separated list - same length as end, chr, and strand)\n";
        print STDERR "\t-end <pos>: Stop position of genomic region (Comma separated list - same length as start, chr, and strand)\n";
        print STDERR "\t-chr <chromosome>:Chromosome for genomic region (Comma separated list - same length as start, end, and strand)\n";
        print STDERR "\t-strand <+|->: Default: + (Comma separated list - same length as start, end, and chr)\n";
	print STDERR "\nAdditional parameters:\n";
        print STDERR "\t-hetero: Data is heterozygous\n";
	print STDERR "\t-data_dir: Path to data directory for individuals - default specified in config\n";
	print STDERR "\t-genome_dir: Path to fastq files for individuals - default specified in config\n";
	print STDERR "\n\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

$_ = "" for my ($genome, $path, $refseq_file, $gene_file, $data, $genome_dir, $refseq_NM, $method, $strand, $chr, $gene, $transcript);
$_ = 0 for my ($hetero, $html, $exon, $gene_found, $allele, $shift, $longest_seq, $filter_no_mut, $start);
$_ = () for my (%strains, @split, @strains, @strains2, $seqs, @parts, @introns, %exons, %peaks, %gene, %codons, %tree, %last_strain, %lookup_strain, $f1, $mut_line, $f1_file, @fileHandles, @exons, %align_nt, %align_nt_seq, @gene, @transcript, @list, @list_id, @start, @end, @chr, @strand, $peak_ref, $exon_ref, %strand);


GetOptions(     "genome=s" => \$genome,
		"refseq_file=s" => \$refseq_file,
		"gene_file=s" => \$gene_file,
                "method=s" => \$method,
                "gene=s{,}" => \@gene,
                "transcript=s{,}" => \@transcript,
                "ind=s{,}" => \@strains,
                "start=s{,}" => \@start,
                "end=s{,}" => \@end,
                "chr=s{,}"=> \@chr,
                "hetero" => \$hetero,
		"data_dir=s" => \$data,
		"genome_dir=s" => \$genome_dir,
		"strand=s{,}" => \@strand, 
                "html" => \$html)
        or die("Error in command line options!\n");

if($method ne "gene" && $method ne "protein" && $method ne "both" && $method ne "genomic" && $method ne "align_gene" && $method ne "align_protein" && $method ne "align_both" && $method ne "align_genomic_gene" && $method ne "align_genomic_protein" && $method ne "align_genomic_both" && $method ne "genomic_protein" && $method ne "genomic_both") {
	&printCMD();
}

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
	($refseq_file, $gene_file) = processing::check_homer_files($refseq_file, $gene_file, $config, scalar @gene, $genome);
}
&write_codons();

for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

if($data eq "") {
	$data = $config->{'data_folder'};
}
if($genome_dir eq "") {
	$genome_dir = $data;
}

print STDERR "Loading shift vectors\n";
for(my $i = 0; $i < @strains; $i++) {
        my($tree_ref, $last, $lookup) = general::read_strains_data($strains[$i], $data, "ref_to_strain");
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
		$transcript = processing::get_refseq_for_gene($gene_file, $gene);
		if($transcript eq "") {
			print STDERR "Could not find Ref-Seq entry for " . $gene . "\n";
			print STDERR "Skip\n";
			next;
		}
		($strand, $chr, $peak_ref, $exon_ref) = processing::save_transcript($refseq_file, $transcript);
		%peaks = %{$peak_ref};
		%exons = %{$exon_ref};
		%strand = %{$strand};
	} elsif($list_id[$entry_id] eq "refseq") {
		$transcript = $list[$entry_id];
		($strand, $chr, $peak_ref, $exon_ref) = processing::save_transcript($refseq_file, $transcript);
		%peaks = %{$peak_ref};
		%exons = %{$exon_ref};
		%strand = %{$strand};
	} else {
		if(substr($chr[$entry_id - $shift], 0, 3) ne "chr") {
			$chr[$entry_id - $shift] = "chr" . $chr[$entry_id - $shift];
		}
		$strand = $strand[$entry_id - $shift];
		$chr = $chr[$entry_id - $shift];
		#print $chr . "\t" . $start[$entry_id - $shift] . "\t" . $end[$entry_id - $shift] . "\t" . $strand . "\n";
		$peaks{substr($chr[$entry_id - $shift], 3)}{$start[$entry_id - $shift]} = $end[$entry_id - $shift];
		$exons{$start[$entry_id - $shift] . "_" . $end[$entry_id - $shift]} = 1;
	}
	$seqs = ();
	if(-e "tmp") { `rm tmp`; }
	($seqs, $longest_seq, $filter_no_mut) = analysis::get_seq_for_peaks("tmp", \%peaks, \@strains, $genome_dir, $allele, $exon, 0, 0, \%tree, \%lookup_strain, \%last_strain, \%strand);
	if($method eq "gene" || $method eq "genomic") {
		&assemble_gene(1);
	} elsif($method eq "protein" || $method eq "genomic_protein") {
		&assemble_gene(0);
		foreach my $s (keys %gene) {
			foreach my $al (keys %{$gene{$s}}) {
				print STDERR $s . " - allele " . $al . "\n";
				if($strand eq "-") {
					$gene{$s}{$al} = analysis::rev_comp($gene{$s}{$al});
				}
				&translate_seq($gene{$s}{$al}, 1);
			}
		}
	} elsif($method eq "both" || $method eq "genomic_both") {
		&assemble_gene(1);
		foreach my $s (keys %gene) {
			foreach my $al (keys %{$gene{$s}}) {
				print STDERR $s . " - allele " . $al . "\n";
				if($strand eq "-") {
					$gene{$s}{$al} = analysis::rev_comp($gene{$s}{$al});
				}
				&translate_seq($gene{$s}{$al}, 1);
			}
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
	} else {
		&printCMD();
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
		for(my $al = 1; $al <= $allele; $al++) {
			$local_seq = $align_nt_seq{$i}{$al};
			$local_seq =~ s/-//g;
			$prot_seq{$i}{$al} = &translate_seq($local_seq, 0);
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			$prot = "";
			$num_gaps = 0;
			$num_gaps_tmp = 0;
			$align_prot = "";
			$prot_pos = 0;
			$local_seq = $align_nt_seq{$i}{$al};
			$local_seq_align = $align_nt{$i}{$al};	
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
						$prot .= substr($prot_seq{$i}{$al}, $prot_pos, 1);
						$align_prot .= "|";
						$prot_pos++;

					}
				} elsif(index(substr($local_seq_align, 0, 3), '.') != -1) {
					$prot .= substr($prot_seq{$i}{$al}, $prot_pos, 1);
					$align_prot .= ".";
					$prot_pos++;
				} else {
		#			print "everything is fine\n";
					$prot .= substr($prot_seq{$i}{$al}, $prot_pos, 1);
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
			$align_prot_seq{$i}{$al} = $prot;
			$align_seq{$i}{$al} = $align_prot;
		}
	}
	#Check for synoymous mutations
	my @alignment_seq;
	my $diff_as = 0;
	my $current_pos = 0;
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			@alignment_seq = split('', $align_seq{$i}{$al});
			$current_pos = 0;
			foreach my $c (@alignment_seq) {
				#Check aa of all strains
				if($c eq ".") {
					$diff_as = 0;
					for(my $k = 0; $k < @strains; $k++) {
						for(my $al_tmp = 1; $al_tmp <= $allele; $al_tmp++) {
							if(substr($align_prot_seq{$i}{$al}, $current_pos, 1) ne substr($align_prot_seq{$k}{$al_tmp}, $current_pos, 1)) {
								$diff_as = 1;
							}
						}	
					}
					if($diff_as == 0) {
						$align_seq{$i}{$al} = substr($align_seq{$i}{$al}, 0, $current_pos) . "*" . substr($align_seq{$i}{$al}, $current_pos + 1);
					}
				}
				$current_pos++;
			}
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			print STDERR $strains[$i] . " - allele - " . $al . "\n";
			print STDERR $align_prot_seq{$i}{$al} . "\n";
			print STDERR $align_seq{$i}{$al} . "\n";
		}
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
		print STDERR $protein_seq . "\n";
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
		for(my $al = 1; $al <= $allele; $al++) {
			@{$orig_seq[$i][$al]} = split('', $gene{uc($strains[$i])}{$al});
			$mut_shift[$i][$al] = 0;
		}
	}
	my $current_mut;
	my $strain_mut;
	my $mut_exists = 0;
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			$run_index = 0;
			$complete_index = 0;
			if(!defined $fileHandles[$i][$al]) { next; }
			$mut_line = read_file_line($fileHandles[$i][$al]);
			@split = split('\t', $mut_line);
			foreach my $e (sort {$a cmp $b} keys %exons) {
				$mut_exists = 0;
				@exons = split("_", $e);
				while($split[0] < $exons[0]) {
					$mut_line = read_file_line($fileHandles[$i][$al]);
					@split = split('\t', $mut_line);
				}
				while($split[0] < $exons[1]) {
					$run_index = $complete_index + ($split[0] - $exons[0]) - 1;
					$mutation_hash{$i}{$al}{$run_index} = length($split[1]) - length($split[2]);
					$mut_line = read_file_line($fileHandles[$i][$al]);
					@split = split('\t', $mut_line);
				}
				$complete_index = $complete_index + $exons[1] - $exons[0];
			}
		}
	}
	#Check if there are mutations in every strain that are always the same
	my $num_mut = 0;
	my $diff = 0;
	foreach my $gene (keys %gene) {
		for(my $al = 1; $al <= $allele; $al++) {
			if(length($gene{$gene}{$al}) > $last_index) {
				$last_index = length($gene{$gene}{$al});
			}
		}
	}
	foreach my $strain (keys %mutation_hash) {
		foreach my $al (keys %{$mutation_hash{$strain}}) {
			foreach my $pos (keys %{$mutation_hash{$strain}{$al}}) {
				$num_mut = 0;
				$diff = 0;
				for(my $i = 0; $i < @strains; $i++) {
					for(my $a1 = 1; $a1 <= $allele; $a1++) {
						if(exists $mutation_hash{$i}{$al}{$pos}) {
							$num_mut++;
							if(exists $mutation_hash{$i}{$a1}{$pos} && $mutation_hash{$i}{$a1}{$pos} ne $mutation_hash{$strain}{$al}{$pos}) {
								$diff = 1;
							}
						}
					}
				}
				if($num_mut == @strains && $diff == 0) {
					for(my $i = 0; $i < @strains; $i++) {
						delete $mutation_hash{$i}{$al}{$pos};
					}
				}
			}
		}
	}
	for(my $run = 0; $run < $last_index; $run++) {
		for(my $i = 0; $i < @strains; $i++) {
			for(my $al = 1; $al <= $allele; $al++) {
				if(!defined $orig_seq[$i][$al][$run - $mut_shift[$i][$al]]) { next; }
				if(exists $mutation_hash{$i}{$al}{$run}) {
					if($mutation_hash{$i}{$al}{$run} == 0) {
						$align[$i][$al][$run] = ".";
						$align_seq[$i][$al][$run] = $orig_seq[$i][$al][$run - $mut_shift[$i][$al]];
					} elsif($mutation_hash{$i}{$al}{$run} > 0) {
						if($strand eq "-") {
							$align[$i][$al][$run] =  "|" . ( '-' x $mutation_hash{$i}{$al}{$run} ); 
							$align_seq[$i][$al][$run] =  $orig_seq[$i][$al][$run - $mut_shift[$i][$al]] . ( '-' x $mutation_hash{$i}{$al}{$run} );

						} else {
							$align[$i][$al][$run] = "|" . ( '-' x $mutation_hash{$i}{$al}{$run} ); 
							$align_seq[$i][$al][$run] = $orig_seq[$i][$al][$run - $mut_shift[$i][$al]] . ( '-' x $mutation_hash{$i}{$al}{$run} );
						}
						$mut_shift[$i][$al] += $mutation_hash{$i}{$al}{$run};
					} else {
						#First check if there are several insertions and get the longest
						my $current = $mutation_hash{$i}{$al}{$run};
						for(my $k = 0; $k < @strains; $k++) {
							for(my $al_tmp = 1; $al_tmp <= $allele; $al_tmp++) {
								if(exists $mutation_hash{$k}{$al_tmp}{$run} && $mutation_hash{$k}{$al_tmp}{$run} < $mutation_hash{$i}{$al_tmp}{$run}) {
									$i = $k;
									$al = $al_tmp;
								}
							}
						}
						for(my $k = 0; $k < @strains; $k++) {
							for(my $al_tmp = 1; $al_tmp <= $allele; $al_tmp++) {
								if($strand eq "-") {
									if(!exists $mutation_hash{$k}{$al_tmp}{$run}) {
										$align[$k][$al_tmp][$run] = "|" . ( '-' x abs($mutation_hash{$i}{$al}{$run}) );
										$align_seq[$k][$al_tmp][$run] = $orig_seq[$k][$al_tmp][$run - $mut_shift[$k][$al_tmp]] . ( '-' x abs($mutation_hash{$i}{$al}{$run} ));
									} elsif(exists $mutation_hash{$k}{$al_tmp}{$run} && $mutation_hash{$k}{$al_tmp}{$run} != $mutation_hash{$i}{$al}{$run}) {
										if($mutation_hash{$k}{$al_tmp}{$run} < 0) {
											$align[$k][$al_tmp][$run] = ( '-' x (abs($mutation_hash{$i}{$al}{$run}) - abs($mutation_hash{$k}{$al_tmp}{$run})));	
											$align_seq[$k][$al_tmp][$run] = ( '-' x (abs($mutation_hash{$i}{$al}{$run}) - abs($mutation_hash{$k}{$al_tmp}{$run})));
										} else {
											$align[$k][$al_tmp][$run] =  ( '-' x (abs($mutation_hash{$i}{$al}{$run}) + $mutation_hash{$k}{$al_tmp}{$run}) ) . "|"; 
											$align_seq[$k][$al_tmp][$run] =  ( '-' x (abs($mutation_hash{$i}{$al}{$run}) + $mutation_hash{$k}{$al_tmp}{$run}) ) . $orig_seq[$k][$al_tmp][$run - $mut_shift[$k][$al_tmp]];
											$mut_shift[$k][$al_tmp] += $mutation_hash{$k}{$al_tmp}{$run};
										}	
									} else {
										$align[$k][$al_tmp][$run] = "|";
										$align_seq[$k][$al_tmp][$run] = $orig_seq[$i][$al][$run - $mut_shift[$k][$al_tmp]];
									}
								} else {
									if(!exists $mutation_hash{$k}{$al_tmp}{$run}) {
										$align[$k][$al_tmp][$run] = "|" . ( '-' x abs($mutation_hash{$i}{$al}{$run}) );
										$align_seq[$k][$al_tmp][$run] = $orig_seq[$k][$al_tmp][$run - $mut_shift[$k][$al_tmp]] . ( '-' x abs($mutation_hash{$i}{$al}{$run} ));
									} elsif(exists $mutation_hash{$k}{$al_tmp}{$run} && $mutation_hash{$k}{$al_tmp}{$run} != $mutation_hash{$i}{$al}{$run}) {
										if($mutation_hash{$k}{$al_tmp}{$run} < 0) {
											$align[$k][$al_tmp][$run] = ( '-' x (abs($mutation_hash{$i}{$al}{$run}) - abs($mutation_hash{$k}{$al_tmp}{$run})));	
											$align_seq[$k][$al_tmp][$run] = ( '-' x (abs($mutation_hash{$i}{$al}{$run}) - abs($mutation_hash{$k}{$al_tmp}{$run})));
										} else {
											$align[$k][$al_tmp][$run] =  "|" . ( '-' x (abs($mutation_hash{$i}{$al}{$run}) + $mutation_hash{$k}{$al_tmp}{$run}) ); 
											$align_seq[$k][$al_tmp][$run] = $orig_seq[$k][$al_tmp][$run - $mut_shift[$k][$al_tmp]] . ( '-' x (abs($mutation_hash{$i}{$al}{$run}) + $mutation_hash{$k}{$al_tmp}{$run}) );
											$mut_shift[$k][$al_tmp] += $mutation_hash{$k}{$al_tmp}{$run};
										}	
									} else {
										$align[$k][$al_tmp][$run] = "|";
										$align_seq[$k][$al_tmp][$run] = $orig_seq[$i][$al][$run - $mut_shift[$k][$al_tmp]];
									}
								}
							}
						}
						#We are through all mutations
						for(my $k = 0; $k < @strains; $k++) {
							for(my $al_tmp = 1; $al_tmp <= $allele; $al_tmp++) {
								delete $mutation_hash{$k}{$al_tmp}{$run};
							}
						}
						$i = @strains + 1;
					}
				} else {
					if(!defined $orig_seq[$i][$al][$run]) {
						#	print "dead\n";
						next;
					}
					$align_seq[$i][$al][$run] = $orig_seq[$i][$al][$run];
					$align[$i][$al][$run] = "|";
				}
			}
		}
	}
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			$align_nt_seq{$i}{$al} = "" . join("", @{$align_seq[$i][$al]});
			$align_nt{$i}{$al} =  "" . join("", @{$align[$i][$al]});
			if($strand eq "-") {
				$align_nt_seq{$i}{$al} = analysis::rev_comp($align_nt_seq{$i}{$al});
				$align_nt{$i}{$al} = reverse($align_nt{$i}{$al});
			}
			if($print == 1) {
				print STDERR $strains[$i] . " - allele " . $al . "\n";
				print STDERR $align_nt_seq{$i}{$al} . "\n";
				print STDERR $align_nt{$i}{$al} . "\n";
			}
		}
	}
}

sub open_filehandles{
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			my $filename = $data . "/" . uc($strains[$i]) . "/chr" . substr($chr, 3) . "_allele_" . $al . ".mut";
			if(-e $filename) {
				open my $fh, "<", $filename or die "Can't open $filename: $!\n";
				$fileHandles[$i][$al] = $fh;
			}
		}
	}
}

sub assemble_gene{
	my $print = $_[0];
	my $assembly = "";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			$assembly = "";
			foreach my $e (sort {$a cmp $b} keys %exons) {
				$assembly .= $seqs->{$chr . "_" . $e . "_" . $strains[$i] . "_" . $al};
			}
			$gene{$strains[$i]}{$al} = $assembly;
			if($print == 1) {
				if($strand eq "-") {
					$assembly = analysis::rev_comp($assembly);
				}
				print STDERR $strains[$i] . " - allele " . $al . "\n";
				print STDERR $assembly;
				print STDERR "\n";
			}
		}
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
