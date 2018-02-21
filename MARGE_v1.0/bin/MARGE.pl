#!/usr/bin/env perl
use strict;
use warnings;

# Copyright Verena M. Link <vlink@ucsd.edu>
# 
# This file is part of MARGE
#
# MARGE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MARGE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


#This is the wrapper script for all scripts that can be used within MARGE
#Used for DOCKER
#Easy access for user

my $command = "";

if(@ARGV < 1) {
	&printCMD();
} elsif($ARGV[0] eq "-h" || $ARGV[0] eq "--help") {
	&printCMD();
} elsif($ARGV[0] eq "prepare_files") {
	$command = "perl CHANGE_TO_WORKDIR/prepare_textfiles.pl";	
	&generate_command();
} elsif($ARGV[0] eq "create_genomes") {
	$command = "perl CHANGE_TO_WORKDIR/create_genome_textfile.pl";
	&generate_command();
} elsif($ARGV[0] eq "update_files") {
	$command = "perl CHANGE_TO_WORKDIR/update_last_shift_and_lookup_table.pl";
	&generate_command();
} elsif($ARGV[0] eq "genome_interactions") {
	$command = "perl CHANGE_TO_WORKDIR/genome_interactions.pl";	
	&generate_command();
} elsif($ARGV[0] eq "mutation_bedfiles") {
	$command = "perl CHANGE_TO_WORKDIR/generate_mutation_bedfile.pl";
	&generate_command();
} elsif($ARGV[0] eq "mutation_info") {
	$command = "perl CHANGE_TO_WORKDIR/count_mutations_per_strain.pl";
	&generate_command();
} elsif($ARGV[0] eq "extract_sequences") {
	$command = "perl CHANGE_TO_WORKDIR/extract_seq_from_peakfiles.pl";
	&generate_command();
} elsif($ARGV[0] eq "annotate_mutations") {
	$command = "perl CHANGE_TO_WORKDIR/add_snps_textfile.pl";
	&generate_command();
} elsif($ARGV[0] eq "shift") {
	$command = "perl CHANGE_TO_WORKDIR/shift_files_to_reference_textfiles.pl";
	&generate_command();
} elsif($ARGV[0] eq "allele_specific_reads") {
	$command = "perl CHANGE_TO_WORKDIR/get_mapped_seq_that_span_mut.pl";
	&generate_command();
} elsif($ARGV[0] eq "annotate_allele_specific") {
	$command = "perl CHANGE_TO_WORKDIR/annotate_heterozygous_data.pl";
	&generate_command();
} elsif($ARGV[0] eq "denovo_motifs") {
	$command = "perl CHANGE_TO_WORKDIR/findMotifsGenome_strains.pl";
	&generate_command();
} elsif($ARGV[0] eq "mutation_analysis") {
	$command = "perl CHANGE_TO_WORKDIR/analyze_ChIP_mutations_tree.pl";
	&generate_command();
} elsif($ARGV[0] eq "summary_heatmap") {
	$command = "perl CHANGE_TO_WORKDIR/generate_heatmap.pl";
	&generate_command();	
} elsif($ARGV[0] eq "download") {
	$command = "perl CHANGE_TO_WORKDIR/update.pl";
	&generate_command();
} else {
	&printCMD();
}
`$command`;

sub generate_command{
	if(@ARGV > 1 && $ARGV[1] ne "--help" && $ARGV[1] ne "-h") {
		for(my $i = 1; $i < @ARGV; $i++) {
			$command  .= " " . $ARGV[$i];
		}
	}

}

sub printCMD {
	print STDERR "MARGE USAGE\n";
	print STDERR "Download pre-processed data\n";
	print STDERR "\tMARGE.pl download\n\n";
	print STDERR "Getting started with data processing\n";
	print STDERR "\tMARGE.pl prepare_files\n";
	print STDERR "\t\tThis script takes VCF input files and generates MARGE mutation files, shifting files and the genomes for the individuals from the VCF file.\n\n";
	print STDERR "\tMARGE.pl create_genomes\n";
	print STDERR "\t\tThis script generates the individual genomes for all individuals specified when called. It only works when the mutation and shifting files were already prepared using MARGE.pl prepare_files.\n\n";
	print STDERR "\tMARGE.pl update_files\n";
	print STDERR "\t\tFor some projects (e.g. the 1000 genome project) mutations are annotated per chromosomes. Therefore, MARGE.pl prepare_files was run per chromosome. In this case some of the MARGE files were overwritten and need to be adjusted after finishing MARGE.pl prepare_files for all chromosomes.\n\n";
	print STDERR "\n\nAccess to mutation data\n";
	print STDERR "\tMARGE.pl genome_interactions\n";
	print STDERR "\t\tThis script allows to extract sequences for different individuals. It also provides protein sequences and alignments. The alignments are based on the VCF file and are not generated using any alignment algorithms.\n\n";
	print STDERR "\tMARGE.pl mutation_bedfiles\n";
	print STDERR "\t\tThis script outputs a bed file per individual per allele that can be uploaded to the UCSC genome browser.\n\n";
	print STDERR "\tMARGE.pl mutation_info\n";
	print STDERR "\t\tThis script counts mutations between different individuals and outputs them in a table format. It can either count all mutations in comparison to the reference or all private mutations\n\n";
	print STDERR "\tMARGE.pl extract_sequences\n";
	print STDERR "\t\tThis script extracts the individual-specific sequences in fasta format.\n\n"; 
	print STDERR "\tMARGE.pl annotate_mutations\n";
	print STDERR "\t\tThis script adds all mutations within specified genomic  forloci for each individual at the end of the input file.\n\n";
	print "\n\nData handling\n";
	print STDERR "\tMARGE.pl shift\n";
	print STDERR "\t\tThis script shifts the data mapped to the individual genome to the reference coordinates. This step is necessary to compare genomic loci between different individuals, as well as for visualization of data in the different genome browsers.\n\n";
	print STDERR "\tMARGE.pl allele_specific_reads\n";
	print STDERR "\t\tFor heterozygous data, this script only keeps reads that were perfectly mapped (to avoid bias introduced by SNPs) and also outputs a file with reads that overlap mutations (to call allele-specific binding\n\n";
	print STDERR "\tMARGE.pl annotate_allele_specific\n";
	print STDERR "\t\tFor heterozygous data, this script annotates all genomic loci between different individuals. It takes reads that overlap mutations for loci with mutations in this individal. In case there is no mutation, it takes the perfectly mapped reads for this locus and divides it by 2 (because of two alleles)\n\n";
	print STDERR "\n\nMutation analysis:\n";
	print STDERR "\tMARGE.pl denovo_motifs\n";
	print STDERR "\t\tHOMER denovo motif analysis extended to integrate the usage of different genomes.\n\n";
	print STDERR "\tMARGE.pl mutation_analysis\n";
	print STDERR "\t\tPairwise or all versus all analysis of the impact of mutations on open chromtin (ATAC-Seq, DNAse Hypersensitivity assays etc.) or transcription factor binding.\n\n";
	print STDERR "\tMARGE.pl summary_heatmap\n";
	print STDERR "\t\tScript to summarize the results from several runs of pairwise mutation analyses or all versus all runs (needs at least two files as input).\n\n";
	print STDERR "For more information call the different scripts either without any parameters or the parameter --help or -h\n\n";
}

