#!/usr/bin/perl


use strict;
use Getopt::Long;
#require '/Users/verenalink/workspace/strains/general/config.pm';
require '../general/config.pm';
require 'analysis.pm';
use Data::Dumper;
#require '../db_part/database_interaction.pm';


$_ = "" for my($genome, $file, $tf, $filename, $last_line);
$_ = () for my(@strains, %promoter, %exp, @split, @name, %mutation_pos, %shift, %current_pos, %save_local_shift, %seq, %PWM, @fileHandlesMotif, %index_motifs, %tag_counts, %fc, @header, %block, %analysis_result, %existance, %diff, %ranked_order, %mut_one, %mut_two, %delta_score, %exp_per_strain, %genename_for_promoter);
$_ = 0 for my($homo, $allele, $region, $motif_score, $motif_start, $more_motifs, $save_pos, $delta, $line_number);
my $data = config::read_config()->{'data_folder'};
#$data = "/Users/verenalink/workspace/strains/data/";

my %comp;
$comp{'A'} = 'T';
$comp{'C'} = 'G';
$comp{'G'} = 'C';
$comp{'T'} = 'A';

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - Order must overlay with order in annotated peak file\n";
        print STDERR "\t-file <file>: RNA-Seq expression file\n";
	print STDERR "\t-TF <transcription factor motif matrix>: Matrix of the TF that was chipped for\n";
        print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-region: Size of promoter (Default: 200)\n";
	print STDERR "\t-delta: Uses motif score differences instead of binary existance\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}

my $param = config::read_config();

$region = 200;

GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
                "strains=s{,}" => \@strains,
                "homo" => \$homo, 
		"-TF=s" => \$tf, 
		"-region=s" => \$region,
		"-delta" => \$delta)
        or die("Error in command line options!\n");


#First step: Get the sequences for the peaks
$allele = 1;
my $ref_save = 0;
$_ = () for my(%current_pos, %mutation_pos, %shift, @split, $filename); 

#Save motif files
my ($index_motif_ref, $PWM_ref, $fileHandlesMotif_ref) = analysis::read_motifs($tf);
%index_motifs = %$index_motif_ref;
%PWM = %$PWM_ref;
@fileHandlesMotif = @$fileHandlesMotif_ref;

#Make sure the reference is in the strains array
for(my $i = 0; $i < @strains; $i++) {
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

print STDERR "Save expression profile\n";

open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 5) eq "Trans") { next; }
	@split = split('\t', $line);
	my $pos;
	if($split[4] eq "+") {
		$pos = $split[2] - $region; 
		$promoter{substr($split[1], 3)}{$split[2] - $region} = $split[2];
	} else {
		$pos = $split[3];
		$promoter{substr($split[1], 3)}{$split[3]} = $split[3] + $region;
	}
	$exp{substr($split[1], 3)}{$pos} = "";
	for(my $i = 0; $i < @strains; $i++) {
		$exp{substr($split[1], 3)}{$pos} .= $split[8 + $i] . "\t";
		$exp_per_strain{$split[0]}{$strains[$i]} = $split[8 + $i];
	}
	$genename_for_promoter{substr($split[1], 3)}{$pos} = $split[0];
	$line_number++;
}

my %peaks;
my $tmp_out = "rna_tmp" . rand(15);
print STDERR "Extracting sequences from strain genomes\n";
print STDERR $tmp_out . "\n";

my ($seq_ref, $save_local_shift_ref) = analysis::get_seq_for_peaks($tmp_out, \%promoter, \@strains, $data, $allele, $line_number);
%seq = %$seq_ref;
%save_local_shift = %$save_local_shift_ref;

my $tmp_out_main_motif = "rna_tmp" . rand(15);
analysis::scan_motif_with_homer($tmp_out, $tmp_out_main_motif, $tf);

analysis::write_header(\@fileHandlesMotif, \@strains, 1);

analysis::analyze_motifs($tmp_out_main_motif, \%seq, \%save_local_shift, \@fileHandlesMotif, \%exp, \@strains, \%index_motifs);

for(my $i = 0; $i < @fileHandlesMotif; $i++) {
	close $fileHandlesMotif[$i];
}

print STDERR "Start analysis!\n";
