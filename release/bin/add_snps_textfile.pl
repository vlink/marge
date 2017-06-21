#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use config;
use processing;
my $config = config::read_config();
use Data::Dumper;

$_ = "" for my($genome, $file, $output, $start, $filename, $gene_file, $refseq_file, $exon_ref, $header, $data, $line);
$_ = 0 for my($chr, $exons, $hetero, $allele);
$_ = () for my(@strains, @split,  %peaks, %exons, @exon_split, %all_exons, @tmp_split);

sub printCMD{
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
        print STDERR "\t-file <input file>: coordinates have to be the reference coordinates\n";
	print STDERR "\t-output <output name>: default <file name>_snps.txt\n";
        print STDERR "\t-ind: one or several individuals (comma separated)\n";
	print STDERR "\n\nOptions for RNA-Seq:\n";
	print STDERR "\tAnnotation RNA-Seq esp. with exons takes a while\n";
        print STDERR "\t-genome <genome> (e.g. mm10/hg38 - to lookup exon and refseq annotations)\n";
	print STDERR "\t-exons: Just looks for mutations within the exons\n";
	print STDERR "\t-refseq_file: File with RefSeq IDs (default in HOMER path defined in config)\n";
	print STDERR "\t-gene_file: File with Gene IDs (default in HOMER path defined in config)\n";
	print STDERR "\t-hetero: Strains are heterozygous\n";
	print STDERR "\t-data_dir: path to mutation data (default specified in config)\n";
        exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-file' => 1, '-ind' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
		"gene_file=s" => \$gene_file,
		"refseq_file=s" => \$refseq_file,
		"output=s" => \$output,
                "ind=s{,}" => \@strains,
		"exons" => \$exons,
		"hetero" => \$hetero,
		"data_dir=s" => \$data)
        or die (&printCMD());

#Set variables
if($data eq "") {
	$data = config::read_config()->{'data_folder'};
} 
if($hetero == 1) {
	$allele = 2;
} else {
	$allele = 1;
}
for(my $i = 0; $i < @strains; $i++) {
        chomp $strains[$i];
        $strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}
if($output eq "") {
	#Check if input filename has a file ending
	@split = split('\.', $file);
	$output = $split[0] . "_snps.txt";
}

#Save each peak
open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if(substr($line, 0, 1) eq "#" || substr($split[1], 0, 3) ne "chr" || length($split[1]) < 4 ) { $header = $line; next; }
	@tmp_split = split("_", $split[1]);
	$chr = substr($tmp_split[0], 3);
	$peaks{$chr}{$split[2]}{'end'} = $split[3];
	$peaks{$chr}{$split[2]}{'id'} = $split[0];
	$peaks{$chr}{$split[2]}{'line'} = $line;
}
close FH;

#Only the exons of the genes from this RNA-Seq experiment should be annotated
if($exons == 1) {
	#Check that genome is specified
	if($genome eq "") {
		print STDERR "Please specify the reference genome to find refseq annotations\n";
		print STDERR "e.g. mm10/hg38\n";
		exit;
	}
	($refseq_file, $gene_file) = processing::check_homer_files($refseq_file, $gene_file, $config, 0, $genome);
	$exon_ref = processing::save_transcript($refseq_file);
	%all_exons = %{$exon_ref};
}
#Start to collect mutations
print STDERR "Processing:\n";
#Iterate over all chrosomes
foreach my $chr (sort {$a cmp $b } keys %peaks) {
	print STDERR "\t\tchromosome " . $chr . " (" . (keys %{$peaks{$chr}}) . " peaks/transcripts) \n";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			#Open file that saves mutation for strain
			$filename = $data ."/" . uc($strains[$i]) . "/chr" . $chr . "_allele_" . $al . ".mut";
			if(!-e $filename) {
				print STDERR "\t\tCan't open $filename\n";
				foreach my $start (keys %{$peaks{$chr}}) {
					for(my $al = 1; $al <= $allele; $al++) {
						$peaks{$chr}{$start}{$strains[$i]}{$al} = "";
					}
				}
				next;
			}
	                open my $fh_mut, "<", "$filename" or die "Can't open $filename\n";
			$line = read_file_line($fh_mut);
			@split = split('\t', $line);
			#Iterate over every mutation - sorted so we don't have to jump around in the mutation file
			foreach my $start (sort {$a <=> $b} keys %{$peaks{$chr}}) {
				$peaks{$chr}{$start}{$strains[$i]}{$al} = "";
				%exons = ();
				#Add exons to hash
				if($exons == 1) {
					%exons = %{$all_exons{$peaks{$chr}{$start}{'id'}}};
				} else {
					$exons{$start . "_" . $peaks{$chr}{$start}{'end'}} = 0;
				}
				#Iterate over all exons - sorted so we don't have to jump around in the mutation file
				foreach my $e (sort {$exons{$a} <=> $exons{$b} } keys %exons) {
					@exon_split = split("_", $e);
					#Run through mutation file till current mutation is greater than start position of exon
					while($line and $split[0] < $exon_split[0]) {
						$line = read_file_line($fh_mut);
						if($line) {
							@split = split('\t', $line);
						}
					}
					#Run through mutation file till current mutation is greater than end position of exon and save all of these mutations
					while($line and $split[0] < $exon_split[1]) {
						$peaks{$chr}{$start}{$strains[$i]}{$al} .= $split[0] . ":" . $split[1] . "->" . $split[2] . ",";
						$line = read_file_line($fh_mut);
						if($line) {
							@split = split('\t', $line); 
						}
					}
				}
			}
			close $fh_mut;
		}
	}
}

#Write output
open OUT, ">$output";
print OUT $header;
for(my $i = 0; $i < @strains; $i++) {
	for(my $al = 1; $al <= $allele; $al++) {
		print OUT "\t" . $strains[$i] . " - " . $al;
	}
}
print OUT "\n";
foreach my $chr (keys %peaks) {
	foreach my $start (keys %{$peaks{$chr}}) {
		print OUT $peaks{$chr}{$start}{'line'};
		for(my $i = 0; $i < @strains; $i++) {
			for(my $al = 1; $al <= $allele; $al++) {
				if(length($peaks{$chr}{$start}{$strains[$i]}{$al}) > 1) {
					chop $peaks{$chr}{$start}{$strains[$i]}{$al};
				}
				print OUT "\t" . $peaks{$chr}{$start}{$strains[$i]}{$al};
			}
		}
		print OUT "\n";
	}
}
close OUT;


#Read in file line by line
sub read_file_line {
        my $fh = shift;
        if ($fh and my $line = <$fh>) {
                chomp $line;
                return $line;
        }
        return;
}
