#!/usr/bin/perl -w

use strict;
use Getopt::Long;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
use Data::Dumper;
my $config = config::read_config();
	
$|=1;

sub printCMD{
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
        print STDERR "\t-genome <genome>\n";
        print STDERR "\t-file <input file>: coordinates have to be the reference coordinates\n";
	print STDERR "\t-output <output name>: default <file name>_snps.txt\n";
        print STDERR "\t-strains: one or several strains (comma separated)\n";
	print STDERR "\n\nOptions for RNA-Seq:\n";
	print STDERR "\tAnnotation RNA-Seq esp. with exons takes a while\n";
	print STDERR "\t-exons: Just looks for mutations within the exons\n";
	print STDERR "\t-refseq_file: File with RefSeq IDs (check HOMER format)\n";
	print STDERR "\t-gene_file: File with Gene IDs (check HOMER format)\n";
	print STDERR "\t-hetero: Strains are heterozygous\n";
        exit;
}

if(@ARGV < 1) {
	&printCMD();
}


$_ = "" for my($genome, $file, $output, $split_sign, $chr_pos, $shift_order, $start, $end, $tmp_start, $transcript, $lin, $line, $filename, $gene_file, $refseq_file, $ref, $exon_ref, $strand, $header, $chr_tmp);
$_ = 0 for my($homer, $format, $chr, $print_header, $exons, $all, $found, $total, $count, $hetero);
$_ = () for my(@strains, @split, $dbh, $sth, $sth_exon,  %lookup_header, $mutations, %peaks, %exons, @exon_split, %all_exons);
my $allele = 2;

my %mandatory = ('-genome' => 1, '-file' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);
my $data = config::read_config()->{'data_folder'};


GetOptions(     "genome=s" => \$genome,
                "file=s" => \$file,
		"gene_file=s" => \$gene_file,
		"refseq_file=s" => \$refseq_file,
		"output=s" => \$output,
                "strains=s{,}" => \@strains,
		"exons" => \$exons,
		"all" => \$all,
		"hetero" => \$hetero,
                "format" => \$format)
        or die ("Error in command line arguments!\n");



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

open FH, "<$file";
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	if(substr($split[1], 0, 3) ne "chr" || length($split[1]) < 4 ) { $header = $line; next; }
	$peaks{substr($split[1], 3)}{$split[2]}{'end'} = $split[3];
	$peaks{substr($split[1], 3)}{$split[2]}{'id'} = $split[0];
	$peaks{substr($split[1], 3)}{$split[2]}{'line'} = $line;
}
close FH;
if($exons == 1) {
	($refseq_file, $gene_file) = processing::check_homer_files($refseq_file, $gene_file, $config, 0, $genome);
	$exon_ref = processing::save_transcript($refseq_file);
	%all_exons = %{$exon_ref};
}
#Start to collect mutations
print STDERR "Processing:\n";
foreach my $chr (sort {$a cmp $b } keys %peaks) {
	print STDERR "\t\tchromosome " . $chr . " (" . (keys %{$peaks{$chr}}) . " peaks/transcripts) \n";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $al = 1; $al <= $allele; $al++) {
			$filename = $data ."/" . uc($strains[$i]) . "/chr" . $chr . "_allele_" . $al . ".mut";
			if(!-e $filename) {
				print STDERR "\t\tCan't open $filename\n";
				foreach my $start (keys %{$peaks{$chr}}) {
					$peaks{$chr}{$start}{$strains[$i]} = "";
				}
				next;
			}
	                open my $fh_mut, "<", "$filename" or die "Can't open $filename\n";
			$line = read_file_line($fh_mut);
			@split = split('\t', $line);
			foreach my $start (sort {$a <=> $b} keys %{$peaks{$chr}}) {
				$peaks{$chr}{$start}{$strains[$i]} = "";
				%exons = ();
				if($exons == 1) {
					%exons = %{$all_exons{$peaks{$chr}{$start}{'id'}}};
				} else {
					$exons{$start . "_" . $peaks{$chr}{$start}{'end'}} = 0;
				}
				foreach my $e (sort {$exons{$a} <=> $exons{$b} } keys %exons) {
					@exon_split = split("_", $e);
					while($line and $split[0] < $exon_split[0]) {
						$line = read_file_line($fh_mut);
						if($line) {
							@split = split('\t', $line);
						}
					}
					while($line and $split[0] < $exon_split[1]) {
						$peaks{$chr}{$start}{$strains[$i]} .= $split[0] . "(" . $al . "):" . $split[1] . "->" . $split[2] . ",";
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

open OUT, ">$output";
print OUT $header;
for(my $i = 0; $i < @strains; $i++) {
	print OUT "\t" . $strains[$i];
}
print OUT "\n";
foreach my $chr (keys %peaks) {
	foreach my $start (keys %{$peaks{$chr}}) {
		print OUT $peaks{$chr}{$start}{'line'};
		for(my $i = 0; $i < @strains; $i++) {
			if(length($peaks{$chr}{$start}{$strains[$i]}) > 1) {
				chop $peaks{$chr}{$start}{$strains[$i]};
			}
			print OUT "\t" . $peaks{$chr}{$start}{$strains[$i]};
		}
		print OUT "\n";
	}
}
close OUT;
exit;

my @fileHandles;
my $outfile_open = 0;

sub open_filehandles{
	for(my $i = 0; $i < @strains - 1; $i++) {
		$filename = $data . "/" . uc($strains[$i]) . "/chr" . $_[0] . "_allele_" . $_[1] . ".mut";
		print STDERR $filename . "\n";
		print $filename . "\n";
		open my $fh_mut, "<", "$filename" or die "Can't open $filename: $!\n";
		push @fileHandles, $fh_mut;
	}
        $outfile_open = 1;
}

sub close_filehandles{
	for(my $i = 0; $i < @strains - 1; $i++) {
                close $fileHandles[$i];
	}
        undef @fileHandles;
}
my @current_line;
print STDERR "Add automated chromosome query\n";
my %save_peaks = ();
for(my $chr = 1; $chr < 19; $chr++) {
	open FH, "<$file";
	if($chr > 1) {
		&close_filehandles($chr, 1);
	}
	&open_filehandles($chr, 1);
	print STDERR "Parsing through file and extracting mutations from database for chromosome " . $chr . "\n";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq "#") {
			if($chr == 1) {
				print OUT $line . "\n";
			}
			next;
		}
		@split = split($split_sign, $line);
		if($split[1] eq "chr" . $chr) {
			my $start = int($split[2]) * 1;
			my $stop = int($split[3]) * 1;
			$save_peaks{$start}->{'stop'} = $stop;
			$save_peaks{$start}->{'line'} = $line;
		}
	}
	foreach my $start_pos (sort {$a <=> $b} keys %save_peaks) {
		for(my $run = 0; $run < @strains - 1; $run++) {
			$save_peaks{$start_pos}->{'line'} .= "\t";
			if(defined read_file_line($fileHandles[$run])) {
				@current_line = split('\t', read_file_line($fileHandles[$run])) if defined read_file_line($fileHandles[$run]);

			} else {
				next;
			} 
			while($current_line[1] < $save_peaks{$start_pos}->{'stop'}) {
				if($current_line[1] > $start_pos) {
					$save_peaks{$start_pos}->{'line'} .= $current_line[1] . ",";
				}
				if(defined read_file_line($fileHandles[$run])) {
					@current_line = split('\t', read_file_line($fileHandles[$run])) if defined read_file_line($fileHandles[$run]);
				} else {
					last;
				}
			}
		}
	}
	foreach my $start_pos (keys %save_peaks) {
		print OUT $save_peaks{$start_pos}->{'line'} . "\n";		
	}
	%save_peaks = ();
}

sub read_file_line {
        my $fh = shift;
        my @split;
        if ($fh and my $line = <$fh>) {
                chomp $line;
                return $line;
        }
        return;
}

sub add_mutations{
	my $start = $_[0];
	my $end = $_[1];
	$mutations = "";
	for(my $i = 0; $i < @strains; $i++) {
		for(my $aa = 0; $aa < $allele; $aa++) {
			$mutations .= "\t";
			if($strains[$i] eq "reference") { next; }
			$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains[$i] . "_allele_" . ($aa + 1) . " WHERE chr = \'" . $chr . "\' AND pos >= " . $start . " AND pos <= " . $end);
#			print "SELECT * FROM " . $genome . "_mutations_" . $strains[$i] . "_allele_" . ($aa + 1) . " WHERE chr = \'" . $chr . "\' AND pos >= " . $start . " AND pos <= " . $end . "\n";
			$sth->execute();
			while(my $mut = $sth->fetchrow_hashref()) {
				$mutations .= $mut->{'pos'} . ":" . $mut->{'reference'} . "->" . $mut->{'strain'} . ";";
			}
		}
	}
	$sth->finish();
}
$dbh->disconnect();

