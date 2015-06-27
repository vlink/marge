#!/usr/bin/perl -w

use strict;
#require '/data/home/vlink/mouse_strains/database/database_interaction.pm';
require 'database_interaction.pm';
use Getopt::Long;
require '../general/config.pm';
use Data::Dumper;

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
	print STDERR "\t-method <gene|protein|both|genomic|align_gene|align_protein|align_both|align_genomic_gene|align_genomic_protein|align_genomic_both|genomic_protein|genomic_both>\n";
        print STDERR "\t-gene <gene name>\n";
        print STDERR "\t-transcript <RefSeq transcript name>\n";
        print STDERR "\t-strains <strains>: Comma-separated list of strains - if not defined only the reference is used - if all is specified all strains are used\n";
        print STDERR "\t-start <pos>: Start position of genomic region\n";
        print STDERR "\t-end <pos>: Stop position of genomic region\n";
        print STDERR "\t-chr <chromosome>:Chromosome for genomic region\n";
	print STDERR "\t-strand <+|->: Default: +\n";
	print STDERR "\t-homo: Data is homozygouse\n";
	print STDERR "\t-html: Outputs a html file - alignment mismatches are color in red\n";
	print STDERR "\n\nInterspecies comparison\n";
	print STDERR "\t-genome2: Genome for second species\n";
	print STDERR "\t-strains2 <strains>: Comma-separated list of strains for second species\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

$_ = "" for my ($genome, $gene, $transcript, $chr, $method, $strand, $genome2);
$_ = 0 for my ($start, $end, $homo, $html);
my (%strains, $dbh, $sth, @strains, @strains2);

my %mandatory = ('-method' => 1, '-genome' => 1);
my %convert = map { $_ => 1 } @ARGV;

config::check_parameters(\%mandatory, \%convert);

GetOptions(	"genome=s" => \$genome,
		"method=s" => \$method,
		"gene=s" => \$gene,
		"transcript=s" => \$transcript,
		"strains=s{,}" => \@strains,
		"start=s" => \$start,
		"end=s" => \$end,
		"chr=s"=> \$chr,
		"homo" => \$homo,
		"strand=s" => \$strand, 
		"genome2=s" => \$genome2,
		"strains2=s{,}" => \@strains2,
		"html" => \$html)
	or die("Error in command line options!\n");


if($strand eq "") {
	$strand = "+";
}

if(@strains < 1) {
	$strains[0] = "reference";
}

if($gene ne "") {
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, $gene);
} elsif($method eq "genomic" || $method eq "align_genomic") {
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, "genomic");
} else {
	database_interaction::set_global_variables(\@strains, $genome, $homo, $html, $transcript);
}

#Strains 1
my $ref_as_strain = 0;
for(my $i = 0; $i < @strains; $i++) {
	chomp $strains[$i];
	$strains[$i] =~ s/,//g;
	if($strains[$i] eq "reference") {
		$ref_as_strain = 1;
	}
}
if(@strains > 0 && $strains[0] eq "all") {
	my $hash = config::get_strains($genome);
	my $i = 0;
	foreach my $key (keys %{$hash}) {
		$strains[$i] = $key;
		$i++;
		if($key eq "reference") {
			$ref_as_strain = 1;
		}
	}
}
if($ref_as_strain == 0) {
	push(@strains, "reference");
}

if($start == 0 && $end == 0 && $gene eq "" && $transcript eq "") {
	print STDERR "Please specify gene or transcript\n";
	exit;
}

if($method =~ /genomic/) {
	if(($start == 0 && $end == 0) || $chr eq "") {
		print STDERR "start, end or chromosome is not specified!\n";
		exit;
	} 
}

if($genome2 eq "") {
	if($method eq "gene") {
		database_interaction::get_gene($gene, $transcript, 1, 0, 0, 0);
	} elsif($method eq "protein") {
		database_interaction::get_gene($gene, $transcript, 0, 0, 1, 0);
	} elsif($method eq "both") {
		database_interaction::get_gene($gene, $transcript, 1, 0, 1, 0);
	} elsif($method eq "genomic") {
		database_interaction::get_genomic_seq($start, $end, $chr, $strand, 1, 0, 0, 0);	
	} elsif($method eq "align_genomic_gene") {
		database_interaction::get_genomic_seq($start, $end, $chr, $strand, 1, 1, 0, 0);	
	} elsif($method eq "align_gene") {
		database_interaction::get_gene($gene, $transcript, 1, 1, 0, 0);
	} elsif($method eq "align_protein") {
		database_interaction::get_gene($gene, $transcript, 0, 0, 1, 1);
	} elsif($method eq "align_both") {
		database_interaction::get_gene($gene, $transcript, 1, 1, 1, 1);
	} elsif($method eq "align_genomic_protein") {
		database_interaction::get_genomic_seq($start, $end, $chr, $strand, 0, 0, 0, 1);
	} elsif($method eq "align_genomic_both") {
		database_interaction::get_genomic_seq($start, $end, $chr, $strand, 0, 1, 0, 1);
	} elsif($method eq "genomic_protein") {
		database_interaction::get_genomic_seq($start, $end, $chr, $strand, 0, 0, 1, 0);
	} elsif($method eq "genomic_both") {
		database_interaction::get_genomic_seq($start, $end, $chr, $strand, 1, 0, 1, 0);
	} else {
		print STDERR "Method unknown!\n";
		exit;
	}
} else {
	#Strains 2
	my $ref_as_strain = 0;
	for(my $i = 0; $i < @strains2; $i++) {
		chomp $strains2[$i];
		$strains2[$i] =~ s/,//g;
		if($strains2[$i] eq "reference") {
			$ref_as_strain = 1;
		}
	}
	if(@strains2 > 0 && $strains2[0] eq "all") {
		my $hash = config::get_strains($genome2);
		my $i = 0;
		foreach my $key (keys %{$hash}) {
			$strains2[$i] = $key;
			$i++;
			if($key eq "reference") {
				$ref_as_strain = 1;
			}
		}
	}
	if($ref_as_strain == 0) {
		push(@strains2, "reference");
	}
	if($method eq "gene") {
		database_interaction::get_gene_interspecies($gene, $transcript, $genome2, \@strains2, 1, 0, 0, 0);
	} elsif($method eq "align_gene") {
		database_interaction::get_gene_interspecies($gene, $transcript, $genome2, \@strains2, 1, 1, 0, 0);
	} elsif($method eq "protein") {
		database_interaction::get_gene_interspecies($gene, $transcript, $genome2, \@strains2, 0, 0, 1, 0);
        } elsif($method eq "align_protein") {
                database_interaction::get_gene_interspecies($gene, $transcript, $genome2, \@strains2, 0, 0, 1, 1);
	} elsif($method eq "both") {
		database_interaction::get_gene_interspecies($gene, $transcript, $genome2, \@strains2, 1, 0, 1, 0);
	} elsif($method eq "align_both") {
		database_interaction::get_gene_interspecies($gene, $transcript, $genome2, \@strains2, 1, 1, 1, 1);
	} elsif($method eq "genomic" || $method eq "align_genomic_gene") {
		database_interaction::get_genomic_seq_interspecies($start, $end, $chr, $strand, $genome2, \@strains2, 1, 1, 0, 0);
	} elsif($method eq "genomic_protein" || $method eq "align_genomic_protein") {
		database_interaction::get_genomic_seq_interspecies($start, $end, $chr, $strand, $genome2, \@strains2, 0, 0, 1, 1);
	} elsif($method eq "genomic_both" || $method eq "align_genomic_both") {
		database_interaction::get_genomic_seq_interspecies($start, $end, $chr, $strand, $genome2, \@strains2, 1, 1, 1, 1);
	} else {
		print STDERR "Method unknown!\n";
		exit;
	}

}
