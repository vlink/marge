#!/usr/bin/perl -w
#
use strict;
use Getopt::Long;
use Storable;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/db_part'};
use processing;
use Set::IntervalTree;
use Data::Dumper;

$_ = 0 for my($hetero, $num_strains, $diff, $allele);
$_ = () for my(@strains, @split, %last_shift_pos_strain, %last_shift_pos_ref, %lookup_table, @tmp, @split_last, @tmp2);
$_ = "" for my($data, $last_line_strains, $last_line_ref, $command, $out);


sub printCMD{
	print STDERR "This script updates the last_shift and lookup tables\n";
	print STDERR "It should be used when a strains was generated out of different files per chromosome, with the parameter -add\n";
	print STDERR "Then the last_shift and lookup tables just contain the last entry, so it needs to be updated and all chromosomes need to be added\n\n\n";
	print STDERR "Usage:\n";
	print STDERR "-strains <list of strains>\n";
	print STDERR "-hetero: Strain is heterozygouse (Default: homozygouse)\n";
	print STDERR "-dir: Directory where mutation files are located - default: folder specified in confing file\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-strains' => 1);
my %convert = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%convert);

#Read in command line arguments
GetOptions(	"-strains=s{,}" => \@strains,
		"-hetero" => \$hetero,
		"-dir=s" => \$data)
or die(&printCMD());

if($data eq "") {
        $data = config::read_config()->{'data_folder'};
}

$num_strains = @strains;
for(my $i = 0; $i < @strains; $i++) {
        $strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

for(my $i = 0; $i < @strains; $i++) {
	my @files = `ls $data/$strains[$i]/*mut`;

	foreach my $f (@files) {
		print $f;
		chomp $f;
		@tmp2 = split("/", $f);
		@tmp = split("_", $tmp2[-1]);
		$tmp[0] = substr($tmp[0], 3);
		#First step: get the last diff
		open FH, "<$f";
		$allele = substr($f, length($f) - 5, 1);
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			if(length($split[1]) != length($split[2])) {
				$diff = length($split[2]) - length($split[1]);
			}
		}
		#Second step: Get the last line in ref to strain
		$command = "tail -n1 " . substr($f, 0, length($f) - 4) . ".ref_to_strain.vector";
		$last_line_ref = `$command`;
		chomp $last_line_ref;
		@split_last = split('\t', $last_line_ref);
		$last_shift_pos_strain{$tmp[0]}{$allele}{'pos'} = ($split_last[-1]+1);
		$last_shift_pos_strain{$tmp[0]}{$allele}{'shift'} = $split_last[0] + $diff;
		#Third step: Get the last line in strain to ref
		$command = "tail -n1 " . substr($f, 0, length($f) - 4) . ".strain_to_ref.vector";
		$last_line_ref = `$command`;
		chomp $last_line_ref;
		@split_last = split('\t', $last_line_ref);
		$last_shift_pos_ref{$tmp[0]}{$allele}{'pos'} = $split_last[1];
		$last_shift_pos_ref{$tmp[0]}{$allele}{'shift'} = $split_last[0];
	}
	$out = $data . "/" . $strains[$i] . "/last_shift_strain.txt";
        store \%last_shift_pos_ref, "$out";
	$out = $data . "/" . $strains[$i] . "/last_shift_ref.txt";
        store \%last_shift_pos_strain, "$out";
	%last_shift_pos_strain = ();
	%last_shift_pos_ref = ();

}
