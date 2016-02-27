#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Storable;
BEGIN {push @INC, '/home/vlink/mouse_strains/marge/general'}
use config;
use Data::Dumper;
use Set::IntervalTree;


$_ = 0 for my ($first, $end, $start, $current_pos, $sam, $peak, $tag, $help, $all, $iterate_chr, $SS, $print_header);
$_ = "" for my ($dir, $organism, $name, $chr, $last, $data_dir, $out_name);
$_ = () for my (@files, @strains, $dbh, $sth, @split, $c, $shift, %shift, %shift_vector, $allele, @a, %last, @shift_vector, @shift_files, @saved_file, %chr, %seen, @tag_files, %tag_files, %lookup, %lookup_chr, @split_tree, %tree_save);
my $print = 1;
$allele = 1;

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
	print STDERR "\t-dir <directory with files to shift>: has to be the same strain\n";
	print STDERR "\t-files <list with files>: comma separated list\n";
	print STDERR "\t-strains: one or several strains - comma separated\n";
	print STDERR "\t-allele: allele that was used: default: 1\n";
	print STDERR "\tIf only one strain is specified all files are shifted from this strain\n";
	print STDERR "\tIf several strains are specified it has to match the number of files specified!\n";
	print STDERR "\t-data_dir <directory>: directory with data: default specified in config\n";
	print STDERR "\n\nFormat to shift (default = sam):\n";
	print STDERR "\t-sam: shifts sam file (result file from mapping with bowtie or STAR\n";
	print STDERR "\t-homer: shifts HOMER peak and RNA files\n";
	print STDERR "\t-tag: shifts tag directory\n";
	print STDERR "\n\n";
        print STDERR "Run modes:\n";
        print STDERR "\t-all: Saves all shifting vectors per strain: needs most memory - default\n";
        print STDERR "\t-chr: shifts files per chromosomes\n";
        print STDERR "\t-SS: saves sam files when shifting perl chromosome - default: iterate over sam file (slower but less memory intensive)\n";
        print STDERR "\n";
	print STDERR "-h | --help: prints help\n";
        exit;
}

print STDERR "TODO: Different input formats need to be supported - so far only sam files!\n";	
print STDERR "Do memory calc\n";

if(@ARGV < 1) {
        &printCMD();
}

#TODO add resume then clean up code and then done
GetOptions(     "dir=s" => \$dir,
                "files=s{,}" => \@files,
                "strains=s{,}" => \@strains,
                "allele=s" => \$allele,
                "data_dir=s" => \$data_dir,
                "sam" => \$sam,
                "homer" => \$peak,
                "tag" => \$tag,
                "all" => \$all,
                "chr" => \$iterate_chr,
                "SS" => \$SS,
                "h" => \$help,
                "help" => \$help)
        or die (&printCMD());

if($help == 1) {
        &printCMD();
}


if($sam == 0 && $peak == 0 && $tag == 0) {
	$sam = 1;
}

my $split_sign;
my $shift_order;
my $chr_pos;

if($iterate_chr == 0) {
	$all = 1;
}

if(@files == 0 && $dir eq "") {
	print STDERR "Files for shifting are missing!\n";
	exit;
}

#Read in files from directory
if(@files == 0) {
	if($sam == 1) {
		@files = `ls $dir/*sam`;
	}
	if($tag == 1) {
		@files = `ls $dir/*tags.tsv`;
	}
	if($peak == 1) {
		@files = `ls $dir/*`;
	}
}

if(@strains == 0) {
	print STDERR "Strains are not specified\n";
	exit;
}
#Check if files and strains have the same length
if(@strains > 1 && @strains != @files) {
	print STDERR "Error: Length of file list and strains list is different\n";
	exit;
}

#Remove commas
for(my $i = 0; $i < @files; $i++) {
	chomp $files[$i];
	$files[$i] =~ s/,//g;
}
for(my $i = 0; $i < @strains; $i++) {
	chomp $strains[$i];
	$strains[$i] =~ s/,//g;
	$strains[$i] = uc($strains[$i]);
}

if($data_dir eq "") {
	$data_dir = config::read_config()->{'data_folder'};
}

#Choose run modes
if($all == 1 && $tag == 0) {
	my $last_strain = "";
	#Save all vectors befor shifting!
	#There is only one strain specified - all files are shifted with the same vector
	for(my $i = 0; $i < @strains; $i++) {
		if($last_strain ne $strains[$i]) {
			$last = $data_dir . "/" . $strains[$i] . "/last_shift_strain.txt"; 
			%last = %{retrieve($last)};
			$last = $data_dir . "/" . $strains[$i] . "/lookup_table_chr.txt";
			%lookup = %{retrieve($last)};
			@shift_files = `ls $data_dir/$strains[$i]/*allele_$allele*strain_to_ref.vector`;
			foreach my $shift (@shift_files) {
				chomp $shift;
				print STDERR "load " . $shift . "\n";
				chomp $shift;
				@split = split("chr", $shift);
				@split = split("_", $split[-1]);
				if(exists $lookup{$split[0]}) {
					$split[0] = $lookup{$split[0]};
				}
				my $tree = Set::IntervalTree->new;
				open SHIFT, "<$shift";
				print STDERR $split[0] . "\n";
				foreach my $line (<SHIFT>) {
					chomp $line;
					@split_tree = split('\t', $line);
					$tree->insert( {'shift'=>$split_tree[0] }, $split_tree[1], $split_tree[2]);	
				}
				close SHIFT;
				$tree_save{$split[0]} = $tree;
			}
		}
		if($sam == 1) {
			$out_name = substr($files[$i], 0, length($files[$i]) - 4) . "_shifted_from_" . $strains[$i] . ".sam";
			open OUT, ">$out_name";
			&shift_sam_file($files[$i], $strains[$i], 0);
			close OUT;
		} elsif($peak == 1) {
			$out_name = $files[$i] . "_shifted_from_" . $strains[$i] . ".txt";
			open OUT, ">$out_name";
			&shift_peak_file($files[$i], $strains[$i], 0);
			close OUT;

		}
		$last_strain = $strains[$i];
	}
} elsif($tag == 0) {
	#Save sam file
	if($SS == 1) {
		for(my $i = 0; $i < @strains; $i++) {
			@saved_file = ();
			$last = $data_dir . "/" . $strains[$i] . "/last_shift_strain.txt";
			%last = %{retrieve($last)};
			$last = $data_dir . "/" . $strains[$i] . "/lookup_table_chr.txt";
			%lookup = %{retrieve($last)};
			#Save sam file
			open FH, "<$files[$i]";
			if($sam == 1) {
				$out_name = substr($files[$i], 0, length($files[$i]) - 4) . "_shifted_from_" . $strains[$i] . ".sam";
				open OUT, ">$out_name";
				foreach my $line (<FH>) {
					chomp $line;
					if(substr($line, 0, 1) eq "@") {
						print OUT $line . "\n";
						next;
					}
					@split = split('\t', $line);
					$chr{substr($split[2], 3)} = 1;
					push(@{$saved_file[substr($split[2], 3)]}, $line);	
				}
				close FH;
			} elsif($peak == 1) {
				print STDERR "save peak files!\n";
				$out_name = $files[$i]. "_shifted_from_" . $strains[$i] . ".txt";
				open OUT, ">$out_name";
				foreach my $line (<FH>) {
					chomp $line;
					if(substr($line, 0, 1) eq "#") {
						print OUT $line . "\n";
					}
					@split = split('\t', $line);
					$chr{substr($split[1], 3)} = 1;
					push(@{$saved_file[substr($split[1], 3)]}, $line);
				}
			}
			my $first = 0;
			foreach my $c (keys %chr) {
				$shift = $data_dir . "/" . $strains[$i] . "/chr" . $c . "_allele_" . $allele . ".strain_to_ref.vector";
				if(exists $lookup{$c}) {
					$c = $lookup{$c};
				}
				my $tree = Set::IntervalTree->new;
				open SHIFT, "<$shift";
				foreach my $line (<SHIFT>) {
					chomp $line;
					@split_tree = split('\t', $line);
					$tree->insert( {'shift'=>$split_tree[0] }, $split_tree[1], $split_tree[2]);
				}
				close SHIFT;
				$tree_save{$c} = $tree;		
				if($sam == 1) {
					&shift_sam_file_saved($c, $first);
				} elsif($peak == 1) {
					&shift_peak_file_saved($c, $first);
				}
				$first++;
			}
			close OUT;
		}
	} else {
		for(my $i = 0; $i < @strains; $i++) {
			$last = $data_dir . "/" . $strains[$i] . "/last_shift_strain.txt";
			%last = %{retrieve($last)};
			$last = $data_dir . "/" . $strains[$i] . "/lookup_table_chr.txt";
			%lookup = %{retrieve($last)};
			#iterate over sam file
			if($sam == 1) {
				$out_name = substr($files[$i], 0, length($files[$i]) - 4) . "_shifted_from_" . $strains[$i] . ".sam";
			} elsif($peak == 1) {
				$out_name = $files[$i] . "_shifted_from_" . $strains[$i] . ".txt";
			}
			open OUT, ">$out_name";
			@shift_files = `ls $data_dir/$strains[$i]/*allele_$allele*strain_to_ref.vector`;
			foreach my $shift (@shift_files) {
				chomp $shift;
				@split = split("chr", $shift);
				@split = split("_", $split[-1]);
				if(exists $lookup{$split[0]}) {
					$split[0] = $lookup{$split[0]};
				}
				my $tree = Set::IntervalTree->new;
				open SHIFT, "<$shift";
				foreach my $line (<SHIFT>) {
					chomp $line;
					@split_tree = split('\t', $line);
					$tree->insert( {'shift'=>$split_tree[0] }, $split_tree[1], $split_tree[2]);
				}
				close SHIFT;
				$tree_save{$split[0]} = $tree;		
				if($sam == 1) {
					&shift_sam_file($files[$i], $strains[$i], $split[0]);
				} elsif($peak == 1) {
					&shift_peak_file($files[$i], $strains[$i], $split[0]);
				}
			}
			close OUT;
		}
	}
} else {
	for(my $i = 0; $i < @files; $i++) {
		chomp $files[$i];
		@tag_files = `ls $files[$i]/*tags.tsv`;		
		for(my $j = 0; $j < @tag_files; $j++) {
			chomp $tag_files[$j];
			$chr = substr($split[-1], 3, length($split[-1]) - 12);
			@split = split("/", $tag_files[$j]);
			$tag_files{$chr}{$tag_files[$j]} = 1;
		}
	}
	print Dumper %tag_files;

}

#Shifts position
sub shift{
	my $chr = $_[0];
	my $chr_num = $chr;
	my $fetch;
	if(exists $lookup{$chr}) {
		$chr_num = $lookup{$chr};
	}
	my $allele = $_[1];
	my $pos_to_shift = $_[2];
	my $pos_shifted;
	if($pos_to_shift > $last{$chr}{$allele}{'pos'}) {
		$pos_shifted = $pos_to_shift + $last{$chr}{$allele}{'shift'};
	} else {
		$fetch = $tree_save{$chr_num}->fetch($pos_to_shift, $pos_to_shift + 1);
		if(scalar(@$fetch) == 0) {
			$pos_shifted = $pos_to_shift;
		} else {
			$pos_shifted = $pos_to_shift + $fetch->[0]->{'shift'};
		}
	}
	return $pos_shifted;
}


sub shift_sam_file_saved{
	#Shift files when all or one chromosomes are saved in memory and sam file is not
	my $chr = $_[0];
	my $first = $_[1];
	print $chr . "\t" . $first . "\n";
	foreach my $line (@{$saved_file[$chr]}) {
		chomp $line;
		@split = split('\t', $line);
		if(!exists $last{$chr}{$allele}) {
			print OUT $line . "\n";
			next;
		}
		print OUT $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . &shift($chr, $allele, $split[3]);
		for(my $j = 4; $j < @split; $j++) {
			print OUT "\t" . $split[$j];
		}
		print OUT "\n";
	}
}

sub shift_peak_file_saved{
	my $chr = $_[0];
	my $first = $_[1];
	foreach my $line (@{$saved_file[$chr]}) {
		chomp $line;
		@split = split('\t', $line);
		if(!exists $last{$chr}{$allele}) {
			print OUT $line . "\n";
			next;
		}
		print OUT $split[0] . "\t" . $split[1] . "\t" . &shift($chr, $allele, $split[2]) . "\t" . &shift($chr, $allele, $split[3]);
		for(my $j = 4; $j < @split; $j++) {
			print OUT "\t" . $split[$j];
		}	
		print OUT "\n";
	}
}


sub shift_sam_file{
	#Shift files when all or one chromosomes are saved in memory and sam file is not
	my $file = $_[0];
	my $strain = $_[1];
	my $chr_to_work_on = $_[2];
	open FH, "<$file";
	print STDERR "Reading in " . $file . "\n";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if(substr($line, 0, 1) eq "@" && $print_header == 0) {
			print OUT $line . "\n";
			next;
		}
		$print_header++;
		if(length($split[2]) < 4) {
			next;
		}
		$chr = substr($split[2], 3);
		if(!exists $last{$chr}{$allele}) {
			if($chr_to_work_on == 0) {
				print OUT $line . "\n";
			} elsif(!exists $seen{$split[0]}) {
				print OUT $line . "\n";
				$seen{$split[0]} = 1;
			}
			next;
		}
		if($all == 0 && $chr != $chr_to_work_on) { next; }
		print OUT $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . &shift($chr, $allele, $split[3]);
		for(my $j = 4; $j < @split; $j++) {
			print OUT "\t" . $split[$j];
		}
		print OUT "\n";
	}
}

sub shift_peak_file{
	my $file = $_[0];
	my $strain = $_[1];
	my $chr_to_work_on = $_[2];
	open FH, "<$file";
	foreach my $line (<FH>) {
		chomp $line;
		@split = split('\t', $line);
		if(!exists $last{$chr}{$allele}) {
			if($chr_to_work_on == 0) {
				print OUT $line . "\n";
			} elsif(!exists $seen{$split[0]}) {
				print OUT $line . "\n";	
				$seen{$split[0]} = 1;
			}
			next;
		}
		if($all == 0 && $chr != $chr_to_work_on) { next; }
		print OUT $split[0] . "\t" . $split[1] . "\t" . &shift($chr, $allele, $split[2]) . "\t" . &shift($chr, $allele, $split[3]);
		for(my $j = 4; $j < @split; $j++) {
			print OUT "\t" . $split[$j];
		}
		print OUT "\n";
	}
}
