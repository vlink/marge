#!/usr/bin/perl -w

package general;

use Storable;
use Set::IntervalTree;

#Read in the lookup and last files for every strain
sub read_strains_data {
	my $strain = $_[0];
	my $data_dir = $_[1];	
	my $allele = $_[2];
	my $shift_direction = $_[3];
	$_ = () for my (%last, %lookup, %tree_strain, @split, @shift_files, @split_tree);
	my $file;
	#Check which shfit direction is specified
	if($shift_direction eq "ref_to_strain") {
		$file = $data_dir . "/" . $strain . "/last_shift_ref.txt";
		if(-e $file) {
			%last = %{retrieve($file)};
		}
	} else {
		$file = $data_dir . "/" . $strain . "/last_shift_strain.txt";
		if(-e $file) {
			%last = %{retrieve($file)};
		}
	}
	$file = $data_dir . "/" . $strain . "/lookup_table_chr.txt";
	if(-e $file) {
        	%lookup = %{retrieve($file)};
	}
	#Read in all shift files and store in Interval Tree
       	@shift_files = `ls $data_dir/$strain/*$shift_direction* 2> /dev/null`;
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
			$tree->insert( {'shift'=>$split_tree[0] }, $split_tree[1], $split_tree[2]+1);
		}
		close SHIFT;
		$tree_strain{$split[0]} = $tree;
        }
	return (\%tree_strain, \%last, \%lookup);	
}

sub wait_10_secs{
	my $file = $_[0];
	print STDERR "$file exists already!\n";
	print STDERR "Waiting for 10 seconds before overwriting!\n";
	print STDERR "Press Ctrl + C to abort\n";
	print STDERR "Wait";
	for(my $i = 0; $i < 10; $i++) {
		print STDERR ".";
		sleep(1);
	}
	print STDERR "\n";
}
1;

