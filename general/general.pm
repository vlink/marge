#!/usr/bin/perl -w

package general;

use Storable;
use Set::IntervalTree;

sub read_strains_data {
	my $strain = $_[0];
	my $data_dir = $_[1];	
	my $allele = $_[2];
	my $shift_direction = $_[3];
	print $strain . "\n";
	print $data_dir . "\n";
	$_ = () for my (%last, %lookup, %tree_strain, @split, @shift_files, @split_tree);
	my $file;
	if($shift_direction eq "ref_to_strain") {
		$file = $data_dir . "/" . $strain . "/last_shift_ref.txt";
        	%last = %{retrieve($file)};
	} else {
		$file = $data_dir . "/" . $strain . "/last_shift_strain.txt";
		%last = %{retrieve($file)};
	}
	$file = $data_dir . "/" . $strain . "/lookup_table_chr.txt";
        %lookup = %{retrieve($file)};
        @shift_files = `ls $data_dir/$strain/*$shift_direction*`;
        foreach my $shift (@shift_files) {
                chomp $shift;
		open FH,"<$shift";
		@split = split("chr", $shift);
		@split = split("_", $split[-1]);
		if(exists $lookup{$split[0]}) {
			$split[0] = $lookup{$split[0]};
		}
		my $tree = Set::IntervalTree->new;
		open SHIFT, "<$shift";
		foreach my $line (<SHIFT>) {
			print $line . "\n";
			chomp $line;
			@split_tree = split('\t', $line);
			$tree->insert( {'shift'=>$split_tree[0] }, $split_tree[1], $split_tree[2]);
		}
		close SHIFT;
		$tree_strain{$split[0]} = $tree;
        }
	return (\%tree_strain, \%last, \%lookup);	
}

1;

