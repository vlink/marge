#!/usr/bin/perl -w

package processing;

sub create_genome {
	my $chr = $_[0];
	my $output = $_[1];
	my $input = $_[2];
	my $mutation_file = $_[3];
	open FH, "<$input";
	#First save the genome in an array
	my @cache = ();
	my @split;
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			next;
		}
		@split = split('', $line);
		foreach my $s (@split) {
			push(@cache, $s);
		}
	}

	#Now open the mutation file and change the mutated sequences
	my @split_mut;
	open FH, "<$mutation_file";
	foreach my $line (<FH>) {
		chomp $line;
		@split_mut = split('\t', $line);
		if($split_mut[0] > @cache) { last; }
		if(length($split_mut[1]) > 1) {
			for(my $i = 1; $i < length($split_mut[1]); $i++) {
				$cache[$split_mut[0] + $i - 1] = "";
			}
		} else {
			$cache[$split_mut[0] - 1] = $split_mut[2];
		}
	}
	close FH;

	#Write output
	my $strain_seq = join("", @cache);
	open OUT, ">$output";
	print OUT ">chr" . $chr . "\n";
	while(length($strain_seq) > 50) {
		print OUT substr($strain_seq, 0, 50) . "\n";
		$strain_seq =~ s/^.{50}//s; 
	}
	print OUT $strain_seq . "\n";
	close OUT;
}
1;
