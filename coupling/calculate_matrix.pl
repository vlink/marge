#!/usr/bin/perl -w
#
use strict;
use Statistics::Basic qw(:all);
use threads;
use Getopt::Long;

$_ = () for my (@files, %save_one, %save_two, @split, @vector, $v, %count, %name_one, %name_two, %data_for_threading, @running, @Threads, $current_thread, $current_thread_level);
$_ = 0 for my ($samples, $core, $first, $two, $s, $count, $all);
$_ = "" for my($output);

sub printCMD {
	print STDERR "\t-files <comma separated list of files> (if one file is specified thecorrelation of the mark is calculated, if two files are specified two marks are considered)\n";
	print STDERR "\t-name <name for output> (default matrix_chr.txt)\n";
	print STDERR "\t-samples <number of samples> (number of samples in the file\n";
	print STDERR "\t-core <number of cores for threading> (Default: 12 cores)\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

GetOptions(	"files=s{,}" => \@files,
		"name=s" => \$output,
		"samples=s" => \$samples, 
		"core=s" => \$core)
	or die ("Error in command line options!\n");

if(@files < 1 || $output eq "" || $samples == 0) {
	&printCMD();
}

if(@files == 2) {
	$two = 1;
}
for(my $i = 0; $i < @files; $i++) {
	$files[$i] =~ s/,//g;
	$files[$i] =~ s/ //g;
}
print STDERR "Read in file number one: " . $files[0] . "\n";
$core = 12;

open(my $fh, "<", $files[0]) or die ("Can't open file $files[0]\n");
while( my $line = <$fh>) {
	chomp $line;
	if($first == 0) { $first++; next; }
	@split = split('\t', $line);
	$count{$split[1]}++;
	@vector = ();
	for(my $i = @split - $samples; $i < @split; $i++) {
		push @vector, $split[$i];
	}
	$v = vector(@vector);
	$save_one{$split[1]}{$split[2]} = $v;
	$name_one{$split[1]}{$split[2]} = $split[0];
}
close $fh;

if($two == 1) {
	$first = 0;
	print STDERR "Read in file number two: " . $files[1] . "\n";
	open(my $fh, "<", $files[1]) or die ("Can't open file $files[1]\n");
	while( my $line = <$fh>) {
		chomp $line;
		if($first == 0) { $first++; next; }
		@split = split('\t', $line);
	#	$count{$split[1]}++;
		@vector = ();
		for(my $i = @split - $samples; $i < @split; $i++) {
			push @vector, $split[$i];
		}
		$v = vector(@vector);
		$save_two{$split[1]}{$split[2]} = $v;
		$name_two{$split[1]}{$split[2]} = $split[0];
	}
	close $fh;
}
print STDERR "done\n";

$all = keys %save_one;
$count = 1;
foreach my $chr (keys %save_one) {
	$data_for_threading{$chr} = 1;
}

while(scalar @Threads < keys %save_one) {
	@running = threads->list(threads::running);
	if(scalar @running < $core) {
		foreach my $key (keys %data_for_threading) {
			$current_thread = $key;
		}
		delete $data_for_threading{$current_thread};
		my $thread = threads->new( sub { thread_routine($current_thread, $two) });
		push(@Threads, $thread);
		my $tid = $thread->tid;
		@running = threads->list(threads::running);
	}
	@running = threads->list(threads::running);
	foreach my $thr (@Threads) {
		if($thr->is_running()) {
			my $tid = $thr->tid;
		} elsif($thr->is_joinable()) {
			my $tid = $thr->tid;
			$thr->join;
		}
	}
	@running = threads->list(threads::running);
}
while(scalar @running > 0) {
	foreach my $thr (@Threads) {
		$thr->join if ($thr->is_joinable());
	}
	@running = threads->list(threads::running);
}

sub thread_routine {
	my $chr = $_[0];
	my $two = $_[1];
	my $o;
	my $calculated_correlation = 0;
	if($output eq "") {
		$o = "matrix_" . $chr . ".txt";
	} else {
		$o = $output . "_" . $chr . ".txt";
	}
	open OUT, ">$o";
	foreach my $key (sort {$a <=> $b} keys %{$save_one{$chr}}) {
		print OUT "\t" . $name_one{$chr}{$key};
	}
	print OUT "\n";
	foreach my $key (sort {$a <=> $b} keys %{$save_one{$chr}}) {
		print OUT $name_one{$chr}{$key};
		if($two == 0) {
			foreach my $key2 (sort {$a <=> $b} keys %{$save_one{$chr}}) {
				$calculated_correlation = correlation($save_one{$chr}{$key}, $save_one{$chr}{$key2});
				if($calculated_correlation eq "n/a") {
					$calculated_correlation = 0;		
				}
				print OUT "\t" . $calculated_correlation;
			}
		} else {
			foreach my $key2 (sort {$a <=> $b} keys %{$save_one{$chr}}) {
				$calculated_correlation = correlation($save_one{$chr}{$key}, $save_two{$chr}{$key2}); 
				if($calculated_correlation eq "n/a") {
					$calculated_correlation = 0;
				}
				print OUT "\t" . $calculated_correlation;
			}

		}
		print OUT "\n";
	}
	close OUT;
	print STDERR "Thread is done for: " . $chr . "\n";
}
