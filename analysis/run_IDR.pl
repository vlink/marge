#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use threads;

$_ = "" for my($output_dir, $command);
$_ = 0 for my($core, $job_id);
$_ = () for my(@tag_dirs, @input_dirs, @split);
my $idr_path = "/gpfs/data01/glasslab/home/jtao/software/anaconda3/bin/idr";
 
sub printCMD{
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
        print STDERR "\t-tag_dirs <list of tag directories>\n";
	print STDERR "\t-input_dirs <list of input tag directores>\n";
	print STDERR "\t-output_dir <name of output director>\n";
	print STDERR "\t-core <number of cores> (default: 4)\n";
        exit;
}


if(@ARGV < 1) {
        &printCMD();
}

GetOptions(     "tag_dirs=s{,}" => \@tag_dirs,
                "input_dirs=s{,}" => \@input_dirs,
                "output_dir=s" => \$output_dir, 
		"core=s" => \$core)
        or die (&printCMD());

if($core == 0) { $core = 4; }
#Check if everything is defined
if(@tag_dirs == 0 || @input_dirs == 0 || $output_dir eq "") {
	&printCMD();
}

if(@tag_dirs == 1) {
	@split = split(' ', $tag_dirs[0]);
	if(@split > @tag_dirs) {
		@tag_dirs = @split;
	}
}
if(@input_dirs == 1) {
	@split = split(' ', $input_dirs[0]);
	if(@split > @input_dirs) {
		@input_dirs = @split;
	}
}
if(-e $output_dir) {
	print STDERR "Output directory exists\n";
	exit;
}
`mkdir $output_dir`;

#First step - pool input directories
$command = "makeTagDirectory " . $output_dir . "/pooled_input -d ";
for(my $i = 0; $i < @input_dirs; $i++) {
	$command .= $input_dirs[$i] . " ";
}	
print STDERR "Pooling input directores\n";
#print $command . "\n";
`$command 2>> $output_dir/error.txt`;


#Second step - pool tag directories and create pseudo replicates per tag directory
$command = "makeTagDirectory " . $output_dir . "/pooled_tag_dirs -d ";
for(my $i = 0; $i < @tag_dirs; $i++) {
	$command .= $tag_dirs[$i] . " ";
	print STDERR "Generate pseudo replicates for " . $tag_dirs[$i] . "\n";
	&generate_pseudo_replicates($tag_dirs[$i]);
}
#print $command . "\n";
print STDERR "Pool all tag directories\n";
`$command 2>> $output_dir/error.txt`;

#Generate pseudo replicates for pooled tag directory
print STDERR "Generate pseudo replicates for " . $output_dir . "/pooled_tag_dirs\n";
&generate_pseudo_replicates($output_dir . "/pooled_tag_dirs");

#Call peaks 
print STDERR "Call peaks on\n";
for(my $i = 0; $i < @tag_dirs; $i++) {
	@split = split("/", $tag_dirs[$i]);
	$command = "findPeaks " . $tag_dirs[$i] . " -L 0 -C 0 -F 0 -fdr 0.9 -i " . $output_dir . "/pooled_input > " . $output_dir . "/peaks_" . $split[-1] . ".txt 2>> $output_dir/error.txt";
#	print $command . "\n";
	print STDERR "\t" . $split[-1] . "\n";
	`$command`;
	$command = "findPeaks " . $output_dir . "/" . $split[-1] . "_pseudo_1 -L 0 -C 0 -F 0 -fdr 0.9 -i " . $output_dir . "/pooled_input > " . $output_dir . "/peaks_" . $split[-1] . "_pseudo_1.txt 2>> $output_dir/error.txt";
#	print $command . "\n";
	print STDERR "\t" . $split[-1] . "_pseudo_1\n";
	`$command`;
	$command = "findPeaks " . $output_dir . "/" . $split[-1] . "_pseudo_2 -L 0 -C 0 -F 0 -fdr 0.9 -i " . $output_dir . "/pooled_input > " . $output_dir . "/peaks_" . $split[-1] . "_pseudo_2.txt 2>> $output_dir/error.txt";
#	print $command . "\n";
	print STDERR "\t" . $split[-1] . "_pseudo_2\n";
	`$command`;
}
$command = "findPeaks " . $output_dir . "/pooled_tag_dirs -L 0 -C 0 -F 0 -fdr 0.9 -i " . $output_dir . "/pooled_input > " . $output_dir . "/peaks_pooled_tag_dirs.txt 2>> $output_dir/error.txt";
#print $command . "\n";
print STDERR "\tpooled_tag_dirs\n";
`$command`;

$command = "findPeaks " . $output_dir . "/pooled_tag_dirs_pseudo_1 -L 0 -C 0 -F 0 -fdr 0.9 -i " . $output_dir . "/pooled_input > " . $output_dir . "/peaks_pooled_tag_dirs_pseudo_1.txt 2>> $output_dir/error.txt";
#print $command . "\n";
print STDERR "\tpooled_tag_dirs_pseudo_1\n";
`$command`;

$command = "findPeaks " . $output_dir . "/pooled_tag_dirs_pseudo_2 -L 0 -C 0 -F 0 -fdr 0.9 -i " . $output_dir . "/pooled_input > " . $output_dir . "/peaks_pooled_tag_dirs_pseudo_2.txt 2>> $output_dir/error.txt";
#print $command . "\n";
print STDERR "\tpooled_tag_dirs_pseudo_2\n";
`$command`;

#Convert all peak files to narrowPeak format
for(my $i = 0; $i < @tag_dirs; $i++) {
	@split = split("/", $tag_dirs[$i]);
	&convert_to_narrowPeaks($output_dir . "/peaks_" . $split[-1] . ".txt");	
	&convert_to_narrowPeaks($output_dir . "/peaks_" . $split[-1] . "_pseudo_1.txt");
	&convert_to_narrowPeaks($output_dir . "/peaks_" . $split[-1] . "_pseudo_2.txt");	
}
&convert_to_narrowPeaks($output_dir . "/peaks_pooled_tag_dirs.txt");
&convert_to_narrowPeaks($output_dir . "/peaks_pooled_tag_dirs_pseudo_1.txt");
&convert_to_narrowPeaks($output_dir . "/peaks_pooled_tag_dirs_pseudo_2.txt");

#Run IDR
#Multithreading might be useful here
#Run IDR on pooled pseudo replciates
my %commands_for_threading;

#Run IDR on pooled pseudo replicates
$command = $idr_path . " --samples " . $output_dir . "/peaks_pooled_tag_dirs_pseudo_1.narrowPeak.txt " . $output_dir . "/peaks_pooled_tag_dirs_pseudo_2.narrowPeak.txt --idr-threshold 0.01 --plot --output-file " . $output_dir . "/idr_pooled_pseudo_replicates.txt";
#print $command . "\n";
$commands_for_threading{$command} = $job_id;
$job_id++;

#Run IDR on pseudo replicates
for(my $i = 0; $i < @tag_dirs; $i++) {
	@split = split("/", $tag_dirs[$i]);
	$command = $idr_path . " --samples " . $output_dir . "/peaks_" . $split[-1] . "_pseudo_1.narrowPeak.txt " . $output_dir . "/peaks_" . $split[-1] . "_pseudo_2.narrowPeak.txt --idr-threshold 0.05 --plot --output-file " . $output_dir . "/idr_" . $split[-1] . ".txt";
#	print $command . "\n";  
	$commands_for_threading{$command} = $job_id;
	$job_id++;
}
$command = $idr_path . " --samples " . $output_dir . "/peaks_" . $split[-1] . "_pseudo_1.narrowPeak.txt " . $output_dir . "/peaks_" . $split[-1] . "_pseudo_2.narrowPeak.txt --idr-threshold 0.05 --plot --output-file " . $output_dir . "/idr_" . $split[-1] . ".txt";

my @split2;
#Run IDR or all replicates versus each other
for(my $i = 0; $i < @tag_dirs - 1; $i++) {
	@split = split("/", $tag_dirs[$i]);
	for(my $j = $i + 1; $j < @tag_dirs; $j++) {
		@split2 = split("/", $tag_dirs[$j]);
		$command = $idr_path . " --samples " . $output_dir . "/peaks_" . $split[-1] . ".narrowPeak.txt " . $output_dir . "/peaks_" . $split2[-1] . ".narrowPeak.txt --idr-threshold 0.05 --plot --output-file " . $output_dir . "/idr_" . $split[-1] . "_vs_" . $split2[-1] . ".txt";
#		print $command . "\n";
		$commands_for_threading{$command} = $job_id;
		$job_id++;
	}
}
$job_id++;

my $commands_to_run = keys %commands_for_threading;

if($core > 1) {
	$_ = () for my (%data_for_threading, @running, @Threads, $current_thread, $current_thread_level1, $current_thread_level2);
#	#Run through queue and write strains data
	while (scalar @Threads < $commands_to_run) {
		@running = threads->list(threads::running);
		if(scalar @running < $core) {
			foreach my $keys (sort {$commands_for_threading{$b} <=> $commands_for_threading{$a} } keys %commands_for_threading) {
					$current_thread = $keys;
					$job_id = $commands_for_threading{$keys};
			}
			delete $commands_for_threading{$current_thread};
			#Start new thread and execute sub-function thread-routine
			my $thread = threads->new( sub { thread_routine($current_thread, $job_id) });
			push (@Threads, $thread);
			my $tid = $thread->tid;
			@running = threads->list(threads::running);
		}
		@running = threads->list(threads::running);
		foreach my $thr (@Threads) {
			if($thr->is_running()){
				my $tid = $thr->tid;
			} elsif($thr->is_joinable()){
				my $tid = $thr->tid;
				$thr->join;
			}
		}
		@running = threads->list(threads::running);
	}
	#Wait for the remaining threads to finish
	while (scalar @running > 0) {
		foreach my $thr (@Threads) {
			$thr->join if ($thr->is_joinable());
		}
		@running = threads->list(threads::running);
	}
} else {
	print STDERR "Running IDR without multi threading\n";
	foreach my $key (keys %commands_for_threading) {
		print STDERR $key . "\n";
		`$key`;
	}
}

print STDERR "IDR analysis finished\n";
print STDERR "Interpretation of results\n";
#Look at consistency of different replicates
my %orig_replicates;
my $num;
my $min_num = 10000000;
my $max_num = 0;
my $min_comp = "";
my $max_comp = "";
my @split_num;
open FLAG, "> $output_dir/flag_experiments.txt";

for(my $i = 0; $i < @tag_dirs - 1; $i++) {
	@split = split("/", $tag_dirs[$i]);
	for(my $j = $i + 1; $j < @tag_dirs; $j++) {
		@split2 = split("/", $tag_dirs[$j]);
		$command = "wc -l " . $output_dir . "/idr_" . $split[-1] . "_vs_" . $split2[-1] . ".txt";
		$num = `$command`;	
		@split_num = split('\s', $num);
		$orig_replicates{$i . "_" . $j} = $split_num[0]; 
		if($min_num > $split_num[0]) {
			$min_num = $split_num[0];
			$min_comp = $split[-1] . " vs " . $split[-2];
		}
		if($max_num < $split_num[0]) {
			$max_num = $split_num[0];
			$max_comp = $split[-1] . " vs " . $split[-2];
		}
	}
} 
my $max_numPeaks_Rep = $max_num;
if($max_num > 2 * $min_num) { 
	print STDERR "We need to flag some replicates here!\n";
	print STDERR "The lowest number of peaks between 2 replicates: " . $min_num . "\n";
	print STDERR "THe highest number of peaks between 2 replicates: " . $max_num . "\n";
	print STDERR "There is a more than 2 fold difference between these numbers!\n";
	print STDERR "Experiments need to be flagged!\n";
	print STDERR "IDR comparison with lowest number of peaks:\n";
	print STDERR "\t$min_comp\n";
	print STDERR "IDR comparison with highest number of peaks:\n";
	print STDERR "\t$max_comp\n";
	print FLAG "Replicates fail consistency test between the different replicates!\n";
	print FLAG "Number of IDR peaks of different replicates are DIFFERENT!\n\n";
	print FLAG "We need to flag some replicates here!\n";
	print FLAG "The lowest number of peaks between 2 replicates: " . $min_num . "\n";
	print FLAG "THe highest number of peaks between 2 replicates: " . $max_num . "\n";
	print FLAG "There is a more than 2 fold difference between these numbers!\n";
	print FLAG "Experiments need to be flagged!\n";
	print FLAG "IDR comparison with lowest number of peaks:\n";
	print FLAG "\t$min_comp\n";
	print FLAG "IDR comparison with highest number of peaks:\n";
	print FLAG "\t$max_comp\n\n";
} else {
	print STDERR "Number of IDR peaks of different replicates are similar!\n";
	print FLAG "Replicates pass consistency test between the different replicates!\n";
	print FLAG "Number of IDR peaks of different replicates are similar!\n\n";
}


#Look at self consistency thresholds
my %self_consistency;
my $self_min = 10000000;
my $self_max = 0;
my $self_min_comp = "";
my $self_max_comp = "";

for(my $i = 0; $i < @tag_dirs; $i++) {
	@split = split("/", $tag_dirs[$i]);
	$command = "wc -l " . $output_dir . "/idr_" . $split[-1] . ".txt";
	$num = `$command`;
	print $num;
	@split = split('\s', $num);
	$self_consistency{$tag_dirs[$i]} = $split[0];
	if($self_min > $split[0]) {
		$self_min = $split[0];
		$self_min_comp = $split[-1];
	}
	if($self_max < $split[0]) {
		$self_max = $split[0];
		$self_max_comp = $split[-1];
	}
}


if($self_max > 2 * $self_min) {
	print STDERR "Replicates fail self consistency test\n";
	print STDERR "We need to flag some replicates - self consistency test failed\n";
	print STDERR "The lowest number of peaks: " . $self_min. "\n";
	print STDERR "The highest number of peaks: " . $self_max . "\n";
	print STDERR "There is a more than 2 fold difference between these numbers\n";
	print STDERR "Experiments need to be flagged!\n";
	print STDERR "Experiment with lowest number of peaks:\n";
	print STDERR "\t$self_min_comp\n";
	print STDERR "Experiment with highest number of peaks:\n";
	print STDERR "\t$self_max_comp\n";
	print FLAG "Replicates fail self consistency test\n";
	print FLAG "Number of peaks between randomly split tag directories are DIFFERENT!\n";
	print FLAG "We need to flag some replicates - self consistency test failed\n";
	print FLAG "The lowest number of peaks: " . $self_min. "\n";
	print FLAG "The highest number of peaks: " . $self_max . "\n";
	print FLAG "There is a more than 2 fold difference between these numbers\n";
	print FLAG "Experiments need to be flagged!\n";
	print FLAG "Experiment with lowest number of peaks:\n";
	print FLAG "\t$self_min_comp\n";
	print FLAG "Experiment with highest number of peaks:\n";
	print FLAG "\t$self_max_comp\n\n";
} else {
	print STDERR "Replicates pass self consistency test\n";
	print FLAG "Replicates pass self consistency test\n";
	print FLAG "Number of peaks between randomly split tag directories are similar!\n\n";
}

#Look at pooled replicates
$command = "wc -l " . $output_dir . "/idr_pooled_pseudo_replicates.txt";
$num = `$command`;
print $num;
@split = split('\s', $num);
my $numPeaks_Rep0 = $split[0];

#Check if max_numPeaks_Rep and numPeaks_Rep0 are similar
if($numPeaks_Rep0/$max_numPeaks_Rep > 2 || $numPeaks_Rep0/$max_numPeaks_Rep < 0.5) {
	print STDERR "Replciates fail pooled consistency test\n";
	print FLAG "Replciates fail pooled consistency test\n";
	print FLAG "Maximal number of peaks between 2 replicates and number of peaks of pooled replicates are DIFFERENT\n";
	print STDERR "max_numPeaks_Rep and numPeaks_Rep are not similar!\n";
	print STDERR "max_numPeaks_Rep: " . $max_numPeaks_Rep . "\n";
	print STDERR "numPeaks_Rep: " . $numPeaks_Rep0 . "\n";
	print STDERR "numPeaks_Rep are the number of peaks from IDR analysis of the pooled pseudo replicates\n";
	print STDERR "max_numPeaks_Rep are the number of peaks from IDR analysis between the replicates\n";
	print FLAG "max_numPeaks_Rep and numPeaks_Rep are not similar!\n";
	print FLAG "max_numPeaks_Rep: " . $max_numPeaks_Rep . "\n";
	print FLAG "numPeaks_Rep: " . $numPeaks_Rep0 . "\n";
	print FLAG "numPeaks_Rep are the number of peaks from IDR analysis of the pooled pseudo replicates\n";
	print FLAG "max_numPeaks_Rep are the number of peaks from IDR analysis between the replicates\n\n";
} else {
	print STDERR "Replciates pass pooled consistency test\n";
	print FLAG "Replciates pass pooled consistency test\n";
	print FLAG "Maximal number of peaks between 2 replicates and number of peaks of pooled replicates pass test\n";
}

print STDERR "Generate output files\n";
print STDERR "Generate conservative and optimal peak file\n";	
open CONS, ">$output_dir/conservative_peakfile.txt";
open OPT, ">$output_dir/optimal_peakfile.txt";
open FH, "<$output_dir/peaks_pooled_tag_dirs.txt";

my $opt = 0;
my $cons = 0;
my $opt_threshold;
if($max_numPeaks_Rep > $numPeaks_Rep0) {
	$opt_threshold = $max_numPeaks_Rep;
} else {
	$opt_threshold = $numPeaks_Rep0;
}

foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "#") {
		next;
	}
	if($cons < $max_numPeaks_Rep) {
		print CONS $line . "\n";
		$cons++;
	}
	if($opt < $opt_threshold) {
		print OPT $line . "\n";
		$opt++;
	}	
}
close FH;
close CONS;
close OUT;
print STDERR "Everything is done!\n";
close FLAG;




sub thread_routine{
	my $command_in_thread = $_[0];
	my $id = $_[1];
	print STDERR "We are in the thread $id\t";
	print STDERR $command_in_thread . "\n";
	`$command_in_thread`;
}

sub generate_pseudo_replicates{
	my $tag_dir = $_[0];
	my @files = `ls $tag_dir/chr*`;
	@split = split("/", $tag_dir);
	my $o1 = $output_dir . "/" . $split[-1] . "_pseudo_1";
	my $o2 = $output_dir . "/" . $split[-1] . "_pseudo_2";
	if(-e $o1 || -e $o2) {
		print STDERR "Pseudo replicate for " . $tag_dir . " exists\n";
		exit;
	}
	`mkdir $o1`;
	`mkdir $o2`;
	$_ = 0 for my ($total_1, $total_2, $unique_1, $unique_2, $all_total_1, $all_total_2, $all_unique_1, $all_unique_2, $number_1, $number_2, $current, $rand);
	$_ = () for my(@name, %tags);
	#Run through all chromosome files
	foreach my $file (@files) {
		chomp $file;
		@name = split('/', $file);
		open OUT1, ">$o1/$name[-1]";
		open OUT2, ">$o2/$name[-1]";
		open FH, "<$file";
		$_ = 0 for($total_1, $total_2, $unique_1, $unique_2);
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$current = $split[4];
			$number_1 = 0;
			$number_2 = 0;
			while($current > 0) {
				$rand = rand(2);
				if($rand < 1) {
					$number_1++;
					$total_1++;
				} else {
					$number_2++;
					$total_2++;
				}
				$current--;
			}
			if($number_1 > 0) {
				print OUT1 $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . $split[3] . "\t" . $number_1 . "\t" . $split[5] . "\n";
				$unique_1++;
			}
			if($number_2 > 0) {
				print OUT2 $split[0] . "\t" . $split[1] . "\t" . $split[2] . "\t" . $split[3] . "\t" . $number_2 . "\t" . $split[5] . "\n";
				$unique_2++;
			}
		}
		$tags{$name[-1]}{'unique'}{1} = $unique_1;
		$tags{$name[-1]}{'unique'}{2} = $unique_2;
		$tags{$name[-1]}{'total'}{1} = $total_1;
		$tags{$name[-1]}{'total'}{2} = $total_2;
		$all_total_1 += $total_1;
		$all_total_2 += $total_2;
		$all_unique_1 += $unique_1;
		$all_unique_2 += $unique_2;
		close FH;
		close OUT1;
		close OUT2;
	}
	open FH, "<$tag_dir/tagInfo.txt";

	open OUT1, ">$o1/tagInfo.txt";
	open OUT2, ">$o2/tagInfo.txt";
	my $count = 0;
	foreach my $line (<FH>) {
		if($count != 1 && $count < 7) {
			print OUT1 $line;
			print OUT2 $line;
		}
		if($count == 1) {
			print OUT1 "genome=mm10\t" . $all_unique_1 . "\t" . $all_total_1 . "\n";	
			print OUT2 "genome=mm10\t" . $all_unique_2 . "\t" . $all_total_2 . "\n";	
		}
		if($count == 7) {
			last;
		}
		$count++;
	}
	close FH;
	foreach my $chr (sort {$a cmp $b} keys %tags) {
		print OUT1 $chr . "\t" . $tags{$chr}{'unique'}{'1'} . "\t" . $tags{$chr}{'total'}{'1'} . "\n";
		print OUT2 $chr . "\t" . $tags{$chr}{'unique'}{'2'} . "\t" . $tags{$chr}{'total'}{'2'} . "\n";
	}
	print OUT1 "cmd=split\n";
	print OUT2 "cmd=split\n";
	close OUT1;
	close OUT2;
}

sub convert_to_narrowPeaks {
	my $file = $_[0];
	open FH, "<$file";
	print STDERR $file . "\n";
	my $o = substr($file, 0, length($file) - 4) . ".narrowPeak.txt";
	my @split_line;
	open OUT, ">$o";
	foreach my $line (<FH>) {
		chomp $line;
		if(substr($line, 0, 1) eq "#") { next; }
		@split_line = split('\t', $line);
		print OUT $split_line[1] . "\t" . $split_line[2] . "\t" . $split_line[3] . "\t" . $split_line[0] . "\t0\t" . $split_line[4] . "\t" . $split_line[7] . "\t-1\t-1\t120\n";
	}
	close FH;
	close OUT;
}
