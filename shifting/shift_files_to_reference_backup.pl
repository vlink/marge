#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use threads;
require '../general/config.pm';
require '../general/system_interaction.pm';

$_ = 0 for my ($multithreading, $user_cores, $first, $end, $start, $current_pos, $sam, $peak, $format, $help);
$_ = "" for my ($dir, $organism, $name);
my (@files, @strains, $dbh, $sth, @split, $c, $shift, %shift_vector, $allele);
my $print = 1;
$allele = 1;

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "General commands:\n";
        print STDERR "\t-organism <organism>\n";
	print STDERR "\t-dir <directory with files to shift>: has to be the same strain\n";
	print STDERR "\t-files <list with files>: comma separated list\n";
	print STDERR "\t-strains: one or several strains - comma separated\n";
	print STDERR "\t-allele: allele that was used: default: 1\n";
	print STDERR "\tIf only one strain is specified all files are shifted from this strain\n";
	print STDERR "\tIf several strains are specified it has to match the number of files specified!\n";
	print STDERR "\n\nMultithreading:\n";
	print STDERR "\t-no-multithreading: No multithreading is used\n";
	print STDERR "\t-core <num>: Number of cores to use (default = 1/4 of all cores)\n";
	print STDERR "\n\nFormat to shift (default = sam):\n";
	print STDERR "\t-sam: shifts sam file (result file from mapping with bowtie or STAR\n";
	print STDERR "\t-homer: shifts HOMER peak and RNA files\n";
	print STDERR "\t-format: checks formats.txt to find the user defined format called shift_format\n";
	print STDERR "\n";
	print STDERR "-h | --help: prints help\n";
        exit;
}


if(@ARGV < 1) {
	&printCMD();
}

#TODO add resume then clean up code and then done
GetOptions(     "organism=s" => \$organism,
                "dir=s" => \$dir,
                "files=s{,}" => \@files,
                "strains=s{,}" => \@strains,
		"allele=s" => \$allele,
		"no-multithreading" => \$multithreading,
		"core=s" => \$user_cores,
		"sam" => \$sam,
		"homer" => \$peak,
		"format" => \$format,
		"h" => \$help,
		"help" => \$help)
	or die (&printCMD());

if($help == 1) {
	&printCMD();
}

if($sam == 0 && $peak == 0 && $format == 0) {
	$sam = 1;
}

if(system_interaction::check_module_installed("DBIx::Threaded") == 1) { $multithreading = 1; }
my $core = system_interaction::check_multithreading_support($multithreading, $user_cores);
my $total_chr = config::chromosome_number($organism);

my $split_sign;
my $shift_order;
my $chr_pos;

open FH, "<../config/formats.txt";

#Read format
foreach my $line (<FH>) {
	chomp $line;
	if($line eq "") { next; }
	@split = split(':', $line);
	if($sam == 1 && $split[0] eq "sam_format") {
		($split_sign, $chr_pos, $shift_order) = config::read_format($split[1]);
		last;
	}
	if($peak == 1 && $split[0] eq "homer_format") {
		($split_sign, $chr_pos, $shift_order) = system_interaction::read_format($split[1]);
		last;

	}
	if($format == 1 && $split[0] eq "user_format") {
		($split_sign, $chr_pos, $shift_order) = system_interaction::read_format($split[1]);
		last;
	}
}
if($chr_pos eq "") {
	print STDERR "No chromosome specified in your format!\n";
	exit;
}
if(@files == 0 && $dir eq "") {
	print STDERR "Files for shifting are missing!\n";
	exit;
}

#Read in files from directory
if(@files == 0) {
	@files = `ls $dir`;
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

#Fill strains array with same strain if there are more than one file
if(@strains == 1) {
	@strains = ($strains[0])x@files;
}

#Remove commas
for(my $i = 0; $i < @files; $i++) {
	chomp $files[$i];
		$files[$i] =~ s/,//g;
		$strains[$i] =~ s/,//g;
}

my $con = config::database_connection();
if($core == 1) {
	for(my $i = 0; $i < @files; $i++) {
		&shift_file($dir . "/" . $files[$i], $strains[$i]);
	}	
} else {
	
        my @running = ();
        my @Threads;
        my $current_thread;
	my $run = 0;
        while (scalar @Threads < (keys @files)) {
                @running = threads->list(threads::running);
                if (scalar @running < $core) {
                        my $thread = threads->new( sub { shift_file($files[$run], $strains[$run]) });
			$run++;
                        push (@Threads, $thread);
                        my $tid = $thread->tid;
                }
                @running = threads->list(threads::running);
                foreach my $thr (@Threads) {
                        if ($thr->is_running()) {
                                my $tid = $thr->tid;
                        }
                        elsif ($thr->is_joinable()) {
                                my $tid = $thr->tid;
                                $thr->join;
                        }
                }

                @running = threads->list(threads::running);
        }
        while (scalar @running > 0) {
                foreach my $thr (@Threads) {
                        $thr->join if ($thr->is_joinable());
                }
                @running = threads->list(threads::running);
        }
	threads->exit;
}


sub shift_file{
	use DBI;
        $dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

	my $f = $_[0];
	my $s = $_[1];
	$name = $f . "_shifted_from_" . $s;
	open OUT, ">$name";
	print STDERR "Shifting $f from $s to reference!\n"; 
	foreach my $c (keys %{$total_chr}) {
		print STDERR "chromosome: " . $c . "\n";
		%shift_vector = ();
		$first = 0;
		$sth = $dbh->prepare("SELECT * FROM offset_" . $organism . "_" . $s . "_allele_" . $allele . " WHERE chr = \'$c\' ORDER by strain_pos");
		$sth->execute();

		while(my $ref = $sth->fetchrow_hashref()) {
			if($first == 0) {
				$start = $ref->{'strain_pos'};
				$current_pos = $ref->{'strain_pos'};
				$shift_vector{$current_pos} = $ref->{'ref_pos'};
				$shift = $ref->{'shift'};
				$first++;
			} else {
				while($current_pos < $ref->{'strain_pos'}) {
					$shift_vector{$current_pos} = $current_pos - $shift;
					$current_pos++;
				}
				$shift_vector{$current_pos} = $ref->{'ref_pos'};
				$current_pos++;
				$end = $ref->{'strain_pos'};
				$shift = $ref->{'shift'};
			}
		}
		$sth->finish();
		$first = 0;
		open FH, "<$f";
		my $save = "";
		foreach my $line (<FH>) {
			chomp $line;
			$save = "";
			if(substr($line, 0, 1) eq "@" || substr($line, 0, 1) eq "#") {
				if($c eq 1) {
					print OUT $line . "\n";
				}
				next;
			}
			@split = split($split_sign, $line);
			if(substr($split[$chr_pos], 3) ne $c) {
				next;
			}
			for(my $j = 0; $j < @{$shift_order}; $j++) {
				if($shift_order->[$j] eq "shift" ) {
					if($split[$j] < $start) {
						$save .= $split[$j] . "\t";
					} elsif($split[$j] > $end) {
						$save .= "". ($split[$j] - $shift) . "\t";
					} else {
						$save .= $shift_vector{$split[$j]} . "\t";
					}
				} else {
					$save .= $split[$j] . "\t";
				}
			}
			for(my $j = @{$shift_order}; $j < @split; $j++) {
				$save .= $split[$j] . "\t";
			}
			chop $save;
			print OUT $save . "\n";
		}
		close FH;
	}
	close OUT;
	$dbh->disconnect();
}	
