#!/usr/bin/perl -w

#Script to interact with the database that saves all the mutation data of the mouse strains

use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';
require '../general/system_interaction.pm';

#Create database connection
my $dbh;
my $sth;

if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "\t-genome <genome>\n";
	print STDERR "\t-strains <strains>: one or several strains (comma separated) - all or not specified: all strains\n";
	print STDERR "\tif this parameter is empty the genome is created for all possible strains!\n";
	print STDERR "\t-out <name of the output directory>: Script automatically creates for each strain a directory with this name + strain\n";
	print STDERR "\n\nMultithreading options\n";
	print STDERR "\t-no-multithreading: Turns multithreading off\n";
	print STDERR "\t-core <num>: Number of cores to use\n";
	print STDERR "\tIf nothing is specified this script uses multithreading with 1/4 of all cores if possible\n";
	print STDERR "\t-resume: Starts from last successful log point!\n";
	print STDERR "\t-homo: Genome is homozygous\n";
	exit;	
}

$_ = () for my($num, %hash, $genome, $mutations, $mut_ref, $end, $seq, $c, $last_pos, $strain, $name, $strains, @split);
$_ = "" for my($organism, $strains_input, $output_dir, $table_name, $dir);
$_ = 0 for my ($multithreading, $user_cores, $resume, $first, $homo);
my $allele = 2;

my %mandatory = ('-genome' => 1, '-out' => 1);
my %convert = map { $_ => 1 } @ARGV;

config::check_parameters(\%mandatory, \%convert);


GetOptions(	"genome=s" => \$organism,
		"strains=s" => \$strains_input,
		"out=s" => \$output_dir,
		"no-multithreading" => \$multithreading,
		"core=s" => \$user_cores,
		"resume" => \$resume, 
		"homo" => \$homo)
or die ("Error in command line arguments!\n");

if($homo == 1) {
	$allele = 1;
}

#Start database connectio
if(system_interaction::check_module_installed("DBIx::Threaded") == 1) { $multithreading = 1;}
my $con = config::database_connection();
my $core = system_interaction::check_multithreading_support($multithreading, $user_cores);
if($core == 1) {
	$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
} else {
	use threads;
	$dbh = DBIx::Threaded->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
}
 
#Check strains
if($strains_input eq "") {
	#now we need to do all strains
	print STDERR "No strains are specified!\n";
	if($resume == 0) {
		print STDERR "Are you sure you want to create the genomes for all strains?\n";
	} else {
		print STDERR "Are you sure you want to create the genomes for all strains that are not generated already?\n";
	}
	print STDERR "Hit Ctrl + C to abort!\n";
	print STDERR "Waiting for 10 seconds\n";
	for(my $i = 1; $i < 10; $i++) {
		print STDERR " . ";
		sleep(1);
	}
	print STDERR "\n";
	print STDERR "Start to generate genomes for all strains!\n";
	$strains = config::get_strains($organism);
}  elsif($strains_input eq "all") {
	$strains = config::get_strains($organism);
} else {
	@split = split(',', $strains_input);
	foreach my $s (@split) {
		$strains->{$s} = 1;
	}
}

if($resume == 1) {
	#Check logfile
	print STDERR "Resume with genomes that were not successfully generated\n";
	if(!-e "genome_create.log") {
		print STDERR "Log file is missing!\n";
		print STDERR "Proceed without resuming\n";
		open LOG, ">genome_create.log";
	} else {
		open IN, "<genome_create.log";
		foreach my $line (<IN>) {
			chomp $line;
			if(exists $strains->{$line}) {
				delete $strains->{$line};
			}
		}
		close IN;
		open LOG, ">>genome_create.log";
		if(keys %{$strains} == 0) {
			print STDERR "All strains were added successfully!\n";
			print STDERR "Nothing to do!\n";
		}	
	}
} else {
	open LOG, ">genome_create.log";
}

my $total_chr = config::chromosome_number($organism);

foreach my $s (keys %{$strains}) {
	$table_name = "offset_" . $organism . "_" . $s . "_allele_1";
	$sth = $dbh->prepare("SELECT COUNT(*) AS c from $table_name");
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();
	$hash{$s} = $num->{'c'};
}

$dbh->disconnect();

if($multithreading == 1) {
	foreach my $s (keys %{$strains}) {
		print $s . "\n";
		&generate_genome($s);
	}
} else {
	my @running = ();
	my @Threads;
	my $current_thread;
	while (scalar @Threads < (keys %{$strains})) {
		@running = threads->list(threads::running);
		if (scalar @running < $core) {
			#Find out which thread we need to start;
			foreach my $keys (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
				$current_thread = $keys;
			}
			delete $hash{$current_thread};
			my $thread = threads->new( sub { generate_genome($current_thread) });
			print STDERR "Fetch genome for " . $current_thread . "\n";
			sleep(4);
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
	while (scalar @running > 1) {
		foreach my $thr (@Threads) {
			$thr->join if ($thr->is_joinable());
		}
		@running = threads->list(threads::running);
	}
}

close LOG;

sub generate_genome {
	my $dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
	$strain = $_[0];
	
	for(my $a = 0; $a < $allele; $a++) {
		$dir = $output_dir . "_" . $strain . "_allele_" . ($a + 1);
		if(!-e $dir) {
			`mkdir $dir`;
		}

		foreach my $c (keys %{$total_chr}) {
			$first = 0;
			$end = 0;
			$name = "chr" . $c . ".fa";
			open OUT, ">$dir/$name";
			print OUT ">chr$c\n";
			$last_pos = -1;
			#Fetch mutations
			$mutations = $dbh->prepare("SELECT * FROM " . $organism . "_mutations_" . $strain . "_allele_" . ($a + 1) . " WHERE chr = \'$c\' ORDER BY pos");
			$mutations->execute();
			$mut_ref = $mutations->fetchrow_hashref();

			#Fetch reference genome
			$genome = $dbh->prepare("SELECT * FROM " . $organism . "_ref_genome where chr = \'$c\' ORDER by pos");
			$genome->execute();

			#Create genome
			while(my $genome_hash = $genome->fetchrow_hashref()) {
				#Create the N's at the beginning - we start with 1 which means 100 so we need to go till last_pos + 1
				while($genome_hash->{'pos'} > $last_pos + 1) {
					for(my $j = 0; $j < 100; $j++) {
						print OUT "N";
					}
					print OUT "\n";
					$last_pos++;
				}
			#	print STDERR "GENOME ($genome_hash->{'pos'} - chr: $c) \n" . $genome_hash->{'seq'} . "\n";
				$last_pos = $genome_hash->{'pos'};
				#No mutation in this line, so just print this genome line
				if((exists $mut_ref->{'pos'} && $mut_ref->{'pos'} > ($genome_hash->{'pos'} * 100 + 100) && $first < length($genome_hash->{'seq'})) || (!exists $mut_ref->{'pos'})) {
					print OUT substr($genome_hash->{'seq'}, $first) . "\n";
					$first = 0;
					next;
				}
				$seq = "";
				#mutation is between beginning and end of this line
				while(exists $mut_ref->{'pos'} && $mut_ref->{'pos'} >= ($genome_hash->{'pos'} * 100 +1) && $mut_ref->{'pos'} <= ($genome_hash->{'pos'} * 100 +100)) {
					#mutation at last position of the line
			#		print STDERR "first: " . $first . " end: " . $end . "\n";
			#		print STDERR "mut : " . $mut_ref->{'pos'} . "\t" . $mut_ref->{'reference'} . "\t" . $mut_ref->{'strain'} . "\n";
					if($first < 100) {
						if(substr($mut_ref->{'pos'}, length($mut_ref->{'pos'}) - 2, 2) eq "00") {
							$end = 100 - $first - 1;
						#mutation within the line
						} else {
							$end = substr($mut_ref->{'pos'}, length($mut_ref->{'pos'}) - 2, 2) - 1 - $first;
						}
			#			print STDERR "NEW: first: " . $first . " end: " . $end . "\n";
			#			print STDERR $genome_hash->{'seq'} . "\n";
			#			print STDERR $first . "\t" . $end . "\n";
						$seq = $seq . substr($genome_hash->{'seq'}, $first, $end);
						$seq = $seq . $mut_ref->{'strain'};
						#Calculate new start position
						$first = $end + length($mut_ref->{'reference'});
						$first = substr($mut_ref->{'pos'}, length($mut_ref->{'pos'}) - 2, 2) + length($mut_ref->{'reference'}) - 1;
			#			print STDERR "Seq so far:\n" . $seq . "\n";	
						#Mutation goes over line	
						if(substr($mut_ref->{'pos'}, length($mut_ref->{'pos'}) - 2, 2) eq "00") {
							$first = 99 + length($mut_ref->{'reference'});
						}
						$mut_ref = $mutations->fetchrow_hashref();
					} else {
						$first = $first - 100;
		#				print STDERR "NEXT\n";
					}
				}
		#		print STDERR "DONE\n\n";
				print OUT $seq ;
				#New line longer than 100 - Recalcualte start point for next line
				if($first > 100) {
					print OUT "\n";
					$first = $first - 100;
				} else {
					print OUT substr($genome_hash->{'seq'}, $first) . "\n";
					$first = 0;
				}
			}
			$mutations->finish();
			$genome->finish();
			close OUT;
		}
	}
	print LOG $strain . "\n";
	$dbh->disconnect();
	return 0;
}
