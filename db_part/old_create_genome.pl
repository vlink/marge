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
	print STDERR "\t-organism <organism>\n";
	print STDERR "\t-strains <comma sepearted list of strains - no spaces!>\n";
	print STDERR "\tif this parameter is empty the genome is created for all possible strains!\n";
	print STDERR "\t-out <name of the output directory>: Script automatically creates for each strain a directory with this name + strain\n";
	print STDERR "\n\nMultithreading options\n";
	print STDERR "\t-no-multithreading: Turns multithreading off\n";
	print STDERR "\t-core <num>: Number of cores to use\n";
	print STDERR "\tIf nothing is specified this script uses multithreading with 1/4 of all cores if possible\n";
	print STDERR "\t-resume: Starts from last successful log point!\n";
	exit;	
}

my $organism = "";
my $strains = "";
my $output_dir = "";
my $multithreading = 0;
my $user_cores = 0;
my $resume = 0;

GetOptions(	"organism=s" => \$organism,
		"strains=s" => \$strains,
		"out=s" => \$output_dir,
		"no-multithreading" => \$multithreading,
		"core=s" => \$user_cores,
		"resume" => \$resume)
or die ("Error in command line arguments!\n");

if($organism eq "" || $output_dir eq "") {
	print STDERR "Missing parameters!\n";
	exit;
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
 
my %strains;
my $command;
my $table_name = $organism . "_mutations";
my @split;

#Check strains
if($strains eq "") {
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
	$command = "SELECT column_name FROM information_schema.columns WHERE table_name=\'$table_name\'";
	$sth = $dbh->prepare($command);
	$sth->execute();
	while(my $s = $sth->fetchrow_arrayref()) {
		if($s->[0] ne "chromosome" && $s->[0] ne "position" && $s->[0] ne "reference") {
			$strains{$s->[0]} = 1;
		}
	}
} else {
	@split = split(',', $strains);
	foreach my $s (@split) {
		$strains{$s} = 1;
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
			if(exists $strains{$line}) {
				delete $strains{$line};
			}
		}
		close IN;
		open LOG, ">>genome_create.log";
		if(keys %strains == 0) {
			print STDERR "All strains were added successfully!\n";
			print STDERR "Nothing to do!\n";
			exit;
		}	
	}
} else {
	open LOG, ">genome_create.log";
}

my $total_chr = config::chromosome_number($organism);
my $num;

my %hash;
if($multithreading == 1) {
	foreach my $s (keys %strains) {
		print $s . "\n";
		&generate_genome($s);
	}
} else {
        foreach my $s (keys %strains) {
                $table_name = "offset_" . $organism . "_" . $s;
                $sth = $dbh->prepare("SELECT COUNT(*) AS c from $table_name");
                $sth->execute();
                $num = $sth->fetchrow_hashref();
                $hash{$s} = $num->{'c'};
        }
        my @running = ();
        my @Threads;
        my $current_thread;
        while (scalar @Threads < (keys %strains)) {
#       while (scalar @Threads < 0) {
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
exit;

close LOG;

#$dbh = DBI->connect("DBI:Pg:dbname=vlink;host=localhost", "vlink", "mousestrands", {'RaiseError' => 1});
#$dbh->do("SET search_path TO inbred, public");

my $genome;
my $mutations;
my $mut_ref;
my $end;
my $seq;
my $c;
my $last_pos;
my $strain;
my $first = 0;
my $name;


sub generate_genome {
	$strain = $_[0];
	$first = 0;
	$last_pos = 0;

	$output_dir = $output_dir . "_" . $strain;
	if(!-e $output_dir) {
		`mkdir $output_dir`;
	}

	for(my $i = 1; $i <= $total_chr; $i++) {
		if($i == $total_chr) {
			$c = "Y";
		} elsif($i == ($total_chr - 1)) {
			$c = "X";
		} else {
			$c = $i;
		}
		$name = "chr" . $c . ".fa";
		open OUT, ">$output_dir/$name";
		print OUT ">chr$c\n";
		$last_pos = 0;
		#Fetch mutations
		$mutations = $dbh->prepare("SELECT position, reference, $strain FROM " . $organism . "_mutations where chromosome = \'$c\' AND reference != $strain ORDER BY position");
		$mutations->execute();
		$mut_ref = $mutations->fetchrow_hashref();

		#Fetch reference genome
		$genome = $dbh->prepare("SELECT * FROM " . $organism . "_ref_genome where chr = \'$c\' ORDER by position");
		$genome->execute();

		#Create genome
		while(my $genome_hash = $genome->fetchrow_hashref()) {
			#Create the N's at the beginning - we start with 1 which means 100 so we need to go till last_pos + 1
			while($genome_hash->{'position'} > $last_pos + 1) {
				for(my $j = 0; $j < 100; $j++) {
					print OUT "N";
				}
				print OUT "\n";
				$last_pos++;
			}
			$last_pos = $genome_hash->{'position'};
			#No mutation in this line, so just print this genome line
			if((exists $mut_ref->{'position'} && $mut_ref->{'position'} > ($genome_hash->{'position'} * 100 + 100)) || (!exists $mut_ref->{'position'})) {
				print OUT substr($genome_hash->{'seq'}, $first) . "\n";
				$first = 0;
				next;
			}
			$seq = "";
			#mutation is between beginning and end of this line
			while(exists $mut_ref->{'position'} && $mut_ref->{'position'} >= ($genome_hash->{'position'} * 100 +1) && $mut_ref->{'position'} <= ($genome_hash->{'position'} * 100 +100)) {
				#mutation at last position of the line
				if(substr($mut_ref->{'position'}, length($mut_ref->{'position'}) - 2, 2) eq "00") {
					$end = 100 - $first - 1;
				#mutation within the line
				} else {
					$end = substr($mut_ref->{'position'}, length($mut_ref->{'position'}) - 2, 2) - 1 - $first;
				}
				$seq = $seq . substr($genome_hash->{'seq'}, $first, $end);
				$seq = $seq . $mut_ref->{$strain};
				
				#Calculate new start position
				$first = substr($mut_ref->{'position'}, length($mut_ref->{'position'}) - 2, 2) + length($mut_ref->{'reference'}) - 1;
				
				#Mutation goes over line	
				if(substr($mut_ref->{'position'}, length($mut_ref->{'position'}) - 2, 2) eq "00") {
					$first = 99 + length($mut_ref->{'reference'});
				}
				$mut_ref = $mutations->fetchrow_hashref();
			}
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
		close OUT;
	}
	print LOG $strain . "\n";
#	$sth->finish();
	return 0;
}
