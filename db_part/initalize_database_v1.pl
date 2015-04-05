#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';
require '../general/system_interaction.pm';
use threads;

my $organism = "";
my $num;
my $snp = "";
my $indel = "";
my $both = "";
my $filter = 0;
my $same = 0;
my @lines;
my $genome_folder = "";
my $genome_file = "";
my $core_user = 0;
my $multithreading = 0;
my $annotation  = "";
my $genes = "";
my $name = "";
my $resume = 0;
if(@ARGV < 1) {
	print STDERR "Usage:\n";
	print STDERR "General commands:\n";
	print STDERR "\t-organism <organism>\n";
	print STDERR "\t-reference_folder <path to folder with reference genome in it split by chromosomes>\n";
	print STDERR "\t-reference_file <path to file with reference genome - all chromosomes in one file\n";
	print STDERR "\n\nInput files:\n";
	print STDERR "\t-snp <file> - file containing only SNPs\n";
	print STDERR "\t-indel <file>  - file containing only indels\n";
	print STDERR "\t-both <file> - file containing SNPs and Indels\n";
	print STDERR "\n\n";
	print STDERR "Filter vcf file.\n";
	print STDERR "\t-filter - Filters out all mutations that did not pass all filters\n"; 	
	print STDERR "\t-same - Filters out all mutations that are heterozyous\n";
	print STDERR "\tif position is heterozygous and -same is not defined the first allele is taken without any further evaluation\n";
	print STDERR "\n\nMultithreading:\n";
	print STDERR "\t-no-multithreading: does not use multithreading to create vectors for database\n";
	print STDERR "\t-core <num>: Number of cores to use for multithreading - if not specified and multithreading is not turned off and possible it uses 1/4 of all cores\n";
	print STDERR "\n\nAnnotation for genome:\n";
	print STDERR "\t-ann <file>: General annotations for parts or the genome downloaded from ucsc\n";
	print STDERR "TODO: Find out what the name is exactly!!!!!\n";
	print STDERR "\t-gene <gene file>: File with annotations for introns and exons\n";
	print STDERR "\t-names <file with gene names>: File that maps transcript names to gene names\n";
	print STDERR "\tIf forrmat different from standard format that is used, please change config file formats.txt\n";
	print STDERR "TODO: Find out what the name of the files is exactly!\n";
	print STDERR "TODO: ADD FORMAT SUPPORT!\n";
	print STDERR "\t-resume: Restart from the last log point\n";
	exit;
}


#TODO add resume then clean up code and then done
GetOptions(	"organism=s" => \$organism,
		"snp=s" => \$snp,
		"indel=s" => \$indel,
		"both=s" => \$both,
		"filter" => \$filter,
		"same" => \$same,
		"reference_folder=s" => \$genome_folder,
		"reference_file=s" => \$genome_file,
		"core=s" => \$core_user,
		"no-multithreading" => \$multithreading,
		"ann=s" => \$annotation,
		"gene=s" => \$genes,
		"names=s" => \$name,
		"resume" => \$resume)
or die("Error in commandline arguments!\n");

#Initialize logpoints
my $check_muts = 0;
my $check_ref = 0;
my %check_offset = ();
my $check_annotation = 0;
my $check_genes = 0;
my %to_delete=();
my @split;

if($resume == 1) {
	print STDERR "Start from last successful logpoint\n";
	if(!-e "database_init.log") {
		print STDERR "Logfile is missing!\n";
		print STDERR "Hit Ctrl + C to interrupt\n";
		print STDERR "Start to initalize the database in 15 seconds!\n";
		for(my $i = 1; $i < 15; $i++) {
			print STDERR " . ";
			sleep(1);
		}
		print STDERR "\n";
		$resume = 0;
	} else {
		open IN, "<database_init.log";
		foreach my $line (<IN>) {
			chomp $line;
			if(substr($line, 0, 8) eq "organism") {
				if(substr($line, 9) ne $organism) {
					print STDERR "Log file is for another organsim!\n";
					print STDERR "Need to restart!\n";
					print STDERR "Hit Ctrl + C to interrupt\n";
					print STDERR "Start to initalize the database in 15 seconds!\n";
					for(my $i = 1; $i < 15; $i++) {
						print STDERR " . ";
						sleep(1);
					}
					print STDERR "\n";
					last;
				}
			}
			if($line eq $organism . "_mutations") {
				$check_muts = 1;
				$to_delete{"tmp_merged.txt"} = 1;
			}	
			if($line eq $organism . "_ref_genome") {
				$check_ref = 1;
				$to_delete{"tmp_ref"} = 1;
			}
			if(substr($line, 0, 6) eq "offset") {
				@split = split('_', $line);
				$check_offset{$split[-1]} = 1;
				$to_delete{"shift_vector_" . $split[-1]} = 1;
			}
			if($line eq $organism . "_annotations") {
				$check_annotation = 1;
				$to_delete{"tmp_ann"} = 1;
			}
			if($line eq $organism . "_genes") {
				$check_genes = 1;
				$to_delete{"tmp_genes"} = 1;
			}
		}
		close IN;
		open LOG, ">>database_init.log";
	}
} else {
	open LOG, ">database_init.log";
}
print LOG "organism:" . $organism  . "\n";
if($genome_folder eq "" && $genome_file eq "") {
	print STDERR "The reference genome is missing!\n";
	exit;
}
#Check multithreading - we need db module for that
if(system_interaction::check_module_installed("DBIx::Threaded") == 1) { $multithreading = 1; }
#Check options

if(($snp ne "" || $indel ne "") && $both ne "") {
	print STDERR "You can not use this combination of files!\n";
	print STDERR "Please specify either two files, one containing SNPs, the other containing Indels\nor one file contaning both!\n";
	exit;
}

my $dbh;
my $sth;

my $header;

my $con = config::database_connection();
if($multithreading == 1) {
	$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
} else {
	$dbh = DBIx::Threaded->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
}

#Merge SNP and indel file if both exists
my $change = 0;
my $snps;
my $indels;
my @s;
my @ref;
my @i;
my $table_name;
my $f1;
my $f2;
my $command;
my @h;
my $current_line = "";
my $last_line = "";
if($check_muts == 0) {
	#open files and filter out header
	if($both ne "") {
		$snp = $both;
	}
	if($snp ne "") {
		open($f1, $snp);
		$snps = read_file_line($f1);
		@s = split('\t', $snps);
		while(substr($s[0], 0, 1) eq "#") {
			$header = $snps;
			$snps = read_file_line($f1);
			@s = split('\t', $snps);
		}

	}
	if($indel ne "") {
		open($f2, $indel);
		$indels = read_file_line($f2);
		@i = split('\t', $indels);
		while(substr($i[0], 0, 1) eq "#") {
			$header = $indels;
			$indels = read_file_line($f2);
			@i = split('\t', $indels);
		}

	}
	#goto DATABASE;
	print STDERR "Merge files and mutations!\n";
	print STDERR "\tchromosome 1\n";


	my $check_snp;
	my $check_indel;
	my $s_c;
	my $i_c;
	#SNP and Indel file are separate - merge the files and also the mutations at the same time
	while($snps and $indels) {
		@s = split('\t', $snps);
		$s_c = &convert_chr_to_num($s[0]);
		$check_snp = &print(\@s);
		if($check_snp eq "") {
			$snps = read_file_line($f1);
			next;
		}
		@i = split('\t', $indels);
		$i_c = &convert_chr_to_num($i[0]);
		$check_indel = &print(\@i);
		if($check_indel eq "") {
			$indels = read_file_line($f2);
			next;
		}
		$last_line = $current_line;
		if($s_c > $i_c) {
			if($change == 0) {
				print STDERR "\tchromosome " . $s[0] . "\n";
			}
			$current_line = $check_indel;
	#		print $indels . "\n";
			&merge($current_line, $last_line);
			$indels = read_file_line($f2);		
			$change = 1;
		}
		if($i_c > $s_c) {
			if($change == 0) {
				print STDERR "\tchromosome " . $i[0] . "\n";
			}
			$current_line = $check_snp;
	#		print $snps . "\n";
			&merge($current_line, $last_line);
			$snps = read_file_line($f1);
			$change = 1;
		}

		if($s_c == $i_c) {
			$change = 0;
			if($s[1] < $i[1]) {
				$current_line = $check_snp;
	#			print $snps . "\n";
				&merge($current_line, $last_line);
				$snps = read_file_line($f1);
			} else {
				$current_line = $check_indel;
	#			print $indels . "\n";
				&merge($current_line, $last_line);
				$indels = read_file_line($f2);
			}
		}
	}
	while($snps) {
		@s = split('\t', $snps);
		$check_snp = &print(\@s);
		$current_line = $check_snp;
		&merge($current_line);
		$snps = read_file_line($f1);
	}
	while($indels) {
		@i = split('\t', $indels);
		$check_indel = &print(\@i);
		$current_line = $check_indel;
		&merge($current_line);
		$indels = read_file_line($f2);
	}
	open OUT, ">tmp_merged.txt";

	for(my $i = 0; $i < @lines; $i++) {
		print OUT $lines[$i] . "\n";
	}
	close OUT;
	print STDERR "ADD TO DATABASE\n";
	$table_name = $organism . "_mutations";
	#Check if table already exists
	if(&check_if_table_exists($table_name) == 1) {
		print STDERR "Table " . $table_name . " already exists!\n";
		print STDERR "Aborting....\n";
		$dbh->do("DROP TABLE " . $table_name);
	#	exit;
	}
	#Create table
	$command = "CREATE TABLE " . $table_name . "(chromosome varchar(2) not null, position bigint not null, reference varchar(116) not null,";
	@h = split('\t', $header);
	my @a;
	print STDERR "Creating table with these strains\n";
	for(my $i = 9; $i < @h; $i++) {
		print STDERR "\t" . lc($h[$i]) . "\n";
		if($h[$i] =~ /^\d/) {
			print STDERR "\tString starts with a number - not a possible naming in postgresql db\n";
			print STDERR "\tRenaming to: ";
			$h[$i] = "s" . lc($h[$i]);
			print STDERR lc($h[$i]) . "\n";
		}
		$command .= lc($h[$i]) . " varchar(116) not null,";
		$h[$i] = lc($h[$i]);
	}
	$command = $command . "PRIMARY KEY(position, chromosome))";
	$dbh->do($command);
	print STDERR "\n\nInserting mutations into database\n";
	print STDERR "This may take a while\n";
	`bash database_insert.sh $table_name tmp_merged.txt`;
	print STDERR "Done!\n\n\n";

	print STDERR "Check if database input was successful\n";
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	if(@lines == $num->{'c'}) {
		print LOG $table_name . "\n";
		$to_delete{'tmp_merged.txt'} = 1;
	} else {
		print STDERR "Input in database was not successful!\n";
		exit;
	}
	print STDERR "Successful\n\n";
} else {
	my $h = `grep #CHROM $snp`;
	chomp $h;
	@h = split('\t', $h);
	print STDERR "Mutations were already successfully added to database\n\n";
	for(my $i = 9; $i < @h; $i++) {
		if($h[$i] =~ /^\d/) {
			$h[$i] = "s" . lc($h[$i]);
		}
		$h[$i] = lc($h[$i]);
	}

}

if($check_ref == 0) {
	print STDERR "Add reference genome!\n";
	my @files;
	if($genome_folder ne "") {
		@files = `ls -1 $genome_folder/*`;
		&ref_genome(@files);
	} 
	if($genome_file ne "") {
		$files[0] = $genome_file;
		&ref_genome(@files);
	}

	open OUT, ">tmp_ref";
	for(my $i = 0; $i < @ref; $i++) {
		print OUT $ref[$i] . "\n";
	}
	close OUT;
	$table_name = $organism . "_ref_genome";
	if(&check_if_table_exists($table_name) == 1) {
		print STDERR "Table " . $table_name . " already exists!\n";
		print STDERR "Aborting...\n";
		$dbh->do("DROP TABLE " . $table_name);
	#	exit;
	}

	#Create table
	$command = "CREATE TABLE " . $table_name . " (position bigint not null, seq varchar(100) not null, chr varchar(20) not null, PRIMARY KEY(position, chr))";
	$dbh->do($command);

	print STDERR "\n\nInserting reference genome into database\n";
	print STDERR "This may take a while\n\n";
	`bash database_insert.sh $table_name tmp_ref`;
	print STDERR "Check if database insert was successful!\n";
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	if(@ref == $num->{'c'}) {
		print LOG $table_name . "\n";
		$to_delete{"tmp_ref"} = 1;
	} else {
		print STDERR "Input in database not successful!\n";
		exit;
	}
	print STDERR "Successful!\n\n";
}

my @new;

use Data::Dumper;

my %strains_to_use = ();

if(keys %check_offset > 0) {
	for(my $i = 9; $i < @h; $i++) {
		if(!exists $check_offset{$h[$i]}) {
			$strains_to_use{$i} = 1;
		}
	}
	print STDERR "" . (@h - 9 - (keys %strains_to_use)) . " strains were already successfully inserted\n";
	print STDERR "Insert " . (keys %strains_to_use) . " strains to database\n\n";
} else {
	for(my $i = 9; $i < @h; $i++) {
		$strains_to_use{$i} = 1;
	}
}
my $core = system_interaction::check_multithreading_support($multithreading, $core_user);
my %hash_copy;
if(keys %strains_to_use > 0) {
	print STDERR "Create shifting vectors!\n";
}
if($core == 1) {
	foreach my $key (keys %strains_to_use) {
		print STDERR "vector for " . $h[$key] . "\n";
		&shift_vector($h[$key]);
	}
} else {
	#First sort strains according to mutations so the threads finish in the right order
	my %hash;
	foreach my $i (keys %strains_to_use) {
		$table_name = $organism . "_mutations"; 
		$sth = $dbh->prepare("SELECT COUNT(*) AS c from $table_name WHERE length(reference) != length($h[$i])");
		$sth->execute();
	       	$num = $sth->fetchrow_hashref();
        	$hash{$i} = $num->{'c'};
		$hash_copy{$i} = $num->{'c'};
	}
	my @running = ();
	my @Threads;
	my $current_thread;
	while (scalar @Threads < (keys %strains_to_use)) {
#	while (scalar @Threads < 0) {
		@running = threads->list(threads::running);
		if (scalar @running < $core) {
			#Find out which thread we need to start;
			foreach my $keys (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
				$current_thread = $keys;
			}
			delete $hash{$current_thread};
			my $thread = threads->new( sub { shift_vector($h[$current_thread]) });
			$to_delete{"shift_vector_" . $h[$current_thread]} = 1;
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
#	threads->exit();
}

if(keys %strains_to_use > 0) {
	print STDERR "\nCheck if database insert was successful!\n";
	foreach my $i (keys %strains_to_use) {
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM offset_" . $organism . "_$h[$i]");
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		if($hash_copy{$i} == $num->{'c'}) {
			print LOG "offset_" . $organism . "_$h[$i]\n";
		} else {
			print STDERR "Not successful!\n";
			exit;
		}
	}
}
print STDERR "Successful!\n\n";

my $file_lines;
if($annotation ne "") {
	if($check_annotation == 0) {
		print STDERR "Loading annotations into database!\n";
		open FH, "<$annotation";
		open OUT, ">tmp_ann";
		my $string;
		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			$string = substr($split[1], 3) . "," . $split[2] . "," . $split[3] . "," . $split[4] . "," . $split[5] . ",";
			if($split[5] eq "TTS") {
				$string .= substr($split[0], 5, length($split[0]) - 6);
			}
			if($split[5] eq "E" || $split[5] eq "I" || $split[5] eq "3UTR" || $split[5] eq "5UTR") {
				@s = split(",", $split[0]);
				@split = split('\(', $s[0]);
				$string .= $split[1];
			}
			print OUT $string . "\n";
		}
		close FH;
		close OUT;
		$table_name = $organism . "_annotations";
		if(&check_if_table_exists($table_name) == 1) {
			print STDERR "Table " . $table_name . " exists!\n";
			print STDERR "Aborting...";
			$dbh->do("DROP TABLE $table_name");
		}
		$command = "CREATE TABLE " . $table_name . "(chr varchar(20) not null, start bigint not null, stop bigint not null, strand varchar(1) not null, annotation varchar(64) not null, gene varchar(116) not null)";
		$dbh->do($command);
		`bash database_insert.sh $table_name tmp_ann`;
		print STDERR "\nCheck if database insert was successful!\n";
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$file_lines = `wc -l tmp_ann`;
		@s = split('\s+', $file_lines);
		if($num->{'c'} == $s[0]) {
			print LOG "$table_name\n"; 
			$to_delete{"tmp_ann"} = 1
		} else {
			print STDERR "Not successful!\n";
			exit;
		}
		print STDERR "Annotations are loaded into database\n\n";
	} else {
		print STDERR "Annotations were already successfully loaded into database\n\n";
	}
}

if($genes ne "") {
	if($check_genes == 0) {
		my %save_names;
		if($name ne "") {
			print STDERR "Names are specified!\n";
			open FH, "<$name";
			my($name, $nf) = config::name_format($organism);
			foreach my $line (<FH>) {
				chomp $line;
				@split = split('\t', $line);
				if(@split < $name) {
					next;
				}
				foreach my $k (keys %{$nf}) {
					$save_names{$split[$k]} = $split[$name];
				}	
			}
			close FH;
			print STDERR "Names are saved!\n";
		} else {
			print STDERR "No names for transcripts are specified\n";
		}
		my $gf = config::gene_format($organism);
		my @introns;
		my @parts;

		my $start;
		my $stop;
		my $exon;
		my $string = "";
		my $blub = 0;
		open OUT, ">tmp_genes";
		if(exists $gf->{'intron_exon'}) {
			print STDERR "Intron and exon are listed together - calculate it!\n";
			open FH, "<$genes";
			foreach my $line (<FH>) {
				chomp $line;
				$start = 0;
				$stop = 0;
				$exon = 1;
				@split = split('\t', $line);
				@introns = split(',', $split[-1]);
				$string = "";
				#Check the strand
				if($split[4] eq "-") {
					for(my $i = @introns - 1; $i > 0; $i--) {
						@parts = split(':', $introns[$i]);
						if(length($parts[0]) > 4) {
							next;
						}
						if(substr($parts[0], 0, 1) eq "E") {
							$start = $parts[1] - 1;
							$i--;
							@parts = split(':', $introns[$i]);
							if(length($parts[0]) < 6 && substr($parts[0], 0, 1) ne "I") {
								print STDERR "Weird!\n";
								exit;
							} 
							$stop = $parts[1] - 1;
							if($name ne "" && exists $save_names{$split[0]}) {
								print OUT $save_names{$split[0]} . ",";
							} else {
								print OUT "Unknown,";
							}
							print OUT $split[0] . "," . $exon . "," . substr($split[1], 3) . "," . $start . "," . $stop . "," . $split[4] . "\n";	
							$exon++;
						}
					}
				} else {
					for(my $i = 0; $i < @introns; $i++) {
						@parts = split(':', $introns[$i]);
						if(length($parts[0]) > 4) {
							next;
						}
						if(substr($parts[0], 0, 1) eq "E") {
							$start = $parts[1] - 1;
							$i++;
							if(@introns > $i) {
								@parts = split(':', $introns[$i]);
								if(length($parts[0]) < 4 && substr($parts[0], 0, 1) ne "I") {
									print STDERR "Weird pos strand\n";
									exit;
								}
								$stop = $parts[1] - 1;
							} else {
								$stop = $split[3];
							}
							if($name ne "" && exists $save_names{$split[0]}) {
								print OUT $save_names{$split[0]} . ",";
							} else {
								print OUT "Unknown,";
							}
							print OUT $split[0] . "," . $exon . "," . substr($split[1], 3) . "," . $start . "," . $stop . "," . $split[4] . "\n";
							$exon++;
						}
					}
				}
			}
		} else{
			print STDERR "Still to figure out what else is possible\n";
		}
		close OUT;
		$table_name = $organism . "_genes";
		if(&check_if_table_exists($table_name) == 1) {
			print STDERR "Table $table_name exists\n";
			print STDERR "Aborting...\n";
			$dbh->do("DROP TABLE $table_name");
		}
		$command = "CREATE TABLE " . $table_name . " (id varchar(40) not null, transcript varchar(20) not null, exon integer not null, chr varchar(20) not null, start bigint not null, stop bigint not null, strand varchar(1) not null)";
		$dbh->do($command);
		$command = "CREATE INDEX  ON " . $table_name . " (chr)";
		$dbh->do($command);
		$command = "CREATE INDEX ON " . $table_name . "(id)";
		$dbh->do($command);
		$command = "CREATE UNIQUE INDEX ON " . $table_name . "(id, transcript, exon)";
		$dbh->do($command);
		$command = "CREATE INDEX ON " . $table_name . "(id, transcript)";
		$dbh->do($command);
		`bash database_insert.sh $table_name tmp_genes`;
		print STDERR "\nCheck if database insert was successful!\n";
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$file_lines = `wc -l tmp_genes`;
		@s = split('\s+', $file_lines);
		if($num->{'c'} == $s[0]) {
			print LOG "$table_name\n"; 
			$to_delete{"tmp_genes"} = 1
		} else {
			exit;
		}
		print STDERR "Genes successfully added to database!\n\n";
	} else {
		print STDERR "Genes were already successfully loaded into database\n\n";
	}
}

print STDERR "Database successfully initalized!\n\n\n";
close LOG;
print STDERR "Delete temporary files\n";
print STDERR "This may take a while\n";
foreach my $f (keys %to_delete) {
	print $f . "\n";
	`rm $f`;
}
print STDERR "Done\n";
print STDERR "Nothing more to do - finishing!\n";
if($multithreading == 0) {
	threads->exit();
}
#SUBFUNCTIONS
my $merge_line = "";
my @len;
my $end = 0;
my $overlap = 0;
my @current;
my @points;
my @last = (0);
my $min_start;
my $number;
my $real_length_last;
my $real_length_current;
my $final_length;
my $n;

my %convert = ();

sub convert_chr_to_num{
	if(exists $convert{$_[0]}) {
		return $convert{$_[0]};
	}
	if($_[0] =~ /^[0-9]+$/) {
		$convert{$_[0]} = $_[0];
		return $convert{$_[0]};
	} else {
		#To make sure that we do not save a number we will have later as chromosome we just add 100
		$convert{$_[0]} = (keys %convert) + 100;
		return $convert{$_[0]};
	}
	return 0;
}

sub check_if_table_exists{
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	return $num->{'c'};
}

#sub check_multithreading_support{
#	#First check if multithreading is not turned off by #user	
#	if($multithreading == 1) #{
#		return 1#;
#	}
#	#Now check if multithreading is compiled into perl
#	my $a = `perl -V | grep useithreads`;
#	if(length($a) > 0) {
#		#Threading is possible - check cores
#		my $number = `nproc`;
#		chomp $number;
#		if($core_user == 0) {
#			return int($number/4);
#		}
#		if($core_user > $number) {
#			print STDERR "User specified more cores than exisiting!\n";
#			print STDERR "Use 1/4 of existing cores\n";
#			return int($number/4);
#		} else {
#			return $core_user;
#		}
#	} else {
#		return 1;
#	}
#	return 0;
#}

sub ref_genome{
	my $chr;
	my $working_seq = "";
	my $first_100 = "";
	my $pos = 0;
	my $save;
	for(my $i = 0; $i < @_; $i++) {
		chomp $_[$i];
		open FH, "<$_[$i]";
		foreach my $line (<FH>) {
			chomp $line;
			if(substr($line, 0, 1) eq ">") {
				$n = () = $working_seq =~ /N/gi;
				if($n < length($working_seq)) {
					$save = $pos . "," . uc($working_seq) . "," . $chr;
					push(@ref, $save);
				}
				$working_seq = "";
				$save = "";
				$chr = substr($line, 4);
				print STDERR "\tchromosome " . $chr . "\n";
				$pos = 0;
				next;
			}
			$working_seq .= $line;
			while(length($working_seq) >= 100) {
				$first_100 = substr($working_seq, 0, 100);
				$working_seq = substr($working_seq, 100);
				$n = () = $first_100 =~ /N/gi;
				if($n < 100) {
					$save = $pos . "," . uc($first_100) . "," . $chr;
					push(@ref, $save);
				}
				$pos++;
			}
		}
		$pos = 0;
		$first_100 = "";
		$save = "";
	}	
}

sub merge{
	my $line = $_[0];
	if($line eq "") {
		return;
	}
	chomp $line;
	#split line to get information for every strain
	@current = split('\t', $line);
	#Save the original annotation to know if there was a . or not and substitute also all . with the reference sequence
	@points = split('\t', $line);
	$merge_line = $current[0] . "," . $current[1] . "," . $current[2];
	for(my $i = 3; $i < @current; $i++) {
		$current[$i] = $current[2] if ($current[$i] eq ".");
		$merge_line = $merge_line . "," . $current[$i];
	} 
	#Save the line look at before, to check if there is an overlap and if it is necessary to merge the two mutations
	$last_line = $lines[-1] if @lines > 0;
	@last = split(',', $last_line) if @lines > 0;
	if(@lines < 1) {
		push(@lines, $merge_line);
		$end = $current[1] + length($current[2]);
		return 0;
	}
	if($current[0] ne $last[0]) { 
		push(@lines, $merge_line);
		$end = 0;
		$overlap = 0;
#		print STDERR "\tchromosome " . $current[0] . "\n";
		return;
	}
	#Start checking for different mutation overlaps
	#Case 1: Both mutations start at the same genomic position
	if($current[1] == $last[1]) {
		#Remove the last line from printing array
		pop(@lines);
		$merge_line = $last[0] . "," . $last[1];
		#The current mutation is longer than the previous one and both start at the same position
		# That means, the current mutation is a SNP. If it were an indel, the annotation in the indel file would be one long indel annotation, not two seperate
		# Therefore, just change the first base of the indel, so the SNP is regarded too
		if(length($current[2]) > length($last[2])) {
			for(my $i = 2; $i < @current; $i++) {
				$merge_line = $merge_line . "," . $last[$i] . substr($current[$i], length($last[$i]), length($current[$i]) - length($last[$i]));
			}
		} else {
		#In this case the current mutation is longer than the previous one which mean the previous one is the indel and the current one is the SNP
			for(my $i = 2; $i < @current; $i++) {
				if(length($last[$i]) < length($current[$i])) { exit; }
				$merge_line = $merge_line . "," . $current[$i] . substr($last[$i], length($current[$i]), length($last[$i]) - length($current[$i]));
			}
		}
		@len = split(',', $merge_line);
		$end = $last[1] + length($len[2]) - 1;
		push(@lines, $merge_line);
	#The new mutation starts at the end of the previous one
	#Possibilities: The last position has a SNP - in this case just substitute the last position as long as it is no gap
	# 		At the last position there starts an indel - substitute the last position as long as there is no gap
	#		If there is a gap at the last position, just add as many gaps as the new mutation is long
	} elsif($current[1] == $end) {
		pop(@lines);
		$merge_line = $last[0] . "," . $last[1];
		for(my $i = 2; $i < @current; $i++) {
			if(substr($last[$i], length($last[$i]) - 1, 1) eq "-") {
				#That should not happen because then either our calculation of the end is wrong or something in the logic of this script is completely screwed up. So in case this happens stop the script
				$merge_line = $merge_line . "," . $last[$i];
			} else {
				#The last position of the last mutation and the first position of the new mutation are different 
				if(substr($last[$i], length($last[$i]) - 1, 1) ne substr($current[$i], 0, 1)) {
					#In the original mutation file is a ., which means that there is no information about this mutation. We always assume that there is the reference, but in this case we assume that the new annotated mutation is more correct and we use the new mutation at this position
					if($points[$i] eq ".") {
						$merge_line = $merge_line . "," . $last[$i] . substr($current[$i], 1, length($current[2]) - 1);
					#There is no point in the current mutation. Therefore, we assume that the previous mutation annotation was incorrect and we use the new mutation
					} else {
						$merge_line = $merge_line . "," . substr($last[$i], 0, length($last[$i]) - 1) . $current[$i];
					}
				#The position at the end of the previous mutation is equal to the start of the new mutation - we must merge the two mutations		
				} else {
					$merge_line = $merge_line . "," . substr($last[$i], 0, length($last[$i]) - 1) . $current[$i];
				}
			}
		}
		my @test = split(',', $merge_line);
		@len = split(',', $merge_line);
		$end = $last[1] + length($len[2]) - 1;
		push(@lines, $merge_line);
	#There is an overlap somewhere in the middle
	} elsif($current[1] < $end) {
		pop(@lines);
		$merge_line = $last[0] . "," . $last[1];
		#Overlap of the mutations
		$overlap = $current[1] - $last[1];
		$real_length_current = length($current[2]);
	
		$real_length_last = length($last[2]);
		#Calcualte real length of the new merge
		$min_start = $last[1];
		$final_length = length($last[2]);
		if($current[1] + $real_length_current > $last[1] + $real_length_last) {
			$final_length = $final_length + ($current[1] + $real_length_current) - ($last[1] + $real_length_last);
		}
		$number = () = $current[2] =~ /-/gi;
		$final_length = $final_length + $number;
		#Now merge the two mutations
		#Check number of gaps in last, because it is possible that we have already merged several indels, so it is possible to have gaps in the last line
		for(my $i = 2; $i < @last; $i++) {
			#We start with the old mutation, because we already know that the old mutation starts earlier in the string than the new now (if they would start at the same position we would go in another if and we know that the mutation file is sorted)
			$merge_line = $merge_line . "," . substr($last[$i], 0, $overlap);
			#After adding the overlapping region we look at the position in the old mutation where the new mutation starts
			#There is a gap - the new mutation can not appear here, so just fill the string with gaps
			# We can fill the string with gaps to the final length, because we know that we just use single mutations. So if there is a gap in the string, it can just be followed by gaps, because an insertion takes place in another strain. 
			#We just do not do anything
			if(($overlap + 1 <= length($last[$i])) || length($last[$i]) >= $overlap + 1) {
				#do nothing
			#there is no special case, so the new mutation is just added to the string at this position
				$merge_line = $merge_line . $current[$i];
			}
			#Now we need to add the rest of the longer mutation at the end of the string, but just if the mutation was added into the string and not gaps - because in the case of gaps we have already filled the string with gaps to max length. 
			#In case the second mutation was actually inserted into the string we need to check if the real length of the first mutation is longer than the real length of the second mutation and the beginning of the first deletion
			# If this is the case it means that the second mutation is a deletion in the first mutation
			#Example:
			#Mutation 1: 100 ACCGTG		Mutation 2: 102 C--- 
			#Result: 100 ACC---GTG
			#Straing 1:  ACC---GTG		Straing 2: ACCTTTGTG
			# If the length of the second mutation and the overlap is greater than the length of the first mutation it means that the second mutation is a deletion in the first mutation, so we do not need to add the rest of the first mutation to the string
			#Example:
			#Mutation 1: 100 ACCGTG		Mutation 2: 102: CGTG
			#Result: 100 ACCGTG 
			#Straing 1: ACCGTG		Strain 2: ACC---
			if($final_length > ($real_length_current + $overlap + $number)) {
				if(length($last[$i]) >= $overlap + 1 && substr($last[$i], $overlap, 1) ne "") {
					$merge_line = $merge_line . substr($last[$i], length($last[$i]) - ($final_length - ($real_length_current + $overlap + $number)));
				}
			}
		}
		my @test = split(',', $merge_line);
		$end = $test[1] + length($test[2]) - 1;
		push(@lines, $merge_line);
		$overlap = 0;
	} else {
		$merge_line = $current[0] . "," . $current[1];
		for(my $i = 2; $i < @current; $i++) {
			$merge_line = $merge_line . "," . $current[$i];
		}
		push(@lines, $merge_line);
		$end = $current[1] + length($current[2]) - 1;
	}
}

#First step filter mutations if just one file is given
#If SNP and Indel file are given - merge the files
sub read_file_line {
	my $fh = shift;
	my @split;

	if ($fh and my $line = <$fh>) {
		chomp $line;
#		return [ split('\t', $line) ];
		return $line;
	}
	return;
}

sub print {
	my $string =  "";
	if($filter == 1 && $_[0]->[6] ne "PASS") {
		return $string;
	}
	my @alter;
	my @allele;
	$string .= $_[0]->[0] . "\t" . $_[0]->[1] . "\t" . $_[0]->[3];
	#Save the mutations
	my @mut = split(',', $_[0]->[4]);
	unshift @mut, ".";

	for(my $i = 9; $i < scalar(@{$_[0]}); $i++) {
		if($_[0]->[$i] eq ".") {
			$string .= "\t" . $_[0]->[$i];
		} else {
			@alter = split(":", $_[0]->[$i]);
			@allele = split('/', $alter[0]);
			if($allele[0] eq $allele[1]) {
				$string .= "\t" . $mut[$allele[0]];
			} else {
				if($same == 1) {
					$string = "";
					return $string;
				} else {
					$string .= "\t" . $mut[$allele[0]];
				}
			}
		}
	}
	return $string;
}

my $query;
my $ref_seq;
my $old_pos;
my $old_seq;

sub shift_vector{
	my %total_chr = config::chromosome_number($organism);
	my $s = $_[0];
	open OUT, ">shift_vector_$s";
#	print OUT "ref pos\tstrain pos\tvector r->s\tchromosome\n";
	#Add here something from config for different genomes!!!!
	foreach my $c (keys %total_chr) {
#	for(my $i = 1; $i <= $total_chr; $i++) {
#		if($i == $total_chr) {
#			$c = "X";
#		} elsif($i == ($total_chr - 1)) {
#			$c = "Y";
#		} else {
#			$c = $i;
#		}
		$query = 0;
		$ref_seq = 0;
		$command = "SELECT * FROM " . $organism . "_mutations WHERE length(reference) != length($s) AND chromosome = \'$c\' ORDER  by position";
		$sth = $dbh->prepare($command);
	#	$sth = $dbh->prepare("SELECT * FROM $organism_mutations WHERE length(reference) != length($s) AND chromosome = \'$c\' ORDER by position");
		$sth->execute();
		while (my $ref = $sth->fetchrow_hashref()) {
                        $old_pos = $ref->{'position'};
                        $old_seq = $ref->{'reference'};
			#iterate over all positions that are the same
			while($ref_seq < $ref->{'position'}) {
				$ref_seq++;
				$query++;
				
			}
			#Reference longer than mutated strain
			if(length($ref->{'reference'}) > length($ref->{$s})) {
				if(length($ref->{$s}) > 1) {
					#we are already at the first position of this deletion, so just add length(strain) - 1
					for(my $i = 1; $i < length($ref->{$s}); $i++) {
						$query++;
						$ref_seq++;
					}
					#Next we jump the in the reference sequence to the end of the insertion in the ref (deletion in strain), therefore we have to add 1 (because we look at the next position) and the lenght of the difference
					$ref_seq++;
					$query++;
					$ref_seq = $ref_seq + (length($ref->{'reference'}) - length($ref->{$s}));
					print OUT $ref_seq . "," . $query . "," . ($query - $ref_seq) . "," . $c . "\n";
				} else {
					#We are already at the current position of the mutation
					$query++;
					$ref_seq++;
#Just delete
					$ref_seq = $ref_seq + (length($ref->{'reference'}) - length($ref->{$s}));
					print OUT $ref_seq . "," . $query . "," . ($query - $ref_seq) . "," . $c . "\n";
				}
			#Reference is shorter than mutated strain
			} else {
				if(length($ref->{'reference'}) > 1) {
					#Add the overlap. We are already at the first postiion of the insertion, so we just add -1 less
					for(my $i = 0; $i < length($ref->{'reference'}); $i++) {
						$query++;
						$ref_seq++;
					}
					#Add the positions of the query and do not change the position of the ref
					for(my $i = 0; $i < (length($ref->{$s}) - length($ref->{'reference'})); $i++) {
						$query++;
					}
					print OUT $ref_seq . "," . $query . "," . ($query - $ref_seq) . "," . $c . "\n"; 
				} else {
					#First position is common - go to next position
					$query++;
					$ref_seq++;
					#Common position is already done, just add insertion
					for(my $i = 0; $i < (length($ref->{$s}) - length($ref->{'reference'})); $i++) {
						$query++;
					}
					print OUT $ref_seq . "," . $query . "," . ($query - $ref_seq) . "," . $c . "\n";
				}	
			}
		}
	}
	close OUT;

	#Create table
	$table_name = "offset_" . $organism . "_" . $s; 
	if(&check_if_table_exists($table_name) == 1) {
		$dbh->do("DROP TABLE " . $table_name);
	}
	$command = "CREATE TABLE " . $table_name . " (ref_pos bigint not null, strain_pos bigint not null, shift bigint not null, chr varchar(20) not null, PRIMARY KEY (ref_pos, chr))";
	$dbh->do($command);
	$command = "CREATE UNIQUE INDEX ON " . $table_name . " (strain_pos, chr)";
	$dbh->do($command);
	`bash database_insert.sh $table_name shift_vector_$s`;
	print STDERR "DONE WITH $s\n";
	$to_delete{"shift_vector_$s"} = 1;
}
