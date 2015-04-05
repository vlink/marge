#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';
require '../general/system_interaction.pm';
use threads;
use Data::Dumper;

$_ = "" for my($organism, $snp, $indel, $both, $genome_folder, $genome_file, $annotation, $genes, $name, $current_line, $last_line);
$_ = 0 for my($filter, $same, $core_user, $multithreading, $resume, $help, $check_muts, $check_ref, $check_annotation, $check_genes, $change);
my($num, @lines, @merge_line, @end, %hash, $dbh, $sth, $header, %check_offset, @split, $snps, $indels, @s, @ref, $ref, @i, $table_name, $f1, $f2, $command, @h, %check_strains, %strains_to_use);

for(my $i = 0; $i < 18; $i++) {
	my @a = ();
	push(@lines, \@a);
}

sub printCMD{
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
	print STDERR "\n\n";
	print STDERR "-h | --help: shows help\n\n\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
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
		"resume" => \$resume, 
		"h" => \$help,
		"help" => \$help)
or &printCMD(); 

if($help == 1) {
	&printCMD();
}

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
			if($line eq $organism . "_ref_genome") {
				$check_ref = 1;
				print STDERR "Reference genome was already added\n";
			}
			if(substr($line, 0, 6) eq "strain") {
				@split = split('_', $line);
				$check_strains{$split[-1]} = 1;
			}
			if($line eq $organism . "_annotations") {
				$check_annotation = 1;
				print STDERR "Annotations were already added\n";
			}
			if($line eq $organism . "_genes") {
				$check_genes = 1;
				print STDERR "Genes were already added\n";
			}
		}
		close IN;
		open LOG, ">>database_init.log";
	}
} else {
	open LOG, ">database_init.log";
	print LOG "organism:" . $organism  . "\n";
}

if($genome_folder eq "" && $genome_file eq "") {
	print STDERR "The reference genome is missing!\n";
	exit;
}

if(system_interaction::check_module_installed("DBIx::Threaded") == 1) { $multithreading = 1; }

if(($snp ne "" || $indel ne "") && $both ne "") {
	print STDERR "You can not use this combination of files!\n";
	print STDERR "Please specify either two files, one containing SNPs, the other containing Indels\nor one file contaning both!\n";
	exit;
}

my $con = config::database_connection();
if($multithreading == 1) {
	$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
} else {
	$dbh = DBIx::Threaded->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
}

if($check_ref == 0) {
	print STDERR "Add reference genome!\n";
	my @files;

	$table_name = $organism . "_ref_genome";
	if(&check_if_table_exists($table_name) == 1) {
		print STDERR "Table " . $table_name . " already exists!\n";
		print STDERR "Aborting...\n";
		$dbh->do("DROP TABLE " . $table_name);
	#	exit;
	}

	#Create table
	$dbh->{ Warn } = 0;
	$command = "CREATE TABLE " . $table_name . " (pos bigint not null, seq varchar(100) not null, chr varchar(20) not null, PRIMARY KEY(pos, chr))";
	$dbh->do($command);
	$dbh->{ Warn } = 1;

	if($genome_folder ne "") {
		@files = `ls -1 $genome_folder/*`;
		&ref_genome(@files);
	} 
	if($genome_file ne "") {
		$files[0] = $genome_file;
		&ref_genome(@files);
	}

	print STDERR "\n\nInserting reference genome into database\n";
	print STDERR "This may take a while\n\n";
	$dbh->{ AutoCommit } = 0;
	$sth = $dbh->prepare("INSERT INTO $table_name (pos, seq, chr) VALUES (?, ?, ?)");
	for(my $i = 0; $i < @ref; $i++) {
		$sth->execute(@{ $ref[$i] });
	}
	$sth->finish();
	$dbh->commit();
	$dbh->{ AutoCommit } = 1;

#	print STDERR "Check if database insert was successful for genome!\n";
#	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
#	$sth->execute();
#	$num = $sth->fetchrow_hashref();
#	$sth->finish();
#	if(@ref == $num->{'c'}) {
		print LOG $table_name . "\n";
#	} else {
#		print STDERR "Input in database not successful!\n";
#		exit;
#	}
	print STDERR "Successful!\n\n";
}

#define header
if($snp ne "") {
	$header = `grep -m 1 \'#CHROM\' $snp`;
} elsif($indel ne "") {
	$header = `grep -m 1 \'#CHROM\' $indel`;
} elsif($both ne "") {
	print STDERR "hmpf\n";
	$header = `grep -m 1 \'#CHROM\' $both`;
}
chomp $header;
@h = split('\t', $header);
for(my $i = 9; $i < @h; $i++) {
	if($h[$i] =~ /^\d/) {
		$h[$i] = "s" . lc($h[$i]);
	}
	$h[$i] = lc($h[$i]);
}

#Merge SNP and indel file if both exists
if(keys %check_strains == 0) {
	#open files and filter out header
	if($both ne "") {
		$snp = $both;
	}
	if($snp ne "") {
		open($f1, $snp);
		$snps = read_file_line($f1);
		while(substr($snps, 0, 1) eq "#") {
			$snps = read_file_line($f1);
			next;
		}
		@s = split('\t', $snps);
	}
	if($indel ne "") {
		open($f2, $indel);
		$indels = read_file_line($f2);
		while(substr($indels, 0, 1) eq "#") {
			$indels = read_file_line($f2);
			next;
		}
		@i = split('\t', $indels);
	}
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
			&merge($current_line, $last_line);
			$indels = read_file_line($f2);		
			$change = 1;
		}
		if($i_c > $s_c) {
			if($change == 0) {
				print STDERR "\tchromosome " . $i[0] . "\n";
			}
			$current_line = $check_snp;
			&merge($current_line, $last_line);
			$snps = read_file_line($f1);
			$change = 1;
		}

		if($s_c == $i_c) {
			$change = 0;
			if($s[1] < $i[1]) {
				$current_line = $check_snp;
				&merge($current_line, $last_line);
				$snps = read_file_line($f1);
			} else {
				$current_line = $check_indel;
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

	#Find the right order to start the threads
	for(my $i = 9; $i < @h; $i++) {
		$hash{$i} = @{$lines[$i - 9]};
		$strains_to_use{$i} = 1;
	}
}
 else {
	for(my $i = 9; $i < @h; $i++) {
		if(!exists $check_strains{$h[$i]}) {
			$strains_to_use{$i} = 1;
		}
	}
	print STDERR "" . (keys %strains_to_use) . " strains will be added\n";
}

my $core = system_interaction::check_multithreading_support($multithreading, $core_user);

if(keys %strains_to_use > 0) {
	print STDERR "Create shifting vectors!\n";
}
if($core == 1) {
	foreach my $key (keys %strains_to_use) {
		print STDERR "vector for " . $h[$key] . "\n";
		&shift_vector($h[$key]);
	}
} else {
	my @running = ();
	my @Threads;
	my $current_thread;
	while (scalar @Threads < (keys %strains_to_use)) {
		@running = threads->list(threads::running);
		if (scalar @running < $core) {
			#Find out which thread we need to start;
			foreach my $keys (sort { $hash{$a} <=> $hash{$b} } keys %hash) {
				$current_thread = $keys;
			}
			delete $hash{$current_thread};
			my $thread = threads->new( sub { shift_vector($h[$current_thread], $current_thread) });
		#	sleep(4);
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


$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

my $file_lines;
if($annotation ne "") {
	if($check_annotation == 0) {
		print STDERR "Loading annotations into database!\n";
	
		$table_name = $organism . "_annotations";
		if(&check_if_table_exists($table_name) == 1) {
			print STDERR "Table " . $table_name . " exists!\n";
			print STDERR "Aborting...\n";
			$dbh->do("DROP TABLE $table_name");
		}
		$dbh->{ Warn } = 0;
		$command = "CREATE TABLE " . $table_name . "(chr varchar(20) not null, start bigint not null, stop bigint not null, strand varchar(1) not null, annotation varchar(64) not null, gene varchar(116) not null)";
		$dbh->do($command);
		$command = "CREATE INDEX  ON " . $table_name . " (gene)";
		$dbh->do($command);
		$command = "CREATE INDEX  ON " . $table_name . " (gene, annotation)";
		$dbh->do($command);
		$dbh->{ Warn } = 1;

		$dbh->{ AutoCommit } = 0;
		$sth = $dbh->prepare("INSERT INTO $table_name (chr, start, stop, strand, annotation, gene) VALUES (?, ?, ?, ?, ?, ?)");

		open FH, "<$annotation";
		my @ann;
		my $entries = 0;

		foreach my $line (<FH>) {
			chomp $line;
			@split = split('\t', $line);
			push(@ann, substr($split[1], 3));
			for(my $i = 2; $i < 6; $i++) {
				push(@ann, $split[$i]);
			}
			if($split[5] eq "TTS") {
				push(@ann, substr($split[0], 5, length($split[0]) - 6));
			} elsif($split[5] eq "E" || $split[5] eq "I" || $split[5] eq "3UTR" || $split[5] eq "5UTR" || $split[5] eq "P") {
				@s = split(",", $split[0]);
				@split = split('\(', $s[0]);
				push(@ann, $split[1]);
			} else {
				push(@ann, "");
			}
			$sth->execute(@ann);
			@ann = ();
			$entries++;
		}
		close FH;
		$sth->finish();
		$dbh->commit();
		$dbh->{ AutoCommit } = 1;

		print STDERR "\nCheck if database insert was successful!\n";
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();

		if($num->{'c'} == $entries) {
			print LOG "$table_name\n"; 
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
		my @genes;

		$dbh->{ Warn } = 0;
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
		$dbh->{ Warn } = 1;

		$dbh->{ AutoCommit } = 0;
		$sth = $dbh->prepare("INSERT INTO $table_name (id, transcript, exon, chr, start, stop, strand) VALUES (?, ?, ?, ?, ?, ?, ?)");
	
		my $occ = 0;
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
								push(@genes, $save_names{$split[0]});
							} else {
								push(@genes, "Unknown");
							}
							push(@genes, $split[0]);
							push(@genes, $exon);
							push(@genes, substr($split[1], 3));
							push(@genes, $start);
							push(@genes, $stop);
							push(@genes, $split[4]);
							$sth->execute(@genes);
							@genes = ();	
							$occ++;
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
								push(@genes, $save_names{$split[0]});
							} else {
								push(@genes, "Unknown");
							}
							push(@genes, $split[0]);
							push(@genes, $exon);
							push(@genes, substr($split[1], 3));
							push(@genes, $start);
							push(@genes, $stop);
							push(@genes, $split[4]);
							$sth->execute(@genes);
							$occ++;
							@genes = ();
							$exon++;
						}
					}
				}
			}
			$sth->finish();
			$dbh->commit();
			$dbh->{ AutoCommit } = 1;
		} else {
			print STDERR "Still to figure out what else is possible\n";
		}
		print STDERR "\nCheck if database insert was successful!\n";
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM $table_name");
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == $occ) {
			print LOG "$table_name\n"; 
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
	$sth->finish();
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
	my $save = [];
#	$dbh->{ AutoCommit } = 0;
#	$sth = $dbh->prepare("INSERT INTO $table_name (pos, seq, chr) VALUES (?, ?, ?)");
#	for(my $i = 0; $i < @ref; $i++) {
#		$sth->execute(@{ $ref[$i] });
#	}
#	$sth->finish();
#	$dbh->commit();
#	$dbh->{ AutoCommit } = 1;

	for(my $i = 0; $i < @_; $i++) {
		chomp $_[$i];
		open FH, "<$_[$i]";
		foreach my $line (<FH>) {
			chomp $line;
			if(substr($line, 0, 1) eq ">") {
				$n = () = $working_seq =~ /N/gi;
				if($n < length($working_seq)) {
			#		$sth->execute($pos, uc($working_seq), $chr);
					push(@{ $save }, $pos);
					push(@{ $save }, uc($working_seq));
					push(@{ $save }, $chr);
					push(@ref, $save);
					$save = [];
				}
				$working_seq = "";
				$chr = substr($line, 4);
				print STDERR "\tchromosome " . $chr . "\n";
			#	$sth->finish();
			#	$dbh->commit();
				$pos = 0;
				next;
			}
			$working_seq .= $line;
			while(length($working_seq) >= 100) {
				$first_100 = substr($working_seq, 0, 100);
				$working_seq = substr($working_seq, 100);
				$n = () = $first_100 =~ /N/gi;
				if($n < 100) {
				#	$sth->execute($pos, uc($first_100), $chr);
					push(@{ $save }, $pos);
					push(@{ $save }, uc($first_100));
					push(@{ $save }, $chr);
					push(@ref, $save);
					$save = [];

				}
				$pos++;
			}
		}
		$pos = 0;
		$first_100 = "";
	}	
#	print STDERR "Commit database\n";
#	$sth->finish();
#	$dbh->commit();
#	$dbh->{ AutoCommit } = 1;
}

sub merge{
	my $line = $_[0];
	if($line eq "") {
		return;
	}
	chomp $line;
	@current = split('\t', $line);
	#Save reference and position
	my $pos = $current[1];
	my $ref = $current[2];
	my $merge_line;
	my @last;
	for(my $i = 3; $i < @current; $i++) {
		@merge_line = ();
		#Check if reference and mutations are different
		if($current[1] eq $current[$i] || $current[$i] eq ".") {
			next;
		} else {
			if(length($current[2]) > 1 && length($current[$i]) > 1) {
				#Remove the first bases if they are the same
				while(substr($current[2], -1) eq substr($current[$i], -1) && length($current[2]) > 1 && length($current[$i]) > 1) {
					chop $current[2];
					chop $current[$i];
				}
				#Now start from the beginning
				my $max = (length($current[2]) > length($current[$i])) ? length($current[$i]) : length($current[2]);	
				for(my $l = 0; $l < $max - 1; $l++) {
					if(length($current[2]) > 1 && length($current[$i]) > 1 && substr($current[2], 0, 1) eq substr($current[$i], 0, 1)) {
						$current[1]++;
						$current[2] = unpack "xA*", $current[2];
						$current[$i] = unpack "xA*", $current[$i];
					}
				}
			}
			#Check if there was already a mutation before this one
			@last = split(',', $lines[$i-3]->[-1]) if @{$lines[$i-3]} > 0;
			if(@{$lines[$i-3]} < 1 || $last[0] ne $current[0]) {
				$merge_line = $current[0] . "," . $current[1] . "," . $current[2] . "," . $current[$i];
				push(@{ $lines[$i-3] }, $merge_line);
				$current[1] = $pos;
				$current[2] = $ref;
				next;
			}
			#Check if there is an overlap between the last and the current mutation
			if($current[1] <= $last[1] + length($last[2]) - 1) {
				#Deletion in the strain so that the new mutation can not exist = ignore second mutation
				if($last[1] + length($last[3]) -1 < $current[1]) {
					$current[1] = $pos;
					$current[2] = $ref;
					next;
				#Start at the same position - Ignore second mutation
				} elsif($last[1] == $current[1]) {
					$current[1] = $pos;
					$current[2] = $ref;
					next;
				} else {
					$merge_line = $last[0] . "," . $last[1] . ",";
					#The current mutation happens within the last mutation
					$merge_line .= substr($last[2], 0, $current[1] - $last[1]);;
					$merge_line .= $current[2];
					$merge_line .= substr($last[2], $current[1] - $last[1] + length($current[2])) . ",";
					$merge_line .= substr($last[3], 0, $current[1] - $last[1]);
					$merge_line .= $current[$i];
					$merge_line .= substr($last[3], $current[1] - $last[1] + length($current[$i]));
					pop(@{ $lines[$i-3]});
					push(@{ $lines[$i-3]}, $merge_line);
					exit;
				}
			} else {
				$merge_line = $current[0] . "," . $current[1] . "," . $current[2] . "," . $current[$i];
				push(@{ $lines[$i-3]}, $merge_line);
				$current[1] = $pos;
				$current[2] = $ref;
				next;

			}
		}
		$current[1] = $pos;
		$current[2] = $ref;
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

	my $dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
	my $sth;
	my %total_chr = config::chromosome_number($organism);
	my $s = $_[0];
	my $i = $_[1];
	$table_name = $organism . "_mutations_" . $s;
	#Check if table already exists
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();
	if($num->{'c'} == 1) {
		print STDERR "Table " . $table_name . " already exists!\n";
		print STDERR "Aborting....\n";
		$dbh->do("DROP TABLE " . $table_name);
	#	exit;
	}
	#Create table
	$dbh->{ Warn } = 0;
	$command = "CREATE TABLE " . $table_name . "(chr varchar(2) not null, pos bigint not null, reference varchar(116) not null, strain varchar(116) not null, index SERIAL PRIMARY KEY)";
	$dbh->do($command);
	$command = "CREATE INDEX  ON " . $table_name . " (chr, pos)";
	$dbh->do($command);
	$dbh->{ Warn } = 1;
	print STDERR "\n\nInserting mutations into database for $s\n";
	print STDERR "This may take a while\n";
	$dbh->{ AutoCommit } = 0;
	$sth = $dbh->prepare("INSERT INTO $table_name (chr, position, reference, strain) VALUES (?, ?, ?, ?)");
	foreach my $l (@{$lines[$i - 9]}) {
		my @a = split(',', $l);
		$sth->execute(@a);
	#	print OUT $l . "\n";
	}
	$sth->finish();
	$dbh->commit();
	$dbh->{ AutoCommit } = 1;
#	`bash database_insert.sh $table_name tmp_merged_$h[$i].txt`;

	#Add here something from config for different genomes!!!!
#	my $c;
	my @vector;
	my $save = [];
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
		$command = "SELECT * FROM " . $organism . "_mutations_" . $s . " WHERE length(reference) != length(strain) AND chr = \'$c\' ORDER  by position";
		$sth = $dbh->prepare($command);
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
			if(length($ref->{'reference'}) > length($ref->{'strain'})) {
				if(length($ref->{'strain'}) > 1) {
					#we are already at the first position of this deletion, so just add length(strain) - 1
					for(my $i = 1; $i < length($ref->{'strain'}); $i++) {
						$query++;
						$ref_seq++;
					}
					#Next we jump the in the reference sequence to the end of the insertion in the ref (deletion in strain), therefore we have to add 1 (because we look at the next position) and the lenght of the difference
					$ref_seq++;
					$query++;
					$ref_seq = $ref_seq + (length($ref->{'reference'}) - length($ref->{'strain'}));
					push($save, $ref_seq);
					push($save, $query);
					push($save, $query - $ref_seq);
					push($save, $c);
					push(@vector, $save);
					$save = [];
				} else {
					#We are already at the current position of the mutation
					$query++;
					$ref_seq++;
#Just delete
					$ref_seq = $ref_seq + (length($ref->{'reference'}) - length($ref->{'strain'}));
					push($save, $ref_seq);
					push($save, $query);
					push($save, $query - $ref_seq);
					push($save, $c);
					push(@vector, $save);
					$save = [];
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
					for(my $i = 0; $i < (length($ref->{'strain'}) - length($ref->{'reference'})); $i++) {
						$query++;
					}
					push($save, $ref_seq);
					push($save, $query);
					push($save, $query - $ref_seq);
					push($save, $c);
					push(@vector,$save);
					$save = [];
				} else {
					#First position is common - go to next position
					$query++;
					$ref_seq++;
					#Common position is already done, just add insertion
					for(my $i = 0; $i < (length($ref->{'strain'}) - length($ref->{'reference'})); $i++) {
						$query++;
					}
					push($save, $ref_seq);
					push($save, $query);
					push($save, $query - $ref_seq);
					push($save, $c);
					push(@vector, $save);
					$save = [];
				}	
			}
		}
		$sth->finish();
	}

	#Create table
	$table_name = "offset_" . $organism . "_" . $s; 
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();
	if($num->{'c'} == 1) {
		print STDERR "Table " . $table_name . " already exists!\n";
		print STDERR "Aborting....\n";
		$dbh->do("DROP TABLE " . $table_name);
	}

	$dbh->{ Warn } = 0;
	$command = "CREATE TABLE " . $table_name . " (ref_pos bigint not null, strain_pos bigint not null, shift bigint not null, chr varchar(20) not null, PRIMARY KEY (ref_pos, chr))";
	$dbh->do($command);
	$command = "CREATE UNIQUE INDEX ON " . $table_name . " (strain_pos, chr)";
	$dbh->do($command);
	$dbh->{ Warn } = 1;
	$dbh->{ AutoCommit } = 0;
	$sth = $dbh->prepare("INSERT INTO " . $table_name . " (ref_pos, strain_pos, shift, chr) VALUES (?, ?, ?, ?)");
	for(my $i = 0; $i < @vector; $i++) {
		$sth->execute(@{ $vector[$i] });
	}
	$sth->finish();
	$dbh->commit();
	$dbh->{ AutoCommit } = 1;
	print LOG "strain_" . $organism . "_" . $s . "\n";
	$dbh->disconnect();
}
