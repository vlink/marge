#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';
require '../general/system_interaction.pm';
use threads;
use Data::Dumper;

$_ = "" for my($genome, $ref_genome, $annotation, $genes, $name, $current_line, $last_line, $snp, $indel, $chr);
$_ = 0 for my($filter, $same, $core_user, $multithreading, $resume, $help, $check_muts, $check_ref, $check_annotation, $check_genes, $change, $homo, $chr_name_length, $length_mut, $length_ref, $force);
my($num, $lines, @merge_line, @end, %hash, $dbh, $sth, $header, %check_offset, @split, $snps, $indels, @s, @ref, $ref, @i, $table_name, $f1, $f2, $command, @h, %check_strains, %strains_to_use, @mut_files, $total_chr, %seen);


sub printCMD{
	print STDERR "Usage:\n";
	print STDERR "General commands:\n";
	print STDERR "\t-genome <genome>\n";
	print STDERR "\t-reference <path to folder or file with reference genome in it split by chromosomes>\n";
	print STDERR "\n\nInput files:\n";
	print STDERR "\t-files <file>: Files with mutations - comma spearated list (max 2 - one file for snps, one file for indels)\n";
	print STDERR "\t-homo - Assumes phenotype is homozygouse - if there are two different variations reported just the first one is taken into acocunt\n";
	print STDERR "\n\n";
	print STDERR "Filter vcf file.\n";
	print STDERR "\t-filter - Filters out all mutations that did not pass all filters\n"; 	
	print STDERR "\t-same - Filters out all mutations that are homozyous\n";
	print STDERR "\tif position is homozygous and -same is not defined the first allele is taken without any further evaluation\n";
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
	print STDERR "\t-force: Overwrites existing tables\n";
	print STDERR "\n\n";
	print STDERR "-h | --help: shows help\n\n\n";
	exit;
}

if(@ARGV < 1) {
	&printCMD();
}

my %mandatory = ('-genome' => 1);
my %commandline = map { $_ => 1 } @ARGV;
config::check_parameters(\%mandatory, \%commandline);

#TODO add resume then clean up code and then done
GetOptions(	"genome=s" => \$genome,
		"files=s{,}" => \@mut_files,
		"filter" => \$filter,
		"same" => \$same,
		"reference=s" => \$ref_genome,
		"core=s" => \$core_user,
		"no-multithreading" => \$multithreading,
		"ann=s" => \$annotation,
		"gene=s" => \$genes,
		"names=s" => \$name,
		"resume" => \$resume, 
		"homo" => \$homo,
		"h" => \$help,
		"force" => \$force,
		"help" => \$help)
or &printCMD(); 

$total_chr = config::chromosome_number($genome);
foreach my $chr (keys %{$total_chr}) {
	if($chr eq "NONE") {
		$chr_name_length = 2;
	} else {
		if(length($chr) > $chr_name_length) {
			$chr_name_length = length($chr);
		}
	}
}

if($ref_genome eq "") {
	$check_ref = 1;
	print STDERR "Reference genome sequence was not specified!\n";
}
if(@mut_files == 0) {
	print STDERR "No mutation files specified. Just adding genome to database!\n";
}

for(my $i = 0; $i < 18; $i++) {
	my @a;
	push(@{$lines->[0]}, \@a);
	if($homo == 0) {
		my @b;
		push(@{$lines->[1]}, \@b);
	}
}


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
		open IN, "<database_init_$genome.log";
		foreach my $line (<IN>) {
			chomp $line;
			if(substr($line, 0, 8) eq "genome") {
				if(substr($line, 9) ne $genome) {
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
			if($line eq $genome . "_ref_genome") {
				$check_ref = 1;
				print STDERR "Reference genome was already added\n";
			}
			if(substr($line, 0, 6) eq "strain") {
				@split = split('_', $line);
				$check_strains{$split[-1]} = 1;
			}
			if($line eq $genome . "_annotations") {
				$check_annotation = 1;
				print STDERR "Annotations were already added\n";
			}
			if($line eq $genome . "_genes") {
				$check_genes = 1;
				print STDERR "Genes were already added\n";
			}
		}
		close IN;
		open LOG, ">>database_init_$genome.log";
	}
} else {
	open LOG, ">database_init.log";
	print LOG $genome . "_ref_genome" . "\n";
}

if(system_interaction::check_module_installed("DBIx::Threaded") == 1) { $multithreading = 1; }

my $con = config::database_connection();
if($multithreading == 1) {
	$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
} else {
	$dbh = DBIx::Threaded->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
}

if($check_ref == 0) {
	print STDERR "Add reference genome!\n";
	open OUT, ">genome_human.txt";
	my @files;

	if(-d $ref_genome) {
		@files = `ls -d -1 $ref_genome/*`;
	} else {
		$files[0] = $ref_genome;
	}
	&ref_genome(@files);

	$table_name = $genome . "_ref_genome";
	if(&check_if_table_exists($table_name) == 1) {
		print STDERR "Table " . $table_name . " already exists!\n";
		if($force == 1) {
			print STDERR "Overwrite table\n";
			$dbh->do("DROP TABLE " . $table_name);
		} else {
			print STDERR "Move to next database insert\n";
			last;
		}
	}

	#Create table
	$dbh->{ Warn } = 0;
	$command = "CREATE TABLE " . $table_name . " (pos bigint not null, seq varchar(100) not null, chr varchar($chr_name_length) not null, PRIMARY KEY(pos, chr))";
	$dbh->do($command);
	$dbh->{ Warn } = 1;

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
	close OUT;
	print LOG $table_name . "\n";
	print STDERR "Successful!\n\n";
}



#define header
if(@mut_files > 0) {
	$header = `grep -m 1 \'#CHROM\' $mut_files[0]`;
	chomp $header;
	if($header eq "") {
		print STDERR "Mutation file does not have a header - naming of the strains not possible!\n";
		exit;
	}
	@h = split /[\t\s]+/, $header;
	for(my $i = 9; $i < @h; $i++) {
		if($h[$i] =~ /^\d/) {
			$h[$i] = "s" . lc($h[$i]);
		}
		$h[$i] = lc($h[$i]);
	}
}
#Merge SNP and indel file if both exists
if(@mut_files > 0 && keys %check_strains < (@h - 9)) {
	#open files and filter out header
	$snp = $mut_files[0];
	open($f1, $snp);
	$snps = read_file_line($f1);
	while(substr($snps, 0, 1) eq "#") {
		$snps = read_file_line($f1);
		next;
	}
	@s = split('\t', $snps);

	if(@mut_files == 2) {
		$indel = $mut_files[1];
		open($f2, $indel);
		$indels = read_file_line($f2);
		while(substr($indels, 0, 1) eq "#") {
			$indels = read_file_line($f2);
			next;
		}
		@i = split('\t', $indels);
	}
	print STDERR "Merge files and mutations!\n";

	my $check_snp;
	my $check_indel;
	my $s_c;
	my $i_c;
	print STDERR "\tchromosome 1\n";
	while($snps and $indels) {
		@s = split('\t', $snps);
		$s_c = &convert_chr_to_num($s[0]);
		@i = split('\t', $indels);
		$i_c = &convert_chr_to_num($i[0]);
		if($s_c > $i_c) {
			if($change == 0) {
				print STDERR "\tchromosome " . $s[0] . "\n";
			}
			&merge($indels);
			$indels = read_file_line($f2);		
			$change = 1;
		}
		if($i_c > $s_c) {
			if($change == 0) {
				print STDERR "\tchromosome " . $i[0] . "\n";
			}
			&merge($snps);
			$snps = read_file_line($f1);
			$change = 1;
		}

		if($s_c == $i_c) {
			$change = 0;
			if($s[1] < $i[1]) {
				&merge($snps);
				$snps = read_file_line($f1);
			} else {
				&merge($indels);
				$indels = read_file_line($f2);
			}
		}
	}
	while($snps) {
		&merge($snps);
		$snps = read_file_line($f1);
	}
	while($indels) {
		&merge($indels);
		$indels = read_file_line($f2);
	}

	#Find the right order to start the threads
	for(my $i = 9; $i < @h; $i++) {
		$hash{$i} = @{$lines->[0]->[$i - 9]};
		$strains_to_use{$i} = 1;
	}
} else {
	if(@mut_files == 0) {
		print STDERR "No mutation files specified!\n";
	} else {
		for(my $i = 9; $i < @h; $i++) {
			if(!exists $check_strains{$h[$i]}) {
				$strains_to_use{$i} = 1;
			}
		}
		print STDERR "" . (keys %strains_to_use) . " strains will be added\n";
	}
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
	my $current_thread = "";
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
	
		$table_name = $genome . "_annotations";
		if(&check_if_table_exists($table_name) == 1) {
			print STDERR "Table " . $table_name . " exists!\n";
			if($force == 1) {
				print STDERR "Overwrite table!\n";
				$dbh->do("DROP TABLE $table_name");
			} else {
				print STDERR "Next step!\n";
				last;
			}
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
	$table_name = $genome . "_genes";
	if(&check_if_table_exists($table_name) == 1) {
		print STDERR "Table $table_name exists\n";
		if($force == 1) {
			print STDERR "Overwrite table!\n";
			$dbh->do("DROP TABLE $table_name");
		} else {
			print STDERR "Move to next insert\n";
			$genes = "";
		}
	}
}
my @rem;

if($genes ne "") {
	if($check_genes == 0) {
		my %save_names;
		if($name ne "") {
			open FH, "<$name";
			my($name, $nf) = config::name_format($genome);
			foreach my $line (<FH>) {
				chomp $line;
				@split = split('\t', $line);
				if(@split < $name) {
					next;
				}
				foreach my $k (keys %{$nf}) {
					@rem = split('\.', $split[$k]);
					$save_names{$rem[0]} = $split[$name];
				}	
			}
			close FH;
			print STDERR "Names are saved!\n";
		} else {
			print STDERR "No names for transcripts are specified\n";
		}

		my $gf = config::gene_format($genome);
		my @introns;
		my @parts;
		my @exons;
		my $exon;
		my $start;
		my $stop;
		my $string = "";
		my $blub = 0;
		my @genes;
		my @add;

		$dbh->{ Warn } = 0;
		$command = "CREATE TABLE " . $table_name . " (id varchar(40) not null, transcript varchar(20) not null, exon integer not null, chr varchar($chr_name_length) not null, start bigint not null, stop bigint not null, strand varchar(1) not null)";
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
				@introns = split(',', $split[$gf->{'intron_exon'}]);
				$string = "";
				#Check the strand
				if($split[$gf->{'strand'}] eq "-") {
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
							if($name ne "" && exists $save_names{$split[$gf->{'transcript'}]}) {
								push(@genes, $save_names{$split[$gf->{'transcript'}]});
							} else {
								push(@genes, "Unknown");
							}
							push(@genes, $split[$gf->{'transcript'}]);
							push(@genes, $exon);
							push(@genes, substr($split[$gf->{'chr'}], 3));
							push(@genes, $start);
							push(@genes, $stop);
							push(@genes, $split[$gf->{'strand'}]);
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
								$stop = $split[$gf->{'stop'}];
							}
							if($name ne "" && exists $save_names{$split[$gf->{'transcript'}]}) {
								push(@genes, $save_names{$split[$gf->{'transcript'}]});
							} else {
								push(@genes, "Unknown");
							}
							push(@genes, $split[$gf->{'transcript'}]);
							push(@genes, $exon);
							push(@genes, substr($split[$gf->{'chr'}], 3));
							push(@genes, $start);
							push(@genes, $stop);
							push(@genes, $split[$gf->{'strand'}]);
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
		} elsif(exists $gf->{'intron'}) {
			print STDERR "Intron and exon are listed separate - calculate it!!\n";
			#Make sure lines are not duplicated
	#		`sort -u $genes > $genes.tmp`;
			open FH, "<$genes";
			foreach my $line (<FH>) {
				chomp $line;
				@split = split('\t', $line);
				if(!exists $total_chr->{substr($split[$gf->{'chr'}], 3)}) {
					next;
				}
				if(exists $seen{$split[$gf->{'transcript'}]}) {
					next;
				} else {
					$seen{$split[$gf->{'transcript'}]} = 1;
				}
				@introns = split(',', $split[$gf->{'intron'}]);
				@exons = split(',', $split[$gf->{'exon'}]);
				#Check the strand
				if($split[$gf->{'strand'}] eq "-") {
					$exon = 1;
					for(my $i = 0; $i < @introns; $i++) {
						if($split[$gf->{'start'}] > $introns[$i]) { 
							$introns[$i] = $split[$gf->{'start'}];
						}
						if($split[$gf->{'stop'}] < $exons[$i]) {
							$exons[$i] = $split[$gf->{'stop'}];
						}
						if($name ne "" && exists $save_names{$split[$gf->{'transcript'}]}) {
							push(@genes, $save_names{$split[$gf->{'transcript'}]});
						} else {
							push(@genes, "Unknown");
						}
						push(@genes, $split[$gf->{'transcript'}]);
						push(@genes, $exon);
						push(@genes, substr($split[$gf->{'chr'}], 3));
						push(@genes, $introns[$i]);
						push(@genes, $exons[$i]);
						push(@genes, $split[$gf->{'strand'}]);
					#	print $split[$gf->{'transcript'}] . "\t" . $exon . "\t" . substr($split[$gf->{'chr'}], 3) . "\t" . $introns[$i] . "\t" . $exons[$i] . "\t" . $split[$gf->{'strand'}] . "\n";
						$sth->execute(@genes);
						@genes = ();	
						$exon++;
					}
				} else {
					for(my $i = 0; $i < @introns; $i++) {
						if($split[$gf->{'start'}] > $introns[$i]) { 
							$introns[$i] = $split[$gf->{'start'}];
						}
						if($split[$gf->{'stop'}] < $exons[$i]) {
							$exons[$i] = $split[$gf->{'stop'}];
						}
						if($name ne "" && exists $save_names{$split[$gf->{'transcript'}]}) {
							push(@genes, $save_names{$split[$gf->{'transcript'}]});
					#		print $save_names{$split[$gf->{'transcript'}]} . "\t";
						} else {
							push(@genes, "Unknown");
					#		print "Unknown\t";
						}
						push(@genes, $split[$gf->{'transcript'}]);
						push(@genes, $i + 1);
						push(@genes, substr($split[$gf->{'chr'}], 3));
						push(@genes, $introns[$i]);
						push(@genes, $exons[$i]);
						push(@genes, $split[$gf->{'strand'}]);
						$sth->execute(@genes);
						@genes = ();
					}
				}
			}
			$sth->finish();
			$dbh->commit();
			$dbh->{ AutoCommit } = 1;
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

sub ref_genome{
	my $chr;
	my $working_seq = "";
	my $first_100 = "";
	my $pos = 0;
	my $save = [];
	open CHR, ">>../config/chromosomes.txt";
	print CHR $genome . "\t";

	my $l = 0;	
	for(my $i = 0; $i < @_; $i++) {
		chomp $_[$i];
		open FH, "<$_[$i]";
		foreach my $line (<FH>) {
			chomp $line;
			if(substr($line, 0, 1) eq ">") {
				$n = () = $working_seq =~ /N/gi;
				if($n < length($working_seq)) {
					$l = length($working_seq);
					while($l < 100) {
						$working_seq .= "N";
						$l++;
					}
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
				if(length($chr) > $chr_name_length) {
					$chr_name_length = length($chr);
				}
				print CHR $chr . ",";
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
	#	$pos = 0;
		$first_100 = "";
	}	
	close CHR;
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
	if(substr($current[0], 0, 3) eq "chr") {
		$current[0] = substr($current[0], 3);
	}
	my $pos = $current[1];
	my $ref = $current[3];
	my $merge_line;
	my @last;
	my @allele;
	my @all;
#	$current[4] =~ s/^\./$current[3]/g;
	$current[4] =~ s/\.//g;
	my @variants = split(',', $current[4]);

	unshift @variants, $current[3];

	if($filter == 1 && $current[6] ne "PASS") {
		return 1;
	}


	my $var;
	for(my $i = 9; $i < @current; $i++) {
		@all = split('/', ((split(':', $current[$i]))[0]));
		if($all[0] eq ".") {
			$all[0] = 0;
		}
		if(@all == 1) { 
			$all[1] = 0; 
		} else {
			if($all[1] eq ".") {
				$all[1] = 0;
			}
		}
		if($homo == 1) {
			shift @all;
		}
		@merge_line = ();
		for(my $h = 0; $h < @all; $h++) {
			$var = $variants[$all[$h]];
			$ref = $current[3];
			$pos = $current[1];
			#Check if reference and mutations are different
			if($ref ne $var) {
				if(length($ref) > 1 && length($var) > 1) {
					#Remove the first bases if they are the same
					while(substr($ref, -1) eq substr($var, -1) && length($ref) > 1 && length($var) > 1) {
						chop $ref;
						chop $var;
					}
					#Now start from the beginning
					my $max = (length($ref) > length($var)) ? length($var) : length($ref);	
					for(my $l = 0; $l < $max - 1; $l++) {
						if(length($ref) > 1 && length($var) > 1 && substr($ref, 0, 1) eq substr($var, 0, 1)) {
							$pos++;
							$ref = unpack "xA*", $ref;
							$var = unpack "xA*", $var;
						}
					}
				}
				if(length($ref) > $length_mut) { $length_mut = length($ref); }
				if(length($var) > $length_mut) { $length_mut = length($var); }

				#Check if there was already a mutation before this one
				@last = split(',', $lines->[$h]->[$i-9]->[-1]) if @{$lines->[$h]->[$i-9]} > 0;
				if(@{$lines->[$h]->[$i-9]} < 1 || $last[0] ne $current[0]) {
					$merge_line = $current[0] . "," . $pos . "," . $ref . "," . $var;
					push(@{ $lines->[$h]->[$i-9] }, $merge_line);
					next;
				}
				#Deletion in the strain so that the new mutation can not exist = ignore second mutation
				if($last[1] == $pos) { next; }
                        	if($current[1] < $last[1] + length($last[2])) {
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
					next;
                                }
				#Deletion in the strain so that the new mutation can not exist = ignore second mutation
				$merge_line = $current[0] . "," . $pos . "," . $ref . "," . $var;
				push(@{ $lines->[$h]->[$i-9]}, $merge_line);
			}
		}
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
#	my $total_chr = config::chromosome_number($genome);
	my $s = $_[0];
	my $i = $_[1];
	my $alleles = 1;
	if($homo == 0) { $alleles = 2; }
	for(my $a = 0; $a < $alleles; $a++) {
		$table_name = $genome . "_mutations_" . $s . "_allele_" . ($a + 1);
		#Check if table already exists
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == 1) {
			print STDERR "Table " . $table_name . " already exists!\n";
			if($force == 1) {
				print STDERR "Overwrite table!\n";
				$dbh->do("DROP TABLE " . $table_name);
			} else {
				print STDERR "Next database insert\n";
				return 1;
			}
		}
		#Create table
		$dbh->{ Warn } = 0;
		$command = "CREATE TABLE " . $table_name . "(chr varchar($chr_name_length) not null, pos bigint not null, reference varchar($length_mut) not null, strain varchar($length_mut) not null, index SERIAL PRIMARY KEY)";
		$dbh->do($command);
		$command = "CREATE INDEX  ON " . $table_name . " (chr, pos)";
		$dbh->do($command);
		$command = "CREATE INDEX ON " . $table_name . " (chr)";
		$dbh->do($command);
		$command = "CREATE INDEX ON " . $table_name . " (pos)";
		$dbh->do($command);
		$dbh->{ Warn } = 1;
		print STDERR "\n\nInserting mutations into database for $s for allele " . ($a + 1) . "\n";
		print STDERR "This may take a while\n";
		$dbh->{ AutoCommit } = 0;
		$sth = $dbh->prepare("INSERT INTO $table_name (chr, pos, reference, strain) VALUES (?, ?, ?, ?)");
		foreach my $l (@{$lines->[$a]->[$i - 9]}) {
			my @com_line = split(',', $l);
			$sth->execute(@com_line);
		}
		$sth->finish();
		$dbh->commit();
		$dbh->{ AutoCommit } = 1;
	#	`bash database_insert.sh $table_name tmp_merged_$h[$i].txt`;
		#Add here something from config for different genomes!!!!
	#	my $c;
		my @vector;
		my $save = [];
		open OUT, ">aaaaa.txt";
		foreach my $c (keys %{$total_chr}) {
			$query = 0;
			$ref_seq = 0;
			$command = "SELECT * FROM " . $genome . "_mutations_" . $s . "_allele_" . ($a + 1) . " WHERE length(reference) != length(strain) AND chr = \'$c\' ORDER  by pos";
			$sth = $dbh->prepare($command);
			$sth->execute();
			while (my $ref = $sth->fetchrow_hashref()) {
				print OUT $ref->{'pos'} . "\t" . $ref->{'reference'} . "\t" . $ref->{'strain'} . "\n";
				$old_pos = $ref->{'pos'};
				$old_seq = $ref->{'reference'};
				#iterate over all positions that are the same
				while($ref_seq < $ref->{'pos'}) {
					$ref_seq++;
					$query++;
					
				}
				#Reference longer than mutated strain
				if(length($ref->{'reference'}) > length($ref->{'strain'})) {
					print OUT "ref longer\n";
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
					print OUT "ref shorted\n";
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
				print OUT $ref_seq . "\t" . $query . "\t" . ($query - $ref_seq) . "\t" . $c . "\n";
			}
			$sth->finish();
		}

		#Create table
		$table_name = "offset_" . $genome . "_" . $s . "_allele_" . ($a + 1); 
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'" . $table_name . "\'") || die $DBI::errstr;
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == 1) {
			print STDERR "Table " . $table_name . " already exists!\n";
			if($force == 1) {
				print STDERR "Overwrite table!\n";
				$dbh->do("DROP TABLE " . $table_name);
			} else {
				print STDERR "Next database insert!\n";
				return 1;
			}
		}

		$dbh->{ Warn } = 0;
		$command = "CREATE TABLE " . $table_name . " (ref_pos bigint not null, strain_pos bigint not null, shift bigint not null, chr varchar($chr_name_length) not null, PRIMARY KEY (ref_pos, chr))";
		$dbh->do($command);
		$command = "CREATE UNIQUE INDEX ON " . $table_name . " (strain_pos, chr)";
		$dbh->do($command);
		$dbh->{ Warn } = 1;
		$dbh->{ AutoCommit } = 0;
		$sth = $dbh->prepare("INSERT INTO " . $table_name . " (ref_pos, strain_pos, shift, chr) VALUES (?, ?, ?, ?)");
		for(my $i = 0; $i < @vector; $i++) {
			for(my $j = 0; $j < @{ $vector[$i] }; $j++) {
				print OUT $vector[$i][$j] . "\t";
			}
			print OUT "\n";
			$sth->execute(@{ $vector[$i] });
		}
		$sth->finish();
		$dbh->commit();
		$dbh->{ AutoCommit } = 1;
	}
	print LOG "strain_" . $genome . "_" . $s . "\n";
	$dbh->disconnect();
}



sub merge{
        my $line = $_[0];
        if($line eq "") {
                return;
        }
        chomp $line;
        @current = split('\t', $line);
        #Save reference and position
        if(substr($current[0], 0, 3) eq "chr") {
                $current[0] = substr($current[0], 3);
        }
        my $pos = $current[1];
        my $ref = $current[3];
        my $merge_line;
        my @last;
        my @allele;
        my @all;
#       $current[4] =~ s/^\./$current[3]/g;
        $current[4] =~ s/\.//g;
        my @variants = split(',', $current[4]);

        unshift @variants, $current[3];

        if($filter == 1 && $current[6] ne "PASS") {
                return 1;
        }


        my $var;
        for(my $i = 9; $i < @current; $i++) {
                @all = split('/', ((split(':', $current[$i]))[0]));
                if($all[0] eq ".") {
                        $all[0] = 0;
                }
                if(@all == 1) {
                        $all[1] = 0;
                } else {
                        if($all[1] eq ".") {
                                $all[1] = 0;
                        }
                }
                if($homo == 1) {
                        shift @all;
                }
                @merge_line = ();
                for(my $h = 0; $h < @all; $h++) {
                        $var = $variants[$all[$h]];
                        $ref = $current[3];
                        $pos = $current[1];
                        #Check if reference and mutations are different
                        if($ref ne $var) {
                                if(length($ref) > 1 && length($var) > 1) {
                                        #Remove the first bases if they are the same
                                        while(substr($ref, -1) eq substr($var, -1) && length($ref) > 1 && length($var) > 1) {
                                                chop $ref;
                                                chop $var;
                                        }
                                        #Now start from the beginning
                                        my $max = (length($ref) > length($var)) ? length($var) : length($ref);
                                        my $max = (length($ref) > length($var)) ? length($var) : length($ref);
                                        for(my $l = 0; $l < $max - 1; $l++) {
                                                if(length($ref) > 1 && length($var) > 1 && substr($ref, 0, 1) eq substr($var, 0, 1)) {
                                                        $pos++;
                                                        $ref = unpack "xA*", $ref;
                                                        $var = unpack "xA*", $var;
                                                }
                                        }
                                }
                                if(length($ref) > $length_mut) { $length_mut = length($ref); }
                                if(length($var) > $length_mut) { $length_mut = length($var); }

                                #Check if there was already a mutation before this one
                                @last = split(',', $lines->[$h]->[$i-9]->[-1]) if @{$lines->[$h]->[$i-9]} > 0;
                                if(@{$lines->[$h]->[$i-9]} < 1 || $last[0] ne $current[0]) {
                                        $merge_line = $current[0] . "," . $pos . "," . $ref . "," . $var;
                                        push(@{ $lines->[$h]->[$i-9] }, $merge_line);
                                        next;
                                }
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
                                if($last[1] == $pos) { next; }
                                if($current[1] < $last[1] + length($last[2])) {
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
                                        next;
                                }
                                #Deletion in the strain so that the new mutation can not exist = ignore second mutation
                                $merge_line = $current[0] . "," . $pos . "," . $ref . "," . $var;
                                push(@{ $lines->[$h]->[$i-9]}, $merge_line);
                        }
                }
        }
}
