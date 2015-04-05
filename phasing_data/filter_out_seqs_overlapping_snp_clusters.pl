#!/usr/bin/perl -w

use strict;
#require '/data/home/vlink/mouse_strains/software/general/config.pm';
require '/data/home/vlink/mouse_strains/software/db_part/database_interaction.pm';
use Data::Dumper;
use DBI;
use Getopt::Long;

sub printCMD {
        print STDERR "Usage:\n";
        print STDERR "\t-genome: Genome\n";
	print STDERR "\t-strain: Strain\n";
	print STDERR "\t-sam: sam file\n";
	print STDERR "\t-mapped: Genome you used to map (either the reference or the strains genome)\n";
	print STDERR "\t-read: Length of the reads you get back from the sequences (default: 50)\n";
        exit;
}

if(@ARGV < 1) {
        &printCMD();
}


my @matrix;
&dynamic_progamming();
exit;
$_ = "" for my ($genome, $gene, $transcript, $chr, $method, $strand, $genome2, $mapped);
$_ = 0 for my ($start, $end, $html, $happens);
my ($dbh, $sth, @strain, @strains2, $sam, $strain);
my $homo = 1;
my $read_length = 50;

my %mandatory = ('-genome' => 1, '-strain' => 1, '-sam' => 1);
my %convert = map { $_ => 1 } @ARGV;

config::check_parameters(\%mandatory, \%convert);

GetOptions(     "genome=s" => \$genome,
                "strain=s" => \@strain,
                "sam=s" => \$sam,
		"read=s" => \$read_length,
                "mapped=s" => \$mapped)
        or die("Error in command line options!\n");

if($mapped eq "") {
	$mapped = $genome;
}
$strain = $strain[0];

database_interaction::set_global_variables(\@strain, $genome, $homo, $html, "genomic");
my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});



#Frist step get cluster of snps if it does not already exist
my $cluster_file = "/data/home/vlink/mouse_strains/software/config/snp_cluster/" . $genome . "_" . $strain . "_" . $read_length . "_mapped_" . $mapped . ".cluster";
print $cluster_file . "\n";
if(-e $cluster_file) {
	my $size = (split('\s+', `ls -lah $cluster_file`))[4];
	if($size eq "0") {
		print STDERR "Find variants clusters!\n";
		&create_cluster_file($cluster_file);
	}
} else {
	print STDERR "Find variants clusters!\n";
	&create_cluster_file($cluster_file);
}

my @split;
my %save_cluster;
my $run = 0;

#Save clusters
print STDERR "Save clusters!\n";
open FH, "<$cluster_file";
foreach my $line (<FH>) {
	chomp $line;
	@split = split('\t', $line);
	for(my $i = 2; $i < @split; $i++) {
		$save_cluster{$split[1]}[$split[$i]] = $run;
	}	
	$run++;
}
close FH;

print STDERR "Done\n";

print STDERR "Read in sam file!\n";
open FH, "<$sam";
my (%save_cands, %save_num, %info, $num);

#open OUT, ">hg19_candidates_jun.txt";
foreach my $line (<FH>) {
	chomp $line;
	if(substr($line, 0, 1) eq "@") {
		next;
	}
	@split = split('\t', $line);
	$chr = substr($split[2], 3);
	#Save cluster information
	
	for(my $i = 0; $i < $read_length; $i++) {
		if(defined $save_cluster{$chr}[$split[3] + $i]) {
			if(exists $save_num{$save_cluster{$chr}[$split[3] + $i]}) {
				$num = $save_num{$save_cluster{$chr}[$split[3] + $i]} + 1;
				$save_num{$save_cluster{$chr}[$split[3] + $i]}++;
				if($split[3] < $info{$save_cluster{$chr}[$split[3] + $i]}->{'min'}) {$info{$save_cluster{$chr}[$split[3] + $i]}->{'min'} = $split[3]; }
				if($split[3] > $info{$save_cluster{$chr}[$split[3] + $i]}->{'max'}) {$info{$save_cluster{$chr}[$split[3] + $i]}->{'max'} = $split[3]; }
			} else {
				$num = 0;
				$save_num{$save_cluster{$chr}[$split[3] + $i]} = 0;
				$info{$save_cluster{$chr}[$split[3] + $i]}->{'min'} = $split[3];
				$info{$save_cluster{$chr}[$split[3] + $i]}->{'max'} = $split[3];
				$info{$save_cluster{$chr}[$split[3] + $i]}->{'strand'} = $split[1];
				$info{$save_cluster{$chr}[$split[3] + $i]}->{'chr'} = substr($split[2], 3);
			}
		#	print OUT $line . "\tCLUSTER: " . $chr . "\t" . ($split[3] + $i) . "\n";;
			$save_cands{$save_cluster{$chr}[$split[3] + $i]}[$num] = $line;
		}
	}		
}
close FH;
#close OUT;

#Create alignments with sam file
my $seq = "";
my $start_ref;
my $start_ali;
my $diff = 0;
my $first = 0;
my($shift_before, $shift_after, $result, $start_index, $stop_index, $real_start, $real_stop);
#Prepare for dynamic programming
foreach my $key (sort {$a <=> $b } keys %save_cands) {
	if($strain ne $mapped) {
		$real_start = $info{$key}->{'min'} - 1;
		$real_stop = $info{$key}->{'max'} + $read_length - 1;
	} else {
		$real_start = $info{$key}->{'min'}- 1;
		$real_stop = $info{$key}->{'max'} + $read_length - 1;

		$sth = $dbh->prepare("SELECT * FROM offset_" . $genome . "_" . $strain . "_allele_1 WHERE strain_pos < " . $real_start . " AND chr = \'" . $info{$key}->{'chr'}. "\' ORDER BY strain_pos DESC LIMIT 1");
		$sth->execute();
		$result = $sth->fetchrow_hashref();
		$real_start = $real_start - $result->{'shift'};
		$sth->finish();
			
		$sth = $dbh->prepare("SELECT * FROM offset_" . $genome . "_" . $strain . "_allele_1 WHERE strain_pos < " . $real_stop . " AND chr = \'" . $info{$key}->{'chr'}. "\' ORDER BY strain_pos DESC LIMIT 1");
		$sth->execute();
		$result = $sth->fetchrow_hashref();
		$real_stop = $real_stop - $result->{'shift'};
		$sth->finish();
	}
	$sth = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strain . "_allele_1 where pos >= " . $real_start . "AND pos <= " . $real_stop . " AND chr ='" . $info{$key}->{'chr'} . "\'");
	$sth->execute();
	while($result = $sth->fetchrow_hashref()) {
		if($first == 0) {
			$start_index = $result->{'index'} - 1;
			$first++;
		}
		$stop_index = $result->{'index'};
		if(length($result->{'reference'}) > length($result->{'strain'})) {
			$real_stop += length($result->{'reference'}) - length($result->{'strain'});
		}
	}
	$sth->finish();
	$stop_index++;
	#Go from here
	$sth = $dbh->prepare("SELECT * from " . $genome . "_mutations_" . $strain. "_allele_1 where index = " . $start_index);
	$sth->execute();
	$result = $sth->fetchrow_hashref();
	
	if($result->{'pos'} + length($result->{'reference'}) >= $real_start) {
			my $a = $result->{'pos'} + length($result->{'reference'}) - $info{$key}->{'min'};
			$real_start = $info{$key}->{'min'} - length($result->{'reference'});
	}
	$sth->finish();

	$sth = $dbh->prepare("SELECT * from " . $genome . "_mutations_" . $strain . "_allele_1 where index = " . $stop_index);
	$sth->execute();
	$result = $sth->fetchrow_hashref();
	if($result->{'pos'} < ($info{$key}->{'max'}  + $read_length - 1)) {
			$real_stop += $info{$key}->{'min'} - $result->{'pos'};
	}
	$sth->finish();
	$seq = "";
	$seq = database_interaction::get_genomic_seq($real_start, $real_stop, $info{$key}->{'chr'}, 0, 0, 0, 0); 
	$seq =~ s/-//g;
	#Initialize matrix
	@matrix = ();
	#create matrix
	my $i = 0;
	foreach my $candidates ($save_cands{$key}) {
		foreach my $cand_seqs (@{$candidates}) {
			@split = split("\t", $cand_seqs);
			$start_ref = $split[3] - $info{$key}->{'min'};
			$start_ali = 0;
			for(my $j = 0; $j < $start_ref; $j++) {
				$matrix[$i][$j] = "-";
			}
			for(my $j = 0; $j < $start_ref; $j++) {
				$matrix[$i][$j] = "-";
			}
			for(my $j = $start_ref; $j < $start_ref + $read_length; $j++) {
				if(substr($seq, $j, 1) ne substr($split[9], $start_ali, 1)) {
					$matrix[$i][$j] = 1;
				} else {
					$matrix[$i][$j] = 0;
				}
				$start_ali++;
			}
			for(my $j = $start_ref + $read_length; $j < $info{$key}->{'max'} - $info{$key}->{'min'} + $read_length; $j++) {
				$matrix[$i][$j] = "-";
			}
			$i++;			
		}
	}
	#Check if all position in the matrix agree
	$diff = 0;
	for(my $j = 0; $j < @{$matrix[0]}; $j++) {
		$first = $matrix[0][$j];
		for(my $i = 0; $i < @matrix; $i++) {
			if($matrix[$i][$j] eq "-") { next; }
			if($first eq "-") { $first = $matrix[$i][$j]; }
			if($matrix[$i][$j] != $first) {
				$diff = 1;
			}
		}
	}
	if($diff == 1) {
		print "Time to check that out in detail!\n";
		for(my $i = 0; $i < @matrix; $i++) {
			for(my $j = 0; $j < @{$matrix[$i]}; $j++) {
				print $matrix[$i][$j];
			}
			print "\n";
		}
	} else {
		print "only one hapotype!\n";
	}
	$first = 0;
	@matrix = ();
}

sub create_cluster_file{
	my $file = $_[0];
	open OUT, ">$file";
	my $old_pos = 0;
	my $total = config::chromosome_number($genome);

	$_ = 0 for my ($cur_pos, $happes, $show_muts, $shift, $second, $shift_vector);
	my %count = ();
	%count = ();
	for(my $j = 1; $j < $read_length; $j++) {
		$count{$j} = 0;
	}
	#Check if there are two alleles
	$sth = $dbh->prepare("SELECT count(*) FROM pg_class c, pg_user u where c.relowner = u.usesysid AND c.relkind = \'r\' AND c.relname !~ '^pg_' AND c.relname LIKE \'" . $genome . "_mutations_" . $strain . "%\'");
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();

	foreach my $c (keys %{$total}) {
		$old_pos = 0;
		if($num->{'count'} == 1) {
			$sth = $dbh->prepare("SELECT pos FROM " . $genome . "_mutations_" . $strain . "_allele_1 AS one WHERE one.chr = \'$c\' ORDER BY index");
		} else {
			$sth = $dbh->prepare("SELECT " . $genome . "_mutations_" . $strain . "_allele_1.pos AS one, " . $genome . "_mutations_" . $strain . "_allele_2.pos AS two FROM " . $genome . "_mutations_" . $strain . "_allele_1 FULL JOIN " . $genome . "_mutations_" . $strain . "_allele_2 ON " . $genome . "_mutations_" . $strain . "_allele_1.pos = " . $genome . "_mutations_" . $strain . "_allele_2.pos WHERE " . $genome . "_mutations_" . $strain . "_allele_1.chr = \'$c\' AND " . $genome . "_mutations_" . $strain . "_allele_2.chr = \'$c\'");
		}
		$sth->execute();
		$happens = 0;
		$show_muts = $strain . "\t" . $c . "\t";
		while(my $mut = $sth->fetchrow_hashref()) {
			$second = 0;
			$cur_pos = $mut->{'one'};
			if($cur_pos eq "") {
				$cur_pos = $mut->{'two'};
				$second = 1;
			}
			if(abs($cur_pos - $old_pos) < $read_length) {
				#Check if it was mapped against the strains genome - in this case we need to shift coordinates!
				if($mapped eq $strain) {
					if($second == 0) {
						$shift = $dbh->prepare("SELECT * from offset_" . $genome . "_" . $strain . "_allele_1 WHERE chr = \'" . $c . "\' AND ref_pos < " . $cur_pos . " ORDER BY ref_pos DESC LIMIT 1");
						$shift->execute();
						$shift_vector = $shift->fetchrow_hashref();
						$shift->finish();
						if(!exists $shift_vector->{'shift'} ) { $shift_vector->{'shift'} = 0; }
						$show_muts .= ($mut->{'one'} + $shift_vector->{'shift'}) . "\t"; 
					} else {
						$shift = $dbh->prepare("SELECT * from offset_" . $genome . "_" . $strain . "_allele_2 WHERE chr = \'" . $c . "\' AND ref_pos < " . $cur_pos . " ORDER BY ref_pos DESC LIMIT 1");
						$shift->execute();
						$shift_vector = $shift->fetchrow_hashref();
						$shift->finish();
						if(!exists $shift_vector->{'shift'} ) { $shift_vector->{'shift'} = 0; }
						$show_muts .= ($mut->{'two'} + $shift_vector->{'shift'}) . "\t"; 
					}
				} else {
					$show_muts .= $mut->{'one'} . "\t";
				}
				$happens++;
			} else {
				if($happens > 1) {
					print OUT $show_muts . "\n";
				}
				$show_muts = $strain . "\t" . $c . "\t";
				$count{$happens}++;
				$happens = 0;
				$old_pos = $cur_pos;

			}
		}
	}
	close OUT;
}

sub dynamic_progamming{
	@matrix = ();
	$matrix[0][0] = 0;
	$matrix[0][1] = 0;
	$matrix[0][2] = '-';
	$matrix[0][3] = '-';
	$matrix[0][4] = '-';
	$matrix[1][0] = '-';
	$matrix[1][1] = '-';
	$matrix[1][2] = 1;
	$matrix[1][3] = 1;
	$matrix[1][4] = '-';
	$matrix[2][0] = 0;
	$matrix[2][1] = 0;
	$matrix[2][2] = 0;
	$matrix[2][3] = 0;
	$matrix[2][4] = '-';
	$matrix[3][0] = '-';
	$matrix[3][1] = '-';
	$matrix[3][2] = 1;
	$matrix[3][3] = 0;
	$matrix[3][4] = 1;

	my @length = ();
	$length[0] = 2;
	$length[1] = 2;
	$length[2] = 4;
	$length[3] = 3;

	my @trace_zero;
	my @trace_one;
	for(my $i = 0; $i < 6; $i++) {
		for(my $j = 0; $j < 6; $j++) {
			$trace_zero[$i][$j] = 0;
		}
	}

	my $max_length = 0;
	#Go through matrix and find start at this position
	my $start_point = 0;
	my @consider = ();

	for(my $i = 0; $i < @matrix; $i++) {
		for(my $j = 0; $j < @{$matrix[0]}; $j++) {	
			print "i: " . $i . "\t" . "j: " . $j . "\n";
			print "value: " . $matrix[$i][$j] . "\n";
			if($matrix[$i][$j] eq "-") { next; }
			if($matrix[$i][$j] == 1) {
				$trace_zero[$j][$i]++;
			}
			print "matrix entry " . $i . " " . ($j + 1) . " is same if value is 0, goes 1 up if value is 1\n";
			print "matrix entry" . ($i + 1) . " " . $j . " is same if value is 1, goes 1 up if value is 0\n";
		}
	#	print "trace one\n";
	#	print Dumper @trace_one;
	#	print "\n\n\n";
		print "trace zero\n";
		for(my $k = 0; $k < @trace_zero; $k++) {
			for(my $k2 = 0; $k2 < @trace_zero; $k2++) {
				print $trace_zero[$k][$k2] . "\t";
			}
			print "\n";
		}
		print "\n\n\n";
	}
}
