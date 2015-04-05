#!/usr/bin/perl -w

use strict;
use DBI;

require '../general/config.pm';

if(@ARGV < 1) {
	print STDERR "specify genome!\n";
	exit;
}

my $genome = $ARGV[0];

my $strains = config::get_strains($genome);
my $con = config::database_connection();

my $dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
my $sth;
my $num;
my $old_pos = 0;

my $total = config::chromosome_number($genome);
my $cur_pos = 0;

my %count;

my $happens = 0;
my $show_muts;

foreach my $s (keys %{$strains}) {
	if($s eq "test") { next; }
	%count = ();
	for(my $j = 1; $j < 50; $j++) {
		$count{$j} = 0;
	}
	#Check if there are two alleles
	$sth = $dbh->prepare("SELECT count(*) FROM pg_class c, pg_user u where c.relowner = u.usesysid AND c.relkind = \'r\' AND c.relname !~ '^pg_' AND c.relname LIKE \'" . $genome . "_mutations_" . $s . "%\'");
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();

	foreach my $c (keys %{$total}) {
		$old_pos = 0;
		if($num->{'count'} == 1) {
			$sth = $dbh->prepare("SELECT pos FROM " . $genome . "_mutations_" . $s . "_allele_1 AS one WHERE one.chr = \'$c\' ORDER BY index");
		} else {
			$sth = $dbh->prepare("SELECT " . $genome . "_mutations_" . $s . "_allele_1.pos AS one, " . $genome . "_mutations_" . $s . "_allele_2.pos AS two FROM " . $genome . "_mutations_" . $s . "_allele_1 FULL JOIN " . $genome . "_mutations_" . $s . "_allele_2 ON " . $genome . "_mutations_" . $s . "_allele_1.pos = " . $genome . "_mutations_" . $s. "_allele_2.pos WHERE " . $genome . "_mutations_" . $s . "_allele_1.chr = \'$c\' AND " . $genome . "_mutations_" . $s . "_allele_2.chr = \'$c\'");
		}
		$sth->execute();
		$happens = 0;
		$show_muts = $s . "\t" . $c . "\t";
		while(my $mut = $sth->fetchrow_hashref()) {
			$cur_pos = $mut->{'one'};
			if($cur_pos eq "") {
				$cur_pos = $mut->{'two'};
			}
			if(abs($cur_pos - $old_pos) < 50) {
			#	$show_muts .= $mut->{'one'} . " (curr: $cur_pos old: $old_pos diff: " . ($cur_pos- $old_pos) . ")\t";;
				$show_muts .= $mut->{'one'} . "\t";;
				$happens++;
			} else {
				if($happens > 1) {
					print $show_muts . "\n";
				}
				$show_muts = $s . "\t" . $c . "\t";
				$count{$happens}++;
				$happens = 0;
				$old_pos = $cur_pos;

			}
		}
	}	
#	foreach my $k (sort {$a <=> $b} keys %count) {
#		print $k . "\t" . $count{$k} . "\n";
#	}
}
