#!/usr/bin/perl -w

use strict;
use DBI;
require '../general/config.pm';

my $dbh;
my $sth;

my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

chomp $ARGV[0];
print "collect data\n";
$sth = $dbh->prepare("SELECT chromosome, position, reference, $ARGV[0] from mouse_mutations where length(reference) != length($ARGV[0])");
$sth->execute();

my $a;

print "find genes\n";
while(my $mut = $sth->fetchrow_hashref()) {
	$a = $dbh->prepare("SELECT * FROM mouse_genes WHERE chr = \'$mut->{'chromosome'}\' AND start < $mut->{'position'} AND stop > $mut->{'position'}");
	$a->execute();
	while(my $b = $a->fetchrow_hashref()) {
		print $b->{'transcript'} . "\n";
	}
}
