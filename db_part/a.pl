#!/usr/bin/perl -w

use strict;
require '../general/config.pm';
use DBI;
#Define variables

my $con = config::database_connection();
my $dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
my $sth;
my $strains;

my $command = "SELECT column_name FROM information_schema.columns WHERE table_name LIKE \'mouse_mutations_%\'";
$sth = $dbh->prepare($command);
$sth->execute();
while(my $s = $sth->fetchrow_arrayref()) {
	$strains = config::get_strains("mouse");
}
$sth->finish();


$sth = $dbh->prepare("SELECT DISTINCT transcript FROM mouse_genes");
$sth->execute();

my $count = 0;

while(my $trans = $sth->fetchrow_hashref()) {
	$command = "perl get_gene.pl -organism mouse -method gene -transcript " . $trans->{'transcript'} . " -strains ";
	foreach my $key (keys %{$strains}) {
		$command .= $key . ", ";
	}
	chop $command;
	print STDERR $command . "\n";
	`$command`;
}
