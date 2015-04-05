#!/usr/bin/perl -w

#Script to interact with the database that saves all the mutation data of the mouse strains

use strict;
use DBI;
use Getopt::Long;
require '../general/config.pm';

my $con = config::database_connection();


print STDERR "Are you sure?\n";
print STDERR"\n\tDATABASE: " . $con->{'dbname'} . "\n";
print STDERR "Waiting for 10 seconds!\n";
for(my $i = 0; $i < 10; $i++) {
	print STDERR " . ";
	sleep(1);
}
print STDERR "\n";


#Create database connection
my $dbh;
my $sth;

$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
$sth = $dbh->prepare("select \'drop table \"\' || tablename || \'\" cascade;\' as del from pg_tables WHERE tablename LIKE \'%mm10%\'"); 
$sth->execute();

while(my $a = $sth->fetchrow_hashref()) {
	print $a->{'del'} . "\n";
	$dbh->do($a->{'del'});
}
