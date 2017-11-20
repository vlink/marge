#!/usr/bin/perl

package config;
use strict;

my @split;

#Read in config file and return different varibales in a hash
sub read_config{
	open CONFIG, "</home/vlink/mouse_strains/marge/general/config.txt";
	my %parameters;
	my @split;
	foreach my $line (<CONFIG>) {
		chomp $line;
		if($line eq "" || substr($line, 0, 1) eq "#") {
			next;
		}
		@split = split('=', $line);
		$parameters{$split[0]} = $split[1];
	}
	close CONFIG;
	return \%parameters;
}

#Check that parameters specified in mand are in the programm call
sub check_parameters{
	my $mand = $_[0];
	my $argv = $_[1];
	my %missing;
	foreach my $key (keys %{$mand}) {
		if(!exists $argv->{$key}) {
			$missing{$key} = 1;
		}
	}

	if(keys %missing == 0) {
		return 0;
	}

	print STDERR "\n\nPlease specify following parameter(s) or use -h or --help to print help\n";
	foreach my $key (keys %missing) {
		print STDERR "\t$key\n";
	}
	print STDERR "\n\n";
	exit;
}
1;
