#!/usr/bin/perl
use strict;

&check_perl_modules();
&check_R_package();

sub check_perl_modules {
	my $load_threads = eval{
		require threads;
		threads->import();
		1;
	};

	my $load_storable = eval{
		require Storable;
		Storable->import();
		1;
	};

	my $load_SetIntervalTree = eval{
		require Set::IntervalTree;
		Set::IntervalTree->import();
		1;	
	};

	my $load_ListUtil = eval{
		require List::Util;
		List::Util->import();
		1;
	};

	my $load_GetoptLong = eval{
		require Getopt::Long;
		Getopt::Long->import();
		1;
	};

	my $load_StatisticsBasic = eval{
		require Statistics::Basic;
		Statistics::Basic->import();
		1;
	};

	my $load_FileCache = eval{
		require FileCache;
		FileCache->import();
		1;
	};

	if($load_threads) {
		print STDERR "Module threads is installed\n";
	} else {
		&install_module("threads");
	}
	if($load_storable) {
		print STDERR "Module Storable is installed\n";
	} else {
		&install_module("Storable");
	}	
	if($load_SetIntervalTree) {
		print STDERR "Module Set::IntervalTree is installed\n";
	} else {
		&install_module("Set::IntervalTree");
	}
	if($load_ListUtil) {
		print STDERR "Module List::Util is installed\n";
	} else {
		&install_module("List::Util");
	}
	if($load_GetoptLong) {
		print STDERR "Module Getopt::Long is installed\n";
	} else {
		&install_module("Getopt::Long");
	}
	if($load_StatisticsBasic) {
		print STDERR "Module Statistics::Basic is installed\n";
	} else {
		&install_module("Statistics::Basic");
	}
	if($load_FileCache) {
		print STDERR "Module FileCache is installed\n";
	} else {
		&install_module("FileCache");
	}

} 
sub install_module{
	my $module = $_[0];
	if($> != 0) {
		print STDERR "User has no root privileges!\n";
		print STDERR "Installing perl modules using cpan is not possible!\n";
		print STDERR "Please execute this script with root priviliges or install perl module $module before restarting!\n";
		exit;
	} else {
		print STDERR "User has root priviliges!\n";
		print STDERR "Installing module $module with cpan\n";
		`perl -MCPAN -e 'install $module'`
	}
	my $install_success = eval{
		require $module;
		($module)->import();
		1;
	};
	if(!$install_success) {
		print STDERR "Could not automatically install Module $module\n";
		print STDERR "Please install module manually and restart script afterwards!\n";
		exit;
	} else {
		print STDERR "Moduel $module successfully installed!\n";
	}
}

sub check_R_package {
	#open R file
	open R, ">test.R";
	print R "list.of.packages <- c(\"seqinr\")\n";
	print R "new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,\"Package\"])]\n";
	print R "if(length(new.packages)) install.packages(new.packages)\n";
	close R;
	my $success = `RScript test.R`;
	if($success == 0) {
		print STDERR "Package seqinr is available on this computer!\n";
	} else {
		print STDERR "R Package seqinr needs to be installed!\n";
	}
	`rm test.R`;
}
