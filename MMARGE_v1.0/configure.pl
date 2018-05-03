#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# Copyright Verena M. Link <vlink@ucsd.edu>
# 
# This file is part of MMARGE
#
# MMARGE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MMARGE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


$_ = 0 for my($help, $return_code);
$_ = "" for my($install_dir, $homer_path, $data_folder, $cmd);

sub printCMD {
	print STDERR "\n\n";
	print STDERR "Script to configure MMARGE installation\n";
	print STDERR "Options:\n";
	print STDERR "\t-install_dir <path to installation directory>\n";
	print STDERR "\t-homer_path <path to HOMER (only necessary when HOMER not part of \$PATH\n";
	print STDERR "\t-data_folder <path to folder used for storing data> (default in install_dir)\n\n";
	print STDERR "\t-h|--help: print this help\n\n\n";
	exit;
}

GetOptions(     "install_dir=s" => \$install_dir,
                "homer_path=s" => \$homer_path,
                "data_folder=s" => \$data_folder,
                "h" => \$help,
                "help" => \$help)
        or die(&printCMD());

if($help == 1) {
	&printCMD();
}	

my $working_dir = `pwd`;
chomp $working_dir;

if($install_dir eq "") {
	$install_dir = $working_dir;
}

if($homer_path eq "") {
	$homer_path = `which homer`;
	chomp $homer_path;
	my @tmp = split("bin", $homer_path);
	$homer_path = $tmp[0];
}

if($data_folder eq "") {
	$data_folder = $install_dir . "/data/";
}
print $install_dir . "\n";
print $homer_path . "\n";
print $data_folder . "\n";

if(!-e $install_dir) {
	$cmd = "mkdir -p $install_dir";
	$return_code = system($cmd);
	&failed($return_code, $install_dir);
}
print STDERR "Installing MMARGE into " . $install_dir. "\n";


if(!-e $data_folder) {
	$cmd = "mkdir -p $data_folder";
	$return_code = system($cmd);
	&failed($return_code, $data_folder);
}
print STDERR "\n";

if($homer_path eq "") {
	print STDERR "Could not find HOMER\n";
	print STDERR "Please specify the path to HOMER with -homer_path\n";
	print STDERR "Abort\n";
	exit;
}

print STDERR "HOMER found in " . $homer_path . "\n\n";

&check_system_programs('tar');

print STDERR "Creating folder structure for MMARGE\n";
$cmd = "mkdir -p $install_dir/config";
$return_code = system($cmd);
&failed($success, $install_dir . "/config");

$cmd = "mkdir -p $install_dir/bin";
$return_code = system($cmd);
&failed($success, $install_dir . "/bin");

$cmd = "mkdir -p $data_folder/motif_distribution_background";
$return_code = system($cmd);
&failed($success, $data_folder . "/motif_distribution_background");

$cmd = "mkdir -p $data_folder/motif_bg";
$return_code = system($cmd);
&failed($success, $data_folder . "/motif_bg");

print STDERR "\n\n";
print STDERR "Writing config file\n";
open CONFIG, ">$install_dir/config/config.txt" or die "Could not open $install_dir/config/config.txt\n";
print CONFIG "marge_folder=" . $install_dir . "/\n";
print CONFIG "data_folder=" . $data_folder . "/\n";
print CONFIG "homer_path=" . $homer_path . "/\n";
print CONFIG "motif_background_path=" . $data_folder . "/motif_distribution_background/\n";
print CONFIG "motif_bg=" . $data_folder . "/motif_bg/\n";
print CONFIG "motif_file=" . $install_dir . "/config/homer_motifs.txt\n";
close CONFIG;

open CONFIG, ">$working_dir/config/config.txt" or die "Could not open $working_dir/config/config.txt\n";
print CONFIG "marge_folder=" . $install_dir . "/\n";
print CONFIG "data_folder=" . $data_folder . "/\n";
print CONFIG "homer_path=" . $homer_path . "/\n";
print CONFIG "motif_background_path=" . $data_folder . "/motif_distribution_background/\n";
print CONFIG "motif_bg=" . $data_folder . "/motif_bg/\n";
print CONFIG "motif_file=" . $install_dir . "/config/homer_motifs.txt\n";
close CONFIG;


sub failed {
	my $return_code = $_[0];
	my $folder = $_[1];
	if($return_code == 0) {
		print STDERR "" . $folder . " was successfully created\n";
	} else {
		print STDERR "Could not create folder " . $folder . "\n";
		print STDERR "Abort\n";
		exit;
	}
}

print STDERR "Check for perl modules\n";
&check_perl_modules();
print STDERR "All necessary perl modules were successfully installed\n";

print STDERR "\n\nCheck for R packages\n";
&check_R_package();
print STDERR "All R packages were successfully installed\n\n\n";
print STDERR "Please run install.pl to finish MMARGE installation!\n\n\n";

sub check_system_programs{
	my $program = $_[0];
	my $check = `which $_[0]`;
	@split = split('\s', $check);
	if($check eq "" || (@split > 1 && $split[1] eq "no")) {
		print STDERR "Program $_[0] not found!\n";
		print STDERR "Please install $_[0] and restart this script\n";
		exit;
	}
}

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

	my $load_SysCPU = eval{
		require Sys::CPU;
		Sys::CPU->import();
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
	if($load_SysCPU) {
		print STDERR "Module Sys::CPU is installed\n";
	} else {
		&install_module("Sys::CPU");
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
