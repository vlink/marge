#!/usr/bin/perl

use Getopt::Long;

$_ = 0 for my($help, $return_code);
$_ = "" for my($install_dir, $homer_path, $data_folder, $cmd);

sub printCMD {
	print STDERR "Script to configure MARGE installation\n";
	print STDERR "Options:\n";
	print STDERR "-install_dir <path to installation directory>\n";
	print STDERR "-homer_path <path to HOMER (only necessary when HOMER not part of \$PATH\n";
	print STDERR "-data_folder <path to folder used for storing data> (default in install_dir)\n";
	print STDERR "-h|--help: print this help\n";
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

if($install_dir eq "") {
	$install_dir = `pwd`;
	chomp $install_dir;
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
print STDERR "Installing MARGE into " . $install_dir. "\n";


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

print STDERR "Creating folder structure for MARGE\n";
$cmd = "mkdir -p $install_dir/config";
$return_code = system($cmd);
&failed($success, $install_dir . "/config");

$cmd = "mkdir -p $install_dir/bin";
$return_code = system($cmd);
&failed($success, $install_dir . "/bin");

$cmd = "mkdir -p $data_folder/motif_distribution_background";
$return_code = system($cmd);
&failed($success, $data_folder . "/motif_distribution_background");


print STDERR "\n\n";
print STDERR "Writing config file\n";
open CONFIG, ">$install_dir/config/config.txt" or die "Could not open $install_dir/config/config.txt\n";
print CONFIG "data_folder=" . $data_folder . "/\n";
print CONFIG "homer_path=" . $homer_path . "/\n";
print CONFIG "motif_background_path=" . $data_folder . "/motif_distribution_background/";
close CONFIG;


print STDERR "Downloading MARGE\n";
print STDERR "Wget...\n";

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

