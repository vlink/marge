#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

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


$_ = () for my(@split, %config, @name, @change, @modules);
$_ = "" for my($command, $line);

print STDERR "Adjust paths in perl files\n";
open CONFIG, "<config/config.txt";
foreach my $line (<CONFIG>) {
	chomp $line;
	@split = split('=', $line);
	$config{$split[0]} = $split[1];
}
#Copy motif file
$command = "cp config/homer_motifs.txt " . $config{'marge_folder'} . "/config/homer.motifs 2> /dev/null";
`$command`;


close CONFIG;
my $working_dir = `pwd`;
chomp $working_dir;

@modules = `ls $working_dir/bin/*pm`;
&change_path(\@modules);

@modules = `ls $working_dir/bin/*pl`;
&change_path(\@modules);


sub change_path {
	my @mod = @{$_[0]};
	my $count = 1;
	foreach my $m (@mod) {
		chomp $m;
		print $m . "\n";
		@name = split('/', $m);
		open OUT, ">.change.pl";
		open FH, "<$m";
		$count = 1;
		if($name[-1] eq "findMotifsGenome_strains.pl") {
			while(<FH>) {
				$line = $_;
				chomp $line;
				if($count == 1) {
					print OUT "#!/usr/bin/env perl\n";
				} elsif($count == 2) {
					print OUT "BEGIN {push \@INC, \'" . $config{'marge_folder'} . "/bin\'}\n";
					print OUT "BEGIN {push \@INC, \'" . $config{'marge_folder'} . "/config\'}\n";
				} elsif($count == 3) {
					print OUT "use warnings;\n";
				} elsif($count == 4) {
					print OUT "use lib \"" . $config{'homer_path'} . "./bin\";\n";
				} elsif($count == 5) {
					print OUT "my \$homeDir = \"" . $config{'homer_path'} . "./\";\n";
				} else {
					if($line =~ m/PATH_TO_HOMER_MOTIF_FILE/g) {
						@change = split("PATH_TO_HOMER_MOTIF_FILE", $line);
						print OUT $change[0] . $config{'marge_folder'} . "/config/homer.motifs" . $change[1] . "\n"; 
					} elsif($line =~ m/PATH_TO_GET_SEQ_FILES/g) {
						@change = split("PATH_TO_GET_SEQ_FILES", $line);
						print OUT $change[0] . $config{'marge_folder'} . "/bin/extract_seq_from_peakfiles.pl" . $change[1] . "\n"; 
					} else {
						print OUT $line . "\n";
					}
				}
				$count++;
			}
		} else {
			while(<FH>) {
				$line = $_;
				chomp $line;
				if($count == 1) {
					print OUT "#!/usr/bin/env perl\n";
				} elsif($count == 2) {
					print OUT "BEGIN {push \@INC, \'" . $config{'marge_folder'} . "/bin\'}\n";
				} else {
					if($line =~ m/CHANGE_TO_CONFIG_FILE/g) {
						@change = split("CHANGE_TO_CONFIG_FILE", $line);
						print OUT $change[0] . $config{'marge_folder'} . "/config/config.txt" . $change[1] . "\n"; 
					} elsif($line =~ m/CHANGE_TO_WORKDIR/g) {
						@change = split("CHANGE_TO_WORKDIR", $line);
						print OUT $change[0] . $working_dir . "/bin/" . $change[1] . "\n";
					} else {
						print OUT $line . "\n";
					}
				}
				$count++;
			}
		}	
		close OUT;
		close FH;
		$command = "mv .change.pl " . $config{'marge_folder'} . "/bin/" . $name[-1];
		`$command`;
		
	}
}

`chmod -R a+x $config{'marge_folder'}/bin`;

print STDERR "\n\n\nMMARGE is installed successfully\n";
print STDERR "To add MMARGE to your path run:\n";
print STDERR "\n\texport PATH=\$PATH:" . $config{'marge_folder'} . "/bin\n";
print STDERR "\nTo add MMARGE path permanantely to your PATH environment run:\n";
print STDERR "\n\techo \"export PATH=\$PATH:" . $config{'marge_folder'} . "/bin && export PATH\" >> ~/.bash_profile\n\n\n";
