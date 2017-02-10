#!/usr/bin/perl
use strict;
use Data::Dumper;

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
$command = "cp config/homer_motifs.txt " . $config{'marge_folder'} . "/config";
`$command`;


close CONFIG;
my $working_dir = `pwd`;
chomp $working_dir;

@modules = `ls $working_dir/bin/*pm`;
&change_path(\@modules);

@modules = `ls $working_dir/bin/*pl`;
&change_path(\@modules);

my @update;
$update[0] = "update.pl";
&change_path(\@update);

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
				if($count == 1) {
					print OUT "#/usr/bin/env perl\n";
				} elsif($count == 2) {
					print OUT "BEGIN {push \@INC, \'" . $config{'marge_folder'} . "/bin\'}\n";
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
						print OUT $change[0] . $config{'marge_folder'} . "/bin/extract_seq_from_peakfiles.pl" . $change[1]; 
					} else {
						print OUT $line;
					}
				}
				$count++;
			}
		} else {
			while(<FH>) {
				$line = $_;
				if($count == 1) {
					print OUT "#/usr/bin/env perl\n";
				} elsif($count == 2) {
					print OUT "BEGIN {push \@INC, \'" . $config{'marge_folder'} . "/bin\'}\n";
				} else {
					if($line =~ m/CHANGE_TO_CONFIG_FILE/g) {
						@change = split("CHANGE_TO_CONFIG_FILE", $line);
						print OUT $change[0] . $config{'marge_folder'} . "/config/config.txt" . $change[1] . "\n"; 
					} else {
						print OUT $line;
					}
				}
				$count++;
			}
		}	
		close OUT;
		close FH;
		if($name[-1] eq "update.pl") {
			$command = "mv .change.pl " . $config{'marge_folder'} . "/" . $name[-1];
			`$command`;
		} else {
			$command = "mv .change.pl " . $config{'marge_folder'} . "/bin/" . $name[-1];
			`$command`;
		}
		
	}
}
