#!/usr/bin/perl

package config;
use strict;
use DBI;

my @split;

sub read_config{
	open CONFIG, "</home/vlink/mouse_strains/marge/config/config.txt";
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
#Do what every necessary in the future
sub database_connection{
	open CONFIG, "</home/vlink/mouse_strains/marge/config/config.txt";
	my $dbname;
	my $host;
	my $user;
	my $pw;
	my %connection;
	foreach my $line (<CONFIG>) {
		chomp $line;
		@split = split('=', $line);
		if($split[0] eq "database") {
			$dbname = $split[1];
			$connection{'dbname'} = $split[1];
		} elsif($split[0] eq "database_user") {
			$user = $split[1];
			$connection{'user'} = $split[1];;
		} elsif($split[0] eq "host") {
			$host = $split[1];
			$connection{'host'} = $split[1];
		} elsif($split[0] eq "database_pw") {
			$pw = $split[1];
			$connection{'pw'} = $split[1];
		}
	}	
	close CONFIG;
	return \%connection;
}

sub chromosome_number{
	my %chr;
	open CHR, "</home/vlink/mouse_strains/marge/config/chromosomes.txt";
	foreach my $line (<CHR>) {
		chomp $line;
		if($line eq "") { next; }
		@split = split('\t', $line);
		if($split[0] eq $_[0]) {
			my @chrom = split(",", $split[1]);
			foreach my $c (@chrom) {
				$chr{$c} = 1;
			}
			return \%chr;

		}
	}
	$chr{'NONE'} = 1;	
	return \%chr;
}

sub gene_format{
	open GENE, "</home/vlink/mouse_strains/marge/config/formats.txt" or die ("formats.txt does not exist!\n");;
	my %gf;
	my @s;
	foreach my $line (<GENE>) {
		chomp $line;
		if($line eq "") { next; }
		@split = split(':', $line);
		if($split[0] eq "$_[0]_genes") {
			@s = split('\t', $split[1]);
			for(my $i = 0; $i < @s; $i++) {
				$gf{$s[$i]} = $i;
			}
		}	
	}
	close GENE;
	return(\%gf);
}

sub name_format{
	open NAME, "</home/vlink/mouse_strains/marge/config/formats.txt" or die("formats.txt does not exist!\n");
	my %nf;
	my $name;
	my $run = 0;
	my @s;
	my $genome = $_[0];
	foreach my $line (<NAME>) {
		chomp $line;
		if($line eq "") { next; }
		@split = split(':', $line);
		if($split[0] eq $genome . "_name") {
			@s = split('\t', $split[1]);
			for(my $i = 0; $i < @s; $i++) {
				if($s[$i] eq "name") {
					$name = $i;
				} 
				if($s[$i] eq "transcript") {
					$nf{$i} = 1;
				}
			}
		}
	}
	close NAME;
	return($name, \%nf);
}	

sub promoter_format{
	open NAME, "</home/vlink/mouse_strains/marge/config/formats.txt" or die("formats.txt does not exist!\n");
	my %nf;
	my $split;
	my $ann;
	my @split;
	my @s;
	foreach my $line (<NAME>) {
		chomp $line;
		if($line eq "") { next; }
		@split = split(':', $line);
		if($split[0] eq "promoter_annotation") {
			@s = split('\t', $split[1]);
			for(my $i = 0; $i < @s; $i++) {
				@split = split(';', $s[$i]);
				if(@split > 1) {
					$nf{$split[0]} = $i;
					if($split[0] eq "name") {
						$split = $split[1] . "\t" . $split[2];
					} elsif($split[0] eq "annotation") {
						$ann = $split[1];
					}
				} else {
					$nf{$s[$i]} = $i;
				}
			}				
		}
	}
	close NAME;
	return($split, $ann, \%nf);
}
sub read_format{
        @split = split('\t', $_[0]);
        my $sign = $split[0];
        shift @split;
        my $chr_pos = "";
        my @vector = @split;
        for(my $i = 0; $i < @split; $i++) {
                if($split[$i] eq "chr") {
                        $chr_pos = $i;
                        last;
                }
        }
        return ($sign, $chr_pos, \@vector);
}

sub read_header{
	@split = split(',', $_[0]);
	for(my $i = 0; $i < @split; $i++) {
		$split[$i] =~ s/^ //g;
	}
	return(\@split);
}

sub get_strains{
	my %strains;
	my @split;
	my $organism = $_[0];

	my $con = &database_connection;
	my $dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});
	my $sth = $dbh->prepare("SELECT relname FROM pg_class c, pg_user u
		WHERE c.relowner = u.usesysid AND c.relkind = 'r' AND c.relname !~ '^pg_'
		AND c.relname LIKE \'" . $organism . "_mutations%allele_1\'");
	$sth->execute();
	while(my $strain = $sth->fetchrow_hashref()) {
		@split = split('_', $strain->{'relname'});
		$strains{$split[2]} = 1;
	}
	return \%strains;
}

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
	print STDERR "\n\nExiting...\n\n";
	exit;
}
1;
