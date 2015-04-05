#!/usr/bin/perl -w

package database_interaction;
use strict;
use DBI;
use Data::Dumper;
require '../general/config.pm';

#Define variables
$_ = () for my (%comp, %nt_to_as, %mut_per_gene);
my ($dbh, $sth, $command, $trans, $id, $ref_index, $comp_seq, @c);

	my $lift = 0;
		my %lift;
my $con = config::database_connection();
$dbh = DBI->connect("DBI:Pg:dbname=$con->{'dbname'};host=$con->{'host'}",  $con->{'user'}, $con->{'pw'}, {'RaiseError' => 1});

#Initialize the hashes for the script
&initialize_data();
our ($strains, $genome, $homo, $allele, $html);

sub set_global_variables{
	$strains = $_[0]; 
	$genome = $_[1];
	$homo = $_[2];
	$html = $_[3];
	my $name = $_[4];
	if($homo == 1) {
		$allele = 1;
	} else {
		$allele = 2;
	}
	if($html == 1) {
		open OUT, ">$name.html";
	}
	return 0;
}

sub get_transcripts{
	my %transcript;
	$trans = $_[0];
	$id = $_[1];
	if($trans eq "") {
		$command = "SELECT DISTINCT transcript FROM " . $genome . "_genes WHERE id=\'$id\'";
		$sth = $dbh->prepare($command);
		$sth->execute();
		while(my $t = $sth->fetchrow_hashref()) {
			$transcript{$t->{'transcript'}} = 1;
		}
		$sth->finish();
	} else {
		$transcript{$trans} = 1;
	}
	return \%transcript;
}

sub get_gene{
        my $id = $_[0];
        my $trans = $_[1];
	my $print_gene = $_[2];
	my $align_gene = $_[3];
	my $print_protein = $_[4];
	my $align_protein = $_[5];
	#Start with the DB
	my $transcript  = &get_transcripts($trans, $id);
        use integer;
	$_ = 0 for my ($exon, $start, $stop, $new_start, $mut_start, $mut_stop, $index, $first, $mut_shift, $run, $align);
	$_ = "" for my ($output_seq, $strand, $current_exon, $seq, $mut, $chr, $as, $print_seq);	
	$_ = () for my (@all_seqs, $ref_seq, $diff, $end, $ass, $align_string, $mut_seq, $last_exon, $general_shift);
	foreach my $trans (keys %{$transcript}) {
		if($print_gene == 1 || $print_protein == 1) {
			if($html == 1) { 
				print OUT $trans . "<br>\n";
			} else {
				print $trans . "\n";
			}
		}
		#Fetch exons from DB
		my $exon_db = $dbh->prepare("SELECT max(exon) FROM " . $genome . "_genes where transcript=\'$trans\'");
		$exon_db->execute();
		$sth = $dbh->prepare("SELECT * FROM " . $genome . "_genes WHERE transcript=\'$trans\'");
		$sth->execute();
		$exon = 0;
		if(!defined ($last_exon = $exon_db->fetchrow_arrayref->[0])) {
			print "Gene not found in database!\n";
			return 1;
		}
		$exon_db->finish();
		#Generate the gene per strain
		while(my $ref = $sth->fetchrow_hashref()) {
			$strand = $ref->{'strand'};
			$chr = $ref->{'chr'};
			#Get the next exon
			if($ref->{'start'} eq $ref->{'stop'}) {
				print "no transcript\n";
				return 1;
			}
			if($exon == $ref->{'exon'} - 1) {
				$exon++;
				#First exon not part of cds - skip
				if($ref->{'start'} > $ref->{'stop'} || $new_start > $ref->{'stop'}) {
					$new_start = $ref->{'start'} if $new_start == 0;
					$start = ($new_start/100);
					$mut_start = $new_start;
					next;
				} else {
					$mut_start = $ref->{'start'};
					$start = ($ref->{'start'}/100);
					if($new_start == 0) {
						$mut_start = $ref->{'start'};
					}
				}
				$stop = ($ref->{'stop'}/100);
				$mut_stop = $ref->{'stop'};
				$new_start = 0;
				#Create reference sequence
				$current_exon = &get_seq($start, $stop, $chr, $mut_start, $mut_stop, $ref->{'exon'}, $last_exon);
				$output_seq .= $current_exon;
				for(my $index = 0; $index < @{$strains}; $index++) {
					if($strains->[$index] eq "reference") { next; }
					for(my $a = 0; $a < $allele; $a++) {
						my ($refs, $muts, $align) = &align_seqs($genome, $current_exon, $strains->[$index], ($a + 1), $mut_start, $mut_stop, $chr, $ref->{'start'}, $ref->{'stop'});
						$ref_seq->[$a]->[$index] .= $refs;
						$mut_seq->[$a]->[$index] .= $muts;
						$align_string->[$a]->[$index] .= $align;
					}
				}
				$current_exon = "";
			}
		}
		if(defined $ref_seq->[0]->[0] && length($ref_seq->[0]->[0]) == 0) {
			print "Gene not found in database!\n";
			return 1;
		}

		if(@{$strains} == 1) {
			if($strand eq "-") {
				$output_seq = reverse($output_seq);
			}
			if($print_gene == 1) {
				if($html == 0) {
					print $output_seq . "\n";
				} else {
					print OUT $output_seq . "<br>";
				}
			}
			if($print_protein == 1) {
				$as = &get_protein_from_seq($output_seq);
				if($html == 0) {
					print $as . "\n";
				} else {
					print OUT $as . "<br>";
				}
			}
			$ref_index = 0;
		}	

		for(my $index = 0; $index < @{$strains}; $index++) {
			for(my $a = 0; $a < $allele; $a++) {
				if($strains->[$index] eq "reference") { 
					$ref_seq->[$a]->[$index] = $output_seq;
					$ref_index = $index;
					$mut_seq->[$a]->[$index] = $output_seq;
					next;
				}
				if($strand eq "-") {
					$comp_seq = "";
					@c = split('', $ref_seq->[$a]->[$index]);
					foreach my $char (@c) {
						$comp_seq = $comp_seq . $comp{$char};
					}
					$ref_seq->[$a]->[$index] = reverse($comp_seq);

					$comp_seq = "";
					@c = split('', $mut_seq->[$a]->[$index]);
					foreach my $char (@c) {
						$comp_seq = $comp_seq . $comp{$char};
					}
					$mut_seq->[$a]->[$index] = reverse($comp_seq);

					$comp_seq = $align_string->[$a]->[$index];
					$align_string->[$a]->[$index] = reverse($comp_seq);
				}
				if($print_gene == 1 && $align_gene == 1) {
					if($html == 0) {
						print "Reference vs " . $strains->[$index] . " for allele " . ($a + 1) . "\n";
						&print_alignment($ref_seq->[$a]->[$index], $mut_seq->[$a]->[$index], $align_string->[$a]->[$index]);
					} else {
						&print_alignment_html($ref_seq->[$a]->[$index], $mut_seq->[$a]->[$index], $align_string->[$a]->[$index], $strains->[$index], $a);
					}
				} elsif($print_gene == 1) {
					$print_seq = $mut_seq->[$a]->[$index];
					$print_seq =~ s/-//g;
					if($html == 0) {
						print $strains->[$index] . " for allele " . ($a + 1) . "\n";
						print $print_seq . "\n";
					} else {
						print OUT $strains->[$index] . " for allele " . ($a + 1) . "<br>";
						print OUT $print_seq . "<br>";
					}
				} 
				if($print_protein == 1 && $align_protein == 1) {
					if($html == 0) {
						print "Reference vs " . $strains->[$index] . " for allele " . ($a + 1) . "\n";
						&align_protein_from_seq($ref_seq->[$a]->[$index], $mut_seq->[$a]->[$index], "reference", $strains->[$index], $a);	
					} else {
						&align_protein_from_seq($ref_seq->[$a]->[$index], $mut_seq->[$a]->[$index], "reference", $strains->[$index], $a, 1);	
					
					}
				 } elsif($print_protein == 1) {
					$as = &get_protein_from_seq($mut_seq->[$a]->[$index]);
					if($html == 0) {
						print "Protein " . $strains->[$index] . " for allele " . ($a + 1) . "\n";
						print $as . "\n";
					} else {
						print OUT $as . "<br>";
					}
					
				}
			}
		}
	}
}

sub get_gene_interspecies{
	my $id = $_[0];
        my $trans = $_[1];
	my $org1 = $genome;
	my $org2 = $_[2];
	my $strains2 = $_[3];
	my $print_gene = $_[4];
	my $align_gene = $_[5];
	my $print_protein = $_[6];
	my $align_protein = $_[7];
	$_ = () for my ($lfit_db, $strand, $inversion, $num);
	#Right order of genomes
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org1 . "_" . $org2 . "\'") || die $DBH::errstr;
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();
	if($num->{'c'} == 0) {
		print "SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org2 . "_" . $org1 . "\'\n";
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org2 . "_" . $org1 . "\'") || die $DBH::errstr;
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == 0) {
			print STDERR "There is no genome alignment available for these two genomes!\n";
			exit;
		} else {
			my $tmp = $org1;
			$org1 = $org2;
			$org2 = $tmp;
			$genome = $org1;
			$tmp = $strains;
			$strains = $strains2;
			$strains2 = $tmp;
		}
	}
	my $start = 0;
	my ($end, $chr1, $chr2, $refs, $muts, $aligns, $final_ref, $final_mut, $final_align, $inter_mut, $lift_db);
	$_ = "" for my($org1_seq, $org2_seq, $inter_ali, $final_org1_seq, $final_org2_seq, $final_inter_ali);

	my $transcript  = &get_transcripts($trans, $id);
	foreach my $trans (keys %{$transcript}) {
		$final_ref = ();
		$final_mut = ();
		$final_align = ();
		if($html == 0) {
			print $trans . "\n";
		} else {
			print OUT "<b>" . $trans . "</b><br><br>";
		}
		#First step - check if the other organism has this gene annotated
		$command = "SELECT DISTINCT id from " . $org2 . "_genes WHERE transcript = \'$trans\'";
		$sth = $dbh->prepare($command);
		$sth->execute();
		if(!defined (my $t = $sth->fetchrow_hashref())) {
			print STDERR "This gene is not annotated in $org2\n";
			print STDERR "Check if genome part exists\n";
		} else {
			print STDER "TODO exists\n";
			exit;
		}
		#Second step - check if other organism even posseses this part of the genome
		$command = "SELECT * from " . $org1 . "_genes WHERE transcript = \'$trans\'";
		$sth = $dbh->prepare($command);
		$sth->execute();
		#Every exon is check individually
		while(my $e = $sth->fetchrow_hashref()) {
			$strand = $e->{'strand'}; 
			my $command2 = "SELECT * from lift_" . $org1 . "_" . $org2 . " WHERE chr1 = \'" . $e->{'chr'} . "\' AND start1 <= " . $e->{'start'} . " AND end1 >= " . $e->{'stop'};
			$lift_db = $dbh->prepare($command2);
			$lift_db->execute();
			if(!defined ($lift = $lift_db->fetchrow_hashref())) {
				print STDERR "The genomic sequence of this gene does not exist in " . $org2 . "\n";
				exit;
			} else {
				print STDERR "Sequence exists for this exon!\n";
				$inversion = $lift->{'inversion'};
				$chr2 = $lift->{'chr2'};
				$chr1 = $lift->{'chr1'};
				my ($start, $mut_start) = $e->{'start'} =~ /^(.*)(.{2})$/s; 
				my ($stop, $mut_stop) = $e->{'stop'} =~ /^(.*)(.{2})$/s; 
				$end = $e->{'start'};	
				my $exon = $e->{'exon'};
				my $s = &get_seq($start, $stop, $e->{'chr'}, $mut_start, $mut_stop, $exon, ($exon + 1));
				$inter_mut = $dbh->prepare("SELECT * from mut_" . $org1 . "_" . $org2 . " WHERE chr1 = \'$e->{'chr'}\' AND pos1 >= $e->{'start'} AND pos1 <= $e->{'stop'}");
				$inter_mut->execute();
				my $ali_start = $e->{'start'};
				$org1_seq = "";
				$org2_seq = "";
				$inter_ali = "";
				my $gap1 = 0;
				my $gap2 = 0;
				while(my $second = $inter_mut->fetchrow_hashref()) {
					$org1_seq .= substr($s, $ali_start - $e->{'start'}, $second->{'pos1'} - $ali_start);
					$org2_seq .= substr($s, $ali_start - $e->{'start'}, $second->{'pos1'} - $ali_start);
					for(my $i = 0; $i < $second->{'pos1'} - $ali_start; $i++) {
						$inter_ali .= "|";
					}

					$org1_seq .= $second->{'mut1'};
					$org2_seq .= $second->{'mut2'};
					if(substr($second->{'mut1'}, 0, 1) eq substr($second->{'mut2'}, 0, 1)) {
						$inter_ali .= "|";
					} else {
						$inter_ali .= ".";
					}

					$end = $second->{'pos1'};
					my $l_mut1 = length($second->{'mut1'});
					while($l_mut1 < length($second->{'mut2'})) {
						$org1_seq .= "-";
						$inter_ali .= "-";
						$l_mut1++;
						$gap1++;
					}
					my $l_mut2 = length($second->{'mut2'});
					while($l_mut2 < length($second->{'mut1'})) {
						$org2_seq .= "-";
						$inter_ali .= "-";
						$l_mut2++;
						$gap1--;
					}
					$ali_start += length($second->{'mut1'});
				}
				$org1_seq .= substr($s, $end - $e->{'start'}, $e->{'stop'} - $e->{'start'});
				$org2_seq .= substr($s, $end - $e->{'start'}, $e->{'stop'} - $e->{'start'});
				for(my $i = $end - $e->{'start'}; $i < $e->{'stop'} - $e->{'start'}; $i++) {
					$inter_ali .= "|";
				}
				
				my $diff = $lift->{'start2'} - $lift->{'start1'};
				if(@{$strains} > 1 && @{$strains2} > 1) {
					for(my $index = 0; $index < @{$strains}; $index++) {
						for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
							for(my $a = 0; $a < $allele; $a++) {
#								print "align: " . $strains->[$index] . "\t" . $strains2->[$index2] . "\n";
								($refs, $muts, $aligns) = &align_seqs($org1, $org1_seq, $strains->[$index], ($a + 1), $e->{'start'}, $e->{'stop'}, $chr1, ($e->{'start'} + $diff), ($e->{'stop'} + $diff), $org2, $org2_seq, $strains2->[$index2], $chr2, $inter_ali);
								$final_ref->[$index]->[$index2]->[$a] .= $refs;
								$final_mut->[$index]->[$index2]->[$a] .= $muts;
								$final_align->[$index]->[$index2]->[$a] .= $aligns;
#								print $strains->[$index] . "\t" . $strains2->[$index2] . "\n";
#								print $refs . "\n" . $muts . "\n" . $aligns . "\n";
							}
						}
					}
				}
				$inter_mut->finish();
			}
		}
		#Check if strand is negative
		if($strand eq "-") {
			for(my $index = 0; $index < @{$strains}; $index++) {
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						$final_ref->[$index]->[$index2]->[$a] = &reverse($final_ref->[$index]->[$index2]->[$a]);
						$final_mut->[$index]->[$index2]->[$a] = &reverse($final_mut->[$index]->[$index2]->[$a]);
						$final_align->[$index]->[$index2]->[$a] = reverse($final_align->[$index]->[$index2]->[$a]);
					}
				}
			}
		}
		if($print_gene == 1 && $align_gene == 0) {
			if($html == 0) {
				for(my $index = 0; $index < @{$strains}; $index++) {
					for(my $a = 0; $a < $allele; $a++) {
						print "Species: " . $org1 . "\nIndividual: " . $strains->[$index] . " (allele: " . ($a + 1) . ")\n";
						my $tmp = $final_mut->[$index]->[0]->[$a];
						$tmp =~ s/-//g;
						print $tmp . "\n";
					}
				}
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						print "Species: " . $org2 . "\nIndividual: " . $strains2->[$index2] . " (allele: " . ($a + 1) . ")\n";
						my $tmp = $final_mut->[0]->[$index2]->[$a];
						$tmp =~ s/-//g;
						print $tmp . "\n";
					}
				}
			} else {
				for(my $index = 0; $index < @{$strains}; $index++) {
					for(my $a = 0; $a < $allele; $a++) {
						my $tmp = $final_mut->[$index]->[0]->[$a];
						$tmp =~ s/-//g;
						&print_html($tmp, "Species: " . $org1 . "<br>Individual: " . $strains->[$index], $a);
					}
				}
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						my $tmp = $final_mut->[0]->[$index2]->[$a];
						$tmp =~ s/-//g;
						&print_html($tmp, "Species: " . $org2 . "<br>Individual: " . $strains2->[$index2], $a);
					}
				}
			}
		} elsif($align_gene == 1) {
			for(my $index = 0; $index < @{$strains}; $index++) {
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						if($html == 0) {
							print STDERR "Species " . $org1 . ": " . $strains->[$index] . " vs. Species " . $org2 . ": " . $strains2->[$index2] . " (allele: " . ($a + 1) . ")\n";
							print "Species " . $org1 . ": " . $strains->[$index] . " vs. Species " . $org2 . ": " . $strains2->[$index2] . " (allele: " . ($a + 1) . ")\n";
							&print_alignment($final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], $final_align->[$index]->[$index2]->[$a]);
						} else {
							&print_alignment_html($final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], $final_align->[$index]->[$index2]->[$a], $strains->[$index], $a, $strains2->[$index2], $org2);
						}
					}
				}
			}	

		}
		if($print_protein == 1 && $align_protein == 0) {
			if($html == 0) {
				for(my $index = 0; $index < @{$strains}; $index++) {
					for(my $a = 0; $a < $allele; $a++) {
						print "Species: " . $org1 . "\nIndividual: " . $strains->[$index] . " (allele: " . ($a + 1) . ")\n";
						print &get_protein_from_seq($final_mut->[$index]->[0]->[$a]) . "\n";
					}
				}
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						print "Species: " . $org2 . "\nIndividual: " . $strains2->[$index2] . " (allele: " . ($a + 1) . ")\n";
						print &get_protein_from_seq($final_mut->[0]->[$index2]->[$a]) . "\n";
					}
				}
			} else {
				for(my $index = 0; $index < @{$strains}; $index++) {
					for(my $a = 0; $a < $allele; $a++) {
						my $tmp = &get_protein_from_seq($final_mut->[$index]->[0]->[$a]);
						&print_html($tmp, "Species: " . $org1 . "<br>Individual: " . $strains->[$index], $a);
					}
				}
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						my $tmp = &get_protein_from_seq($final_mut->[0]->[$index2]->[$a]);
						&print_html($tmp, "Species: " . $org2 . "<br>Individual: " . $strains2->[$index2], $a);
					}
				}
			}	
		} elsif($align_protein == 1) {
			for(my $index = 0; $index < @{$strains}; $index++) {
				for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
					for(my $a = 0; $a < $allele; $a++) {
						if($html == 0) {
							&align_protein_from_seq($final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], "Species: " . $org1 . " - Individual: " . $strains->[$index], "Species: " . $org2 . " - Individual: " . $strains2->[$index2], $a);
						} else {
							&align_protein_from_seq($final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], "Species: " . $org1 . " - Individual: " . $strains->[$index], "Species: " . $org2 . " - Individual: " . $strains2->[$index2], $a, 1);	
						}	
					}
				}
			}
		} 
		$sth->finish();
	}
}

sub print_alignment{
	my $s1 = $_[0];
	my $s2 = $_[1];
	my $align = $_[2];
	my $length = length($s1);
	my $two = length($s2);
	$_ = 0 for my($run, $pos_s1, $pos_s2, $c);
	while($run < $length) {
		if(length($s1) - $run > 100) {
			print "" . ($pos_s1 + 1) . "\t" .  substr($s1, $run, 100) . "\n";
			print "\t" . substr($align, $run, 100) . "\n";
			print "" . ($pos_s2 + 1) . "\t" . substr($s2, $run, 100) . "\n";
			print "\n\n";
			$c = () = substr($s1, $run, 100) =~ /-/gi;
			$c = () = substr($s2, $run, 100) =~ /-/gi;
		} else {
			print "" . ($pos_s1 + 1) . "\t" .  substr($s1, $run) . "\n";
			print "\t" . substr($align, $run) . "\n";
			print "" . ($pos_s2 + 1) . "\t" . substr($s2, $run) . "\n";
			print "\n\n";

		}
		$pos_s1 += (100 - $c);
		$pos_s2 += (100 - $c);
		$run += 100;
	}
}

sub print_alignment_html{
	my $s1 = $_[0];
	my $s2 = $_[1];
	my $align = $_[2];
	my $length = length($s1);
	my $strain = $_[3];
	my $a = $_[4];
	my $strain2 = "Reference";
	my $org2 = "";
	if(@_ > 5) {
		$strain2 = $_[5];
	}
	if(@_ > 6) {
		$org2 = $_[6];
	}
	$_ = 0 for my ($run, $pos_s1, $pos_s2, $c, $col_s1, $col_s2, $col_align);
	$_ = () for my (@work_s1, @work_s2, @work_align);
	if($org2 ne "") {
		print OUT "<b>Species: " . $genome . ": " . $strain2 . " vs. Species: " . $org2 . ": " . $strain . " (allele " . ($a + 1) . ")</b><br>\n";
	} else {
		print OUT "<b>$strain2 vs. " . $strain . " (allele " . ($a + 1) . ")</b><br>\n";

	}
	while($run < $length) {
		if(length($s1) > 100) {
			@work_s1 = split('', substr($s1, $run, 100));
			@work_s2 = split('', substr($s2, $run, 100));
			@work_align = split('', substr($align, $run, 100));
		} else {
			@work_s1 = split('', substr($s1, $run));
			@work_s2 = split('', substr($s2, $run));
			@work_align = split('', substr($align, $run));
		}
		$col_s1 = "";
		$col_s2 = "";
		$col_align = "";

		for(my $i = 0; $i < @work_s1; $i++) {
			if($work_align[$i] ne "|") {
				$col_s1 .= "<font color=\"red\">$work_s1[$i]</font>";
				$col_s2 .= "<font color=\"red\">$work_s2[$i]</font>";
				$col_align .= "<font color=\"red\">$work_align[$i]</font>";
			} else {
				$col_s1 .= $work_s1[$i];
				$col_s2 .= $work_s2[$i];
				$col_align .= $work_align[$i];
			}
		}
		print OUT "<pre>" . ($pos_s1 + 1) . "\t" .  $col_s1 . "</pre>";
		print OUT "<pre>\t" . $col_align .  "</pre>";
		print OUT "<pre>" . ($pos_s2 + 1) . "\t" . $col_s2 . "</pre>";
		print OUT "<br>";
		$c = () = substr($s1, $run, 100) =~ /-/gi;
		$pos_s1 += (100 - $c);
		$c = () = substr($s2, $run, 100) =~ /-/gi;
		$pos_s2 += (100 - $c);
		$run += 100;
	}
	print OUT "<br>\n";
}

sub align_protein_from_seq {
	my $ref_seq = $_[0];
	my $mut_seq = $_[1];
	my $strain1 = $_[2];
	my $strain2 = $_[3];
	my $allele = $_[4];
#	my $html = 0;
#	if(@_ > 4) { $html = 1; }
	my $align_string;
	$_ = 0 for my($kr, $km, $k_ref, $k_mut, $sum_ref, $sum_mut);
	$align_string = "";
	#Save the original sequence with gaps - we always know if there was a gap in the original nucleotide sequence
	my $work_mut = $mut_seq;
	my $work_ref = $ref_seq;
	#Save the seq without gaps so we can translate it into as
	$work_mut =~ s/-//gi;
	$work_ref =~ s/-//gi;
	#Split seq with gaps in triplets - keeping track of gaps
	my @split_mut = unpack("(a3)*", $work_mut);
	my @split_ref = unpack("(a3)*", $work_ref);
	#Split seq without gaps in triplets - translating to as
	my @split_work_mut = unpack("(a3)*", $mut_seq);
	my @split_work_ref = unpack("(a3)*", $ref_seq);
	my $max = (@split_ref > @split_mut) ? @split_ref : @split_mut;
	my $as_mut;
	my $as_ref;
	my ($c_ref, $c_mut);
	for(my $k = 0; $k < $max; $k++) {
		#Count number of gaps so we now when to add a gap in the seq
		$c_ref = () = $split_work_ref[$kr] =~ /\-/gi;
		$c_mut = () = $split_work_mut[$km] =~ /\-/gi;
		$sum_ref += $c_ref;
		$sum_mut += $c_mut;
		if(@split_ref > $k_ref && @split_mut > $k_mut) {
			#Check if working seq is the same
			if($split_work_ref[$kr] ne $split_work_mut[$km]) {
				#Working seq is  a complete gap - add gap to seq
				if($split_work_ref[$kr] eq "---") {
					$align_string .= "-";
					$as_ref .= "-";
					$as_mut .= $nt_to_as{$split_mut[$k_mut]};
					$k_mut++;
					$km++;
					$kr++;
					$sum_ref = $sum_ref - 3;
					next;
				} elsif($split_work_mut[$km] eq "---") {
					$align_string .= "-";
					$as_mut .= "-";
					$as_ref .= $nt_to_as{$split_ref[$k_ref]};
					$k_ref++;
					$kr++;
					$km++;
					$sum_mut = $sum_mut - 3;
					next;
				} else {
					#There were 3 gaps in the nt seq - time to add a gap in the as seq
					if($sum_ref == 3) {
						$align_string .= "-";
						$as_ref .= "-";
						$as_mut .= $nt_to_as{$split_mut[$k_mut]};
						$k_mut++;
						$sum_ref = 0;
					}
					if($sum_mut == 3) {
						$align_string .= "-";
						$as_mut .= "-";
						$sum_mut = 0;
						$as_ref .= $nt_to_as{$split_ref[$k_ref]};
						$k_ref++;
					}
				}
				#Working seq is not the same - but the as did not change - add as to string 
				if($nt_to_as{$split_mut[$k_mut]} eq $nt_to_as{$split_ref[$k_ref]}) {
					$as_ref .= $nt_to_as{$split_ref[$k_ref]};
					$as_mut .= $nt_to_as{$split_mut[$k_mut]};
					$align_string .= "|";
					$k_ref++;
					$k_mut++;
					$kr++;
					$km++;
				#Working seq is not the same - as did change - add to string
				} else {
					$align_string .= ".";
					$as_ref .= $nt_to_as{$split_ref[$k_ref]};
					$as_mut .= $nt_to_as{$split_mut[$k_mut]};
					$k_ref++;
					$k_mut++;
					$kr++;
					$km++;

				}
			#Working seq is the same
			} else {
				#Check if it is still a triplet (often at the end with gaps it does not work out)
				if(!exists $nt_to_as{$split_mut[$k_mut]} && !exists $nt_to_as{$split_ref[$k_ref]}) {
					$as_ref .= "#";
					$as_mut .= "#";
					$align_string .= "-";
					$k_ref++;
					$k_mut++;
					next;
				}
				if(!exists $nt_to_as{$split_mut[$k_mut]}) {
					$as_ref .= $nt_to_as{$split_ref[$k_ref]};
					$as_mut .= "#";
					$align_string .= "-";
					$k_ref++;
					next;
				}
				if(!exists $nt_to_as{$split_ref[$k_ref]}) {
					$as_mut .= $nt_to_as{$split_mut[$k_mut]};
					$as_ref .= "#";
					$align_string .= "-";
					$k_mut++;
					next;
				}
				#As are the same
				if($nt_to_as{$split_mut[$k_mut]} eq $nt_to_as{$split_ref[$k_ref]}) {
					$as_ref .= $nt_to_as{$split_ref[$k_ref]};
					$as_mut .= $nt_to_as{$split_mut[$k_mut]};
					$align_string .= "|";
				} else {
				#AS different - possible because working seq can be the same, but AS not because there were gaps
				#Example: Working seq old: TTG	TT-
				#	  Working seq now: CCA	CCA
				#	  Split seq old:   TTG	TTC
				# 	  Split seq now:   CCA	CAX
				#C => heck if AS are the same
					$align_string .= ".";
					$as_ref .= $nt_to_as{$split_ref[$k_ref]};
					$as_mut .= $nt_to_as{$split_mut[$k_mut]};
				}
				$k_ref++;
				$k_mut++;
				$kr++;
				$km++;
			}
		#Fill as seq with the longer one
		} elsif(@split_mut > $k_mut && length($split_mut[$k_mut]) == 3) {
			$as_mut .= $nt_to_as{$split_mut[$k_mut]};
			$k_mut++;
			next;
		} elsif(@split_ref > $k_ref && length($split_ref[$k_ref]) == 3) {
			$as_ref .= $nt_to_as{$split_ref[$k_ref]};
			$k_ref++;
			next;
		}
	}

	if($html == 1) {
		if($strain1 eq "reference" && $allele == 0) {
			print OUT $trans . "<br>\n";
		}
		&print_alignment_html($as_ref, $as_mut, $align_string, $strain1, $allele, $strain2);
	} else {
		if($strain1 eq "reference") {
			print $trans . "\n"; 
		}
		print STDERR "Reference vs " . $strain1 . " for allele " . ($allele + 1) . "\n";
		&print_alignment($as_ref, $as_mut, $align_string, 0);
	}
	return 1;
}

sub get_protein_from_seq {
	my $seq = $_[0];
	my $as;
	my $mut_seq;
	my $num;
	my $trans_number = 0;

	$mut_seq = $seq;
	$as = "";
	$mut_seq =~ s/-//g;
	my @split_seq = unpack("(a3)*", $mut_seq);
	foreach my $s (@split_seq) {
		if(length($s) == 3) {
			$as = $as . $nt_to_as{$s};
		}
	}
	return $as;
}


sub get_genomic_seq{
	my $start_genome = substr($_[0], 0, length($_[0]) - 2);
	my $stop_genome = substr($_[1], 0, length($_[1]) - 2);
	my $start = $_[0];
	my $stop  = $_[1];
	my $chr = $_[2];
	my $strand = $_[3]; 
	my $align = $_[6];
	my $first = 0;
	my $mut_shift = 0;
	my @general_shift;
	my $alleles;

	if($homo == 1) {
		$alleles = 1;
	} else {
		$alleles = 2;
	}
	my $ref = "";	
	#Get genomic sequence
	my $current_exon = &get_seq($start_genome, $stop_genome, $chr, $start, $stop, 1, 2);
	for(my $index = 0; $index < @{$strains}; $index++) {
		if($strains->[$index] eq "reference") { next; }
		for(my $a = 0; $a < $allele; $a++) {
			my ($refs, $muts, $align) = &align_seqs($genome, $current_exon, $strains->[$index], ($a + 1), $start, $stop, $chr, $start, $stop);
			print $refs . "\n";
			print $muts . "\n";
			print $align . "\n";
		#	$ref_seq->[$a]->[$index] .= $refs;
		#	$mut_seq->[$a]->[$index] .= $muts;
		#	$align_string->[$a]->[$index] .= $align;
		}
	}
	exit;
		#	$current_exon = "";
		#}
#	}


	my $seq;
	my ($old_pos, $exact_start, $exact_stop); 
	while(my $s = $seq->fetchrow_hashref()) {
		while($old_pos + 1 < $s->{'pos'}) {
			for(my $n = 0; $n < 100; $n++) {
				$ref .= "N";
			}
			$old_pos++;
		}
		if($s->{'pos'} == $start_genome) {
			if($s->{'pos'} == $stop_genome) {
				$ref = $ref . substr($s->{'seq'}, $exact_start, $exact_stop - $exact_start);
			} else {
				$ref = $ref . substr($s->{'seq'}, $exact_start);
			}
		} elsif($s->{'pos'} == $stop_genome) {
			$ref = $ref . substr($s->{'seq'}, 0, $exact_stop);
		} else {
			if($old_pos + 1 == $s->{'pos'}) {
				$ref = $ref . $s->{'seq'};
			} else {
				for(my $n = 0; $n < 100; $n++) {
					$ref .= "N";
				}
			}
		}
		$old_pos = $s->{'pos'};
	}
	$seq->finish();

	if($old_pos == 0) {
		for(my $i = $start; $i < $stop; $i++) {
			$ref .= "N";
		}
	}

	$_ = () for my(@output, @strain, @ref, @align, @ass_start);
	for(my $index = 0; $index < @{$strains}; $index++) {
		for(my $a = 0; $a < $alleles; $a++) {
			if($strains->[$index] eq "reference") { $ref[$a][$index] = $ref; $ref_index = $index; next; }
			$ass_start[$a][$index] = 0;
			$output[$a][$index] = "";
			$ref[$a][$index] = "";
			$align[$a][$index] = "";
			$general_shift[$a][$index] = 0;
			my $mut = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains->[$index] . "_allele_" . ($a + 1) . " WHERE pos >= " . $start . " AND pos <= "  . $stop . " AND chr=\'" . $chr . "\'");
			$mut->execute();
			#Insert mutations into genomic sequence
			while(my $m = $mut->fetchrow_hashref()) {
				if($first == 0) {
					my $last = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strains->[$index] . "_allele_" . ($a + 1) . " WHERE index = " . ($m->{'index'} - 1));
					$last->execute();
					my $t = $last->fetchrow_hashref();
					if($t->{'pos'} + length($t->{'reference'}) - 1 >= $start) {
						$mut_shift = length($t->{'reference'}) - length($t->{'strain'});
					}
					$last->finish();
					$first++;
				}
				for(my $i = 0; $i < $mut_shift; $i++) {
					$align[$a][$index] .= "-";
					$output[$a][$index] .= "-";
				}
				$output[$a][$index] .= substr($ref, $ass_start[$a][$index] + $mut_shift, $m->{'pos'} - ($start + $ass_start[$a][$index] + $mut_shift + 1));
				$ref[$a][$index] .= substr($ref, $ass_start[$a][$index], $m->{'pos'} - $start - $ass_start[$a][$index] - 1);
				for(my $i = 0; $i < $m->{'pos'} - $start - $ass_start[$a][$index] - $mut_shift - 1; $i++) {
					$align[$a][$index] .= "|";
				}
				$general_shift[$a][$index] = $mut_shift;
				$mut_shift = 0;
				my $run = length($m->{'reference'}) > length($m->{'strain'}) ? length($m->{'strain'}) : length($m->{'reference'});
				for(my $k = 0; $k < $run; $k++) {
					if(substr($m->{'reference'}, $k, 1) ne substr($m->{'strain'}, $k, 1)) {
						$align[$a][$index] .= ".";
					} else {
						$align[$a][$index] .= "|";
					}
				}
				$output[$a][$index] = $output[$a][$index] . $m->{'strain'};
				$ref[$a][$index] = $ref[$a][$index] . $m->{'reference'};

				$ass_start[$a][$index] = $m->{'pos'} - 1 - $start + length($m->{'reference'});
			#	if($align == 1) {
					if(length($m->{'reference'}) < length($m->{'strain'})) {
						for(my $i = 0; $i < length($m->{'strain'}) - length($m->{'reference'}); $i++) {
							$ref[$a][$index] .= "-";
							$align[$a][$index] .= "-";
						}
			#		} else {
			#			for(my $i = 0; $i < length($m->{'reference'}) - length($m->{'strain'}); $i++) {
			#				$output[$a][$index] .= "-";
			#				$align[$a][$index].= "-";
			#			}
			#		}
				}
			}
			$mut->finish();
			if($strains->[$index] ne "reference") {
				if(length($ref) > $ass_start[$a][$index] + $general_shift[$a][$index]) {
					#Add end of sequence
					$ref[$a][$index] = $ref[$a][$index] . substr($ref, $ass_start[$a][$index]);
					$output[$a][$index] = $output[$a][$index] . substr($ref, $ass_start[$a][$index] + $general_shift[$a][$index]);
					for(my $i = 0; $i < length($ref) - $ass_start[$a][$index]; $i++) {
						$align[$a][$index] .= "|";
					}
				}
			}
	#	}

	#	for(my $index = 0; $index < @{$strains}; $index++) {
			if($strand eq "-") {
				if($strains->[$index] eq "reference") {
					$comp_seq = "";
					@c = split('', $ref[$a][$index]);
					foreach my $char (@c) {
						$comp_seq = $comp_seq . $comp{$char};
					}
					$ref[$a][$index] = reverse($comp_seq);
				} else {
				
					@c = split('', $output[$a][$index]);
					$comp_seq = "";
					foreach my $char (@c) {
						$comp_seq = $comp_seq . $comp{$char};
					}
					$output[$a][$index] = reverse($comp_seq);
					@c = split('', $ref[$a][$index]);
					$comp_seq = "";
					foreach my $char (@c) {
						$comp_seq = $comp_seq . $comp{$char};
					}
					$ref[$a][$index] = reverse($comp_seq);
					$align[$a][$index] = reverse($align[$a][$index]);
				}
			}
		}
	}

	for(my $index = 0; $index < @{$strains}; $index++) {
		for(my $a = 0; $a < $alleles; $a++) {
			if($align == 0 || ($align == 1 && @{$strains} == 1)) {
				if($html == 0) {
					if($index == 0) {
						print "Reference (allele " . ($a + 1) . "):\n";
						print $ref[$a][$ref_index] . "\n";
					}
					if($strains->[$index] eq "reference") { next; }
					print "$strains->[$index] (allele " . ($a + 1) . "):\n";
					print $output[$a][$index] . "\n";
				} else {
					if($index == 0) {
						&print_html($ref[$a][$ref_index], "Reference", $a);
					}
					if($strains->[$index] eq "reference") { next; }
					&print_html($output[$a][$index], $strains->[$index], $a);
				}
			} else {
				if($strains->[$index] eq "reference") { next; }
				if($html == 0) {
					print "Reference vs $strains->[$index] (allele " . ($a + 1) . ")\n";
					&print_alignment($ref[$a][$index], $output[$a][$index], $align[$a][$index]);
				} else {
					&print_alignment_html($ref[$a][$index], $output[$a][$index], $align[$a][$index], $strains->[$index], $a);
				}
			}
		}
	}
	return 1;
}

sub get_genomic_seq_interspecies {
	my $start_genome = substr($_[0], 0, length($_[0]) - 2);
	my $stop_genome = substr($_[1], 0, length($_[1]) - 2);
	my $start = $_[0];
	my $stop  = $_[1];
	my $chr = $_[2];
	my $strand = $_[3]; 
	my $org2 = $_[4];
	my $align = $_[6];
	my $first = 0;
	my $mut_shift = 0;
	my @general_shift;
	my $alleles;
	my $org1 = $genome;
	my $strains2;
	$_ = () for my ($lfit_db, $inversion, $num);
	#Right order of genomes
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org1 . "_" . $org2 . "\'") || die $DBH::errstr;
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();
	if($num->{'c'} == 0) {
		print "SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org2 . "_" . $org1 . "\'\n";
		$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org2 . "_" . $org1 . "\'") || die $DBH::errstr;
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == 0) {
			print STDERR "There is no genome alignment available for these two genomes!\n";
			exit;
		} else {
			my $tmp = $org1;
			$org1 = $org2;
			$org2 = $tmp;
			$genome = $org1;
			$tmp = $strains;
			$strains = $strains2;
			$strains2 = $tmp;
		}
	}
	return 1;
}

sub print_html {
	my $run = 0;
	my $seq = $_[0];
	my $strain = $_[1];
	my $a = $_[2];
	print OUT "<FONT FACE=\"Courier\"><font size = 3><b>" . $strain . " allele " . ($a+1) . "</b></font><br>\n";
	while($run < length($seq) - 100) {
		print OUT substr($seq, $run, 100) . "<br>\n";
		$run = $run + 100;
	}
	print OUT substr($seq, $run) . "<br><br></FONT>\n";
}

sub align_seqs {
	my $org1 = $_[0];
	my $seq1 = $_[1];
	my $strain = $_[2];
	my $allele = $_[3];
	my $mut_start = $_[4];
	my $mut_stop = $_[5];
	my $chr1 = $_[6];
	my $ref_start = $_[7]; #ref->{'start'}
	my $ref_stop = $_[8];
	my $org2 = "";
	my $seq2 = "";
	my $chr2 = "";
	my $strain2 = "";
	my $ali_seq = "";
	if(@_ > 9) {
		$org2 = $_[9];
		$seq2 = $_[10];
		$strain2 = $_[11];
		$chr2 = $_[12];
		$ali_seq = $_[13];
	}
	
	$_ = 0 for my ($exon, $start, $stop, $new_start, $index, $first, $mut_shift, $run);
	$_ = "" for my ($output_seq, $strand, $seq, $mut, $chr);	
	$_ = () for my (@all_seqs, $ref_seq, $diff, $end, $ass, $align_string, $mut_seq, $last_exon, $general_shift, @w1, @w2);
	if($org2 ne "" && $org1 ne $org2) {
		if($strain ne "reference") {
			$mut = $dbh->prepare("SELECT * FROM " . $org1 . "_mutations_" . $strain . "_allele_" . $allele . " WHERE pos >=" . $mut_start . " AND pos <= " . $mut_stop . " AND chr=\'" . $chr1 . "\'");
			$mut->execute();
			my $mut_pos = 0;
			my $work_pos = 0;
			while(my $m = $mut->fetchrow_hashref()) {
				$mut_pos = $m->{'pos'} - $mut_start - $work_pos;
				$ref_seq .= substr($seq1, $work_pos, $mut_pos) . $m->{'reference'};
				$align_string .= substr($ali_seq, $work_pos, $mut_pos);
			#	$mut_seq .= substr($seq2, $work_pos, $mut_pos) . substr($seq2, $mut_pos + $work_pos, 1);
				$mut_seq .= substr($seq2, $work_pos, $mut_pos) . $m->{'strain'};
				if(substr($m->{'strain'}, 0, 1) eq substr($seq2, $mut_pos + $work_pos, 1)) {
					$align_string .= "|";
				} else {
					$align_string .= ".";
				}
				my $ref_l = 1;
				my $str_l = 1;
				while(length($m->{'reference'}) > $ref_l) {
					$ref_seq .= "-";
					$align_string .= "-";
					$mut_seq .= substr($seq2, $mut_pos + $work_pos, 1);
					$mut_pos++;
					$ref_l++;
				}
				while(length($m->{'strain'}) > $str_l) {
					$mut_seq .= "-";
					$align_string .= "-";
					$str_l++;
				}
				$work_pos = $mut_pos + 1;
			}
			$ref_seq .= substr($seq1, $mut_pos);
			$align_string .= substr($ali_seq, $mut_pos);
			$mut_seq .= substr($seq2, $mut_pos);
			$mut->finish();
		} elsif($strain2 ne "reference") {
			#Get mutations for organism two
			$mut = $dbh->prepare("SELECT * FROM " . $org2 . "_mutations_" . $strain2 . "_allele_" . $allele . " WHERE pos >=" . $ref_start . " AND pos <= " . $ref_stop . " AND chr=\'" . $chr2 . "\'");
			$mut->execute();
			my $mut_pos = 0;
			my $work_pos = 0;
			while(my $m = $mut->fetchrow_hashref()) {
				@w2 = split('', $seq2);
				$mut_pos =  $m->{'pos'} - $ref_start - $work_pos;
				my $run = 0;
				my $num_mut1 = 0;
				my $num_mut2 = 0;
				#Get number of gaps in this shift before the mutation happens, because the mutation in the strains is annotated from the reference of this organism, not from the other organism
				#TODO use shift vector
				my $gaps = $dbh->prepare("SELECT * FROM mut_" . $org1 . "_" . $org2 . " WHERE chr1 = \'$chr1\' AND chr2 = \'$chr2\' AND pos1 > " . $lift->{'start1'} . " AND pos1 < " . $mut_start . " AND (length(mut1) > 1 OR length(mut2) > 1)");
				$gaps->execute();
				while(my $g = $gaps->fetchrow_hashref()) {
					$num_mut1 += length($g->{'mut1'});
					$num_mut2 += length($g->{'mut2'});
				}
				#Add beginning of seq and mutation
				$ref_seq .= substr($seq1, $work_pos, $mut_pos) . $m->{'reference'};
				$align_string .= substr($ali_seq, $work_pos, $mut_pos);
				$mut_seq .= substr($seq2, $work_pos, $mut_pos) . $m->{'strain'};
				#Create alignment sequence
				if(substr($m->{'strain'}, 0, 1) eq substr($seq1, $mut_pos, 1)) {
					$align_string .= "|";
				} else {
					$align_string .= ".";
				}
				#Add gaps
				my $ref_l = 1;
				my $str_l = 1;
				while(length($m->{'reference'}) > $ref_l) {
					$mut_seq .= "-";
					$align_string .= "-";
					$ref_l++;
				}
				while(length($m->{'strain'}) < $str_l) {
					$ref_seq .= "-";
					$align_string .= "-";
					$str_l++;
				}
				#Shift to new position
				$work_pos  = $mut_pos + length($m->{'reference'});
			}
			#Add the end
			$ref_seq .= substr($seq1, $work_pos);
			$align_string .= substr($ali_seq, $work_pos);
			$mut_seq .= substr($seq2, $work_pos);
		} else {
			$ref_seq = $seq1;
			$mut_seq = $seq2;
			$align_string = $ali_seq;
		}
	} else {
		$first = 0;
		$mut = $dbh->prepare("SELECT * FROM " . $org1 . "_mutations_" . $strain . "_allele_" . $allele . " WHERE pos >" . $mut_start . " AND pos <= " . $mut_stop . " AND chr=\'" . $chr1 . "\'");
		$mut->execute();
		$ass = 0;
		$diff = 0;
		while(my $m = $mut->fetchrow_hashref()) {
			if($m->{'pos'} == $ref_stop) {
				next;
			}
			if($first == 0) {
				my $last = $dbh->prepare("SELECT * FROM " . $genome . "_mutations_" . $strain . "_allele_" . $allele . " WHERE index = " . ($m->{'index'} - 1));
				$last->execute();
				my $t = $last->fetchrow_hashref();
				if($t->{'pos'} + length($t->{'reference'}) - 1 >= $mut_start) {
					$mut_shift = length($t->{'reference'}) - length($t->{'strain'});
				}
				$last->finish();
				$first++;
			}
			$diff = $m->{'pos'} - $ref_start - 1;
			if($diff == -1) { $diff++; }
			if($diff + length($m->{'reference'}) > length($seq1)) {
				$end = length($seq1) - $diff;
			} else {
				$end = length($m->{'reference'})
			}

			for(my $k = 0; $k < $mut_shift; $k++) {
				$align_string .= "-";
				$mut_seq .= "-";
			}
			$mut_seq = $mut_seq. substr($seq1, $ass + $mut_shift, $diff - ($ass + $mut_shift));
			$ref_seq = $ref_seq . substr($seq1, $ass, $diff - $ass);
			for(my $k = 0; $k < ($diff - $ass - $mut_shift); $k++) {
					$align_string .= "|";
			}
			$ref_seq = $ref_seq. $m->{'reference'};
			$mut_seq = $mut_seq . $m->{'strain'};
			my $run = length($m->{'reference'}) > length($m->{'strain'}) ? length($m->{'strain'}) : length($m->{'reference'});

			for(my $k = 0; $k < $run; $k++) {
				if(substr($m->{'reference'}, $k, 1) eq substr($m->{'strain'}, $k, 1)) {
					$align_string .= "|";
				} else {
					$align_string .= ".";
				}
			}
			#Reference seq from DB is longer than strain - Deletion in strain - add gaps to strain
			if(length($m->{'strain'}) < length($m->{'reference'})) {
				for(my $k = 0; $k < (length($m->{'reference'}) - length($m->{'strain'})); $k++) {
					$mut_seq .= "-";
					$align_string .= "-";
				#	$ass++;
				}
			}
			#Strain seq from DB is longer than reference - Deletion in reference - add gaps to reference
			if(length($m->{'reference'}) < length($m->{'strain'})) {
				for(my $j = 0; $j < (length($m->{'strain'}) - length($m->{'reference'})); $j++) {
					$ref_seq .= "-";
					$align_string .= "-";
				}
			}
			$ass = $diff + length($m->{'reference'});
			if($ass > length($seq1)) {
				$ass = length($seq1);
			}
			$general_shift = $mut_shift;
			$mut_shift = 0;
		}
		$mut->finish();

		$mut_seq = $mut_seq . substr($seq1, $ass);
		$ref_seq .= substr($seq1, $ass);
		for(my $k = $ass; $k < length($seq1); $k++) {
			$align_string .= "|";
		}
		if($mut_seq eq "") {
			$mut_seq = $seq1;
			$ref_seq = $seq1;
		}
	}
	return($ref_seq, $mut_seq, $align_string);
}


sub get_seq{
	my $start = $_[0];
	my $stop = $_[1];
	my $chr = $_[2];
	my $mut_start = $_[3];
	my $mut_stop = $_[4];
	my $exon = $_[5];
	my $last_exon = $_[6];
	my $seq = $dbh->prepare("SELECT * FROM " . $genome . "_ref_genome WHERE pos >= " . $start . " AND pos <= " . $stop . " AND chr=\'" . $chr . "\'");
	$seq->execute();
	#Get all the mutations which are in this part of the sequence
	#Create reference sequence
	my $current_exon = "";
	while(my $s = $seq->fetchrow_hashref()) {
		#Get the first part, if the complete sequence lies within one sequence part     
		if($s->{'pos'} eq $start && $start == $stop) {
			$current_exon = $current_exon . substr($s->{'seq'}, $mut_start%100, $mut_stop - $mut_start);
		#Get the first part from the start of the exon to the end of the first sequence line
		} elsif($s->{'pos'} eq $start) {
			if ($exon == 1) {
				$current_exon = $current_exon . substr($s->{'seq'}, $mut_start%100, 100 - $mut_start%100 + 1);

			} else {
				$current_exon = $current_exon . substr($s->{'seq'}, $mut_start%100, 100 - $mut_start%100);
			}
		#Get the last part from the beginning of the sequence line to the end of the exon
		} elsif($s->{'pos'} eq $stop) {
			if($exon == $last_exon) {
				$current_exon = $current_exon . substr($s->{'seq'}, 0, ($mut_stop%100));
			} else {
				$current_exon = $current_exon . substr($s->{'seq'}, 0, ($mut_stop%100));
			}
		#Add the lines in the middle of the exon that span a complete sequence line
		} else {
			$current_exon = $current_exon . $s->{'seq'};
		}
	}
	return $current_exon;
}

sub reverse{
	my $comp_seq = "";
	my @c = split('', $_[0]);
	foreach my $char (@c) {
		$comp_seq = $comp_seq . $comp{$char};
	}
	my $seq = reverse($comp_seq);
	return $seq;
}

sub initialize_data{
        (%comp) = (     'A' => 'T',
                        'T' => 'A',
                        'C' => 'G',
                        'G' => 'C',
			'-' => '-');

        (%nt_to_as) = ('TCA'=>'S', #Serine
                'TCC'=>'S', #Serine
                'TCG'=>'S', #Serine
                'TCT'=>'S', #Serine
                'TTC'=>'F', #Phenylalanine
                'TTT'=>'F', #Phenylalanine
                'TTA'=>'L', #Leucine
                'TTG'=>'L', #Leucine
                'TAC'=>'Y', #Tyrosine
                'TAT'=>'Y', #Tyrosine
                'TAA'=>'_', #Stop
                'TAG'=>'_', #Stop
                'TGC'=>'C', #Cysteine
                'TGT'=>'C', #Cysteine
                'TGA'=>'_', #Stop
                'TGG'=>'W', #Tryptophan
                'CTA'=>'L', #Leucine
                'CTC'=>'L', #Leucine
                'CTG'=>'L', #Leucine
                'CTT'=>'L', #Leucine
                'CCA'=>'P', #Proline
                'CAT'=>'H', #Histidine
                'CAA'=>'Q', #Glutamine
                'CAG'=>'Q', #Glutamine
                'CGA'=>'R', #Arginine
                'CGC'=>'R', #Arginine
                'CGG'=>'R', #Arginine
                'CGT'=>'R', #Arginine
                'ATA'=>'I', #Isoleucine
                'ATC'=>'I', #Isoleucine
                'ATT'=>'I', #Isoleucine
                'ATG'=>'M', #Methionine
                'ACA'=>'T', #Threonine
                'ACC'=>'T', #Threonine
                'ACG'=>'T', #Threonine
                'ACT'=>'T', #Threonine
                'AAC'=>'N', #Asparagine
                'AAT'=>'N', #Asparagine
                'AAA'=>'K', #Lysine
                'AAG'=>'K', #Lysine
                'AGC'=>'S', #Serine#Valine
                'AGT'=>'S', #Serine
                'AGA'=>'R', #Arginine
                'AGG'=>'R', #Arginine
                'CCC'=>'P', #Proline
                'CCG'=>'P', #Proline
                'CCT'=>'P', #Proline
                'CAC'=>'H', #Histidine
                'GTA'=>'V', #Valine
                'GTC'=>'V', #Valine
                'GTG'=>'V', #Valine
                'GTT'=>'V', #Valine
                'GCA'=>'A', #Alanine
                'GCC'=>'A', #Alanine
                'GCG'=>'A', #Alanine
                'GCT'=>'A', #Alanine
                'GAC'=>'D', #Aspartic Acid
                'GAT'=>'D', #Aspartic Acid
                'GAA'=>'E', #Glutamic Acid
                'GAG'=>'E', #Glutamic Acid
                'GGA'=>'G', #Glycine
                'GGC'=>'G', #Glycine
                'GGG'=>'G', #Glycine
                'GGT'=>'G', #Glycine
                );
}
1;
