#!/usr/bin/perl -w

package database_interaction;
use strict;
use DBI;
use Data::Dumper;
require '../general/config.pm';

#Define variables
$_ = () for my (%comp, %nt_to_as, %mut_per_gene);
my ($dbh, $sth, $command, $trans, $id, $ref_index, $comp_seq, @c, %seen);

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
		if($name ne "") {
			open OUT, ">$name.html";
			print STDERR "Writing output to $name.html\n";
		} else {
			open OUT, ">genomic_sequence.html";
			print STDERR "Writing output to genomic_sequence.html\n";
		}
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
  	my ($refs, $muts, $align); 
	$_ = 0 for my ($exon, $start, $stop, $new_start, $mut_start, $mut_stop, $index, $first, $mut_shift, $run);
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
					for(my $a = 0; $a < $allele; $a++) {
						if($strains->[$index] eq "reference") { 
							$ref_seq->[$a]->[$index] .= $current_exon;
							$mut_seq->[$a]->[$index] .= $current_exon;
						} else {
							($refs, $muts, $align) = &align_seqs($genome, $current_exon, $strains->[$index], ($a + 1), $mut_start, $mut_stop, $chr, $ref->{'start'}, $ref->{'stop'});
							$ref_seq->[$a]->[$index] .= $refs;
							$mut_seq->[$a]->[$index] .= $muts;
							$align_string->[$a]->[$index] .= $align;
						}
					}
				}
				$current_exon = "";
			}
		}
		if(defined $ref_seq->[0]->[0] && length($ref_seq->[0]->[0]) == 0) {
			print "Gene not found in database!\n";
			return 1;
		}
	}

	for(my $index = 0; $index < @{$strains}; $index++) {
		for(my $a = 0; $a < $allele; $a++) {
			if($strand eq "-") {
				($refs, $muts, $align) = &reverse($strains->[$index], $ref_seq->[$a]->[$index], $mut_seq->[$a]->[$index], $align_string->[$a]->[$index]);
			}
			&print_routine($index, $a, $print_gene, $align_gene, $print_protein, $align_protein, $refs, $muts, $align);
		}
	}
	return 1;
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
	my ($end, $chr1, $chr2, $refs, $muts, $aligns, $final_ref, $final_mut, $final_align, $inter_mut, $lift_db, $align);
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
			#	print STDERR "Sequence exists for this exon!\n";
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
								($refs, $muts, $aligns) = &align_seqs($org1, $org1_seq, $strains->[$index], ($a + 1), $e->{'start'}, $e->{'stop'}, $chr1, ($e->{'start'} + $diff), ($e->{'stop'} + $diff), $org2, $org2_seq, $strains2->[$index2], $chr2, $inter_ali);
								$final_ref->[$index]->[$index2]->[$a] .= $refs;
								$final_mut->[$index]->[$index2]->[$a] .= $muts;
								$final_align->[$index]->[$index2]->[$a] .= $aligns;
							}
						}
					}
				}
				$inter_mut->finish();
			}
		}

		for(my $index = 0; $index < @{$strains}; $index++) {
			for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
				for(my $a = 0; $a < $allele; $a++) {
					if($strand eq "-") {
						($final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], $final_align->[$index]->[$index2]->[$a]) = &reverse($strains->[$index], $final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], $final_align->[$index]->[$index2]->[$a]);
					}
					&print_routine_interspecies($index, $a, $print_gene, $align_gene, $print_protein, $align_protein, $final_ref->[$index]->[$index2]->[$a], $final_mut->[$index]->[$index2]->[$a], $final_align->[$index]->[$index2]->[$a], $org1, $org2, $strains2->[$index2]);
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
	if(@_ > 6) {
		$strain2 = $_[6];
	}
	if(@_ > 5) {
		$org2 = $_[5];
	}
	$_ = 0 for my ($run, $pos_s1, $pos_s2, $c, $col_s1, $col_s2, $col_align);
	$_ = () for my (@work_s1, @work_s2, @work_align);
	if($org2 ne "") {
		print OUT "<b>Species: " . $genome . " - " . $strain . " vs. Species: " . $org2 . " - " . $strain2 . " (allele " . ($a + 1) . ")</b><br>\n";
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
#		print @work_s1 . "\t" . @work_align . "\t" . @work_s2 . "\n";
		if(@work_s1 != @work_align) { 
			print $s1 . "\n" . $align . "\n" . $s2 . "\n";
			print "WHAT THE FUCK\n";
			exit; 
		}
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
					if(($sum_ref > 0 || $sum_mut > 0) && ($k_mut + 3 > @split_mut || $k_ref + 3 > @split_ref)) {
						$align_string .= "-";
						if($sum_ref > 0 && $k_ref + 3 > @split_ref) {
							$as_ref .= "-";
						} else {
							$as_ref .= $nt_to_as{$split_ref[$k_ref]};
						}
						if($sum_mut > 0 && $k_mut + 3 > @split_mut) {
							$as_mut .= "-";
						} else {
	                                                $as_mut .= $nt_to_as{$split_mut[$k_mut]};
						}
						
						return ($as_ref, $as_mut, $align_string);
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
	return ($as_ref, $as_mut, $align_string);
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
	my $print_gene = $_[4];
	my $align_gene = $_[5];
	my $print_protein = $_[6];
	my $align_protein = $_[7];
	if(@_ > 7) {
		$strains = $_[8];
	}
	my $alleles;
	my $return;

	if($homo == 1) {
		$alleles = 1;
	} else {
		$alleles = 2;
	}
	my ($refs, $muts, $align) = "";	
	#Get genomic sequence
	my $current_exon = &get_seq($start_genome, $stop_genome, $chr, $start, $stop, 1, 2);
	for(my $index = 0; $index < @{$strains}; $index++) {
		for(my $a = 0; $a < $allele; $a++) {
			if($strains->[$index] eq "reference") {  
				$refs = $current_exon;
			} else {
				($refs, $muts, $align) = &align_seqs($genome, $current_exon, $strains->[$index], ($a + 1), $start, $stop, $chr, $start, $stop);
			}
		}
		if($strand eq "-") {
			($refs, $muts, $align) = &reverse($strains->[$index], $refs, $muts, $align);
		}
		for(my $a = 0; $a < $alleles; $a++) {
			&print_routine($index, $a, $print_gene, $align_gene, $print_protein, $align_protein, $refs, $muts, $align);
		}
		if($align_gene == 1) {
			$return = "reference\n" . $refs . "\n" . $strains->[$index] . "\n" . $muts . "\n";
		} else {
			if($strains->[$index] eq "reference") {
				$return .= $strains->[$index] . "\n" . $refs . "\n";
			} else {
				$return .= $strains->[$index] . "\n" . $muts . "\n";
			}
		}
	}
	return $return;
}





sub print_routine{
	my $index = $_[0];
	my $a = $_[1];
	my $print_gene = $_[2];
	my $align_gene = $_[3];
	my $print_protein = $_[4];
	my $align_protein = $_[5];
	my $ref = $_[6];
	my $mut = $_[7];
	my $align = $_[8];
	my $org2;
	my $strains2;
	my ($print_seq, $as, $as_2, $as_align);

	if(@{$strains} == 1) {
		if($print_gene == 1) {
			if($html == 0) {
				print "Reference (allele: " . ($a + 1) . ")\n";
				print $ref . "\n";
			} else {
				print OUT "Reference (allele: " . ($a + 1) . ")<br>";
				print OUT $ref . "<br>";
			}
		}
		if($print_protein == 1) {
			$as = &get_protein_from_seq($ref);
			if($html == 0) {
				print "Reference (allele: " . ($a + 1) . ")\n";
				print $as . "\n";
			} else {
				print "Reference (allele: " . ($a + 1) . ")<br>";
				print OUT $as . "<br>";
			}
		}
		$ref_index = 0;
		return 1;
	}	

	if($print_gene == 1 && $align_gene == 0) {
		if($strains->[$index] eq "reference") {
			$print_seq = $ref;
		} else {
			$print_seq = $mut;
		}
		$print_seq =~ s/-//g;
		if($html == 0) {
			print $strains->[$index] . " (allele " . ($a + 1) . ")\n";
			print $print_seq . "\n";
		} else {
			print OUT $strains->[$index] . " (allele " . ($a + 1) . ")<br>";
			print OUT $print_seq . "<br>";
		}
	} elsif($print_gene == 1 && $align_gene == 1) {
		if($strains->[$index] eq "reference") { return 1; }
		if($html == 0) {
			print "Reference vs " . $strains->[$index] . " for allele " . ($a + 1) . "\n";
			&print_alignment($ref, $mut, $align);
		} else {
			&print_alignment_html($ref, $mut, $align, $strains->[$index], $a);
		}
	}

	if($print_protein == 1 && $align_protein == 0) { 
		$as = &get_protein_from_seq($mut);
		if($html == 0) {
			print "Protein " . $strains->[$index] . " for allele " . ($a + 1) . "\n";
			print $as . "\n";
		} else {
			print OUT $as . "<br>";
		}
	} elsif($print_protein == 1 && $align_protein == 1) {
		if($strains->[$index] eq "reference") { return 1; }
		my @a = &align_protein_from_seq($ref, $mut, "reference", $strains->[$index], $a);     
		$as = $a[0];
		$as_2 = $a[1];
		$as_align = $a[2];
		if($html == 0) {
			print "Reference vs " . $strains->[$index] . " for allele " . ($a + 1) . "\n";
			&print_alignment($as, $as_2, $as_align);	
		} else {
			&print_alignment_html($as, $as_2, $as_align, $strains->[$index], $a);
		}
	}
}

sub print_routine_interspecies{
	my $index = $_[0];
	my $a = $_[1];
	my $print_gene = $_[2];
	my $align_gene = $_[3];
	my $print_protein = $_[4];
	my $align_protein = $_[5];
	my $ref = $_[6];
	my $mut = $_[7];
	my $align = $_[8];
	my $org1 = $_[9];
	my $org2 = $_[10];
	my $strain2 = $_[11];
	my ($print_seq, $as, $as_2, $as_align);
	if($print_gene == 1 && $align_gene == 0) {
		if(!exists $seen{$org1 . "_" . $strains->[$index]}) {
			$print_seq = $ref;
			$print_seq =~ s/-//g;
			if($html == 0) {
				&print_seqs($print_seq, "Species: " . $org1 . " - " . $strains->[$index], $a);
			} else {
				&print_html($print_seq, "Species: " . $org1 . " - " . $strains->[$index], $a);
			}
			$seen{$org1 . "_" . $strains->[$index]} = 1;
		}
		if(!exists $seen{$org2 . "_" . $strain2}) {
			$print_seq = $mut;
			$print_seq =~ s/-//g;
			if($html == 0) {
				&print_seqs($print_seq, "Species: " . $org2 . " - " . $strain2, $a);
			} else {
				&print_html($print_seq, "Species: " . $org2 . " - " . $strain2, $a);
			}
			$seen{$org2 . "_" . $strain2} = 1;
		}

	} elsif($print_gene == 1 && $align_gene == 1) {
		if($html == 0) {
			print "Species " . $org1 . " - " . $strains->[$index] . " vs Species " . $org2 . " - " . $strain2 . " (allele " . ($a + 1) . ")\n";
			&print_alignment($ref, $mut, $align);
		} else {
			&print_alignment_html($ref, $mut, $align, $strains->[$index], $a, $org2, $strain2);
		}
	}

	if($print_protein == 1 && $align_protein == 0) { 
		$as = &get_protein_from_seq($mut);
		if($html == 0) {
			&print_seqs($as, "Protein species " . $org1 . " - " . $strains->[$index], $a);
		} else {
			&print_html($as, "Protein species " . $org1 . " - " . $strains->[$index], $a);
		}
	} elsif($print_protein == 1 && $align_protein == 1) {
	#	if($strains->[$index] eq "reference") { return 1; }
		my @a = &align_protein_from_seq($ref, $mut, $strains->[$index], $strain2, $a);     
		$as = $a[0];
		$as_2 = $a[1];
#		print "\n\nSpecies " . $org1 . " - " . $strains->[$index] . " vs species " . $org2 . " - " . $strain2 . " (allele " . ($a + 1) . ")\n";
#		print $ref . "\n" . $align . "\n" . $mut . "\n";
		$as_align = $a[2];
#		print $as . "\n" . $as_align . "\n" . $as_2 . "\n";
#		if(length($as) != length($as_align)) { 
#			print length($as) . "\t" . length($as_align) . "\t" . length($as_2) . "\n\n\n";
#			print "EXIT\n"; exit;
#		}
		if($html == 0) {
			print "Species " . $org1 . " - " . $strains->[$index] . " vs species " . $org2 . " - " . $strain2 . " (allele " . ($a + 1) . ")\n";
			&print_alignment($as, $as_2, $as_align);	
		} else {
			&print_alignment_html($as, $as_2, $as_align, $strains->[$index], $a, $org2, $strain2);
		}
	}
}

sub reverse {
	my $strain = $_[0];
	my $ref = $_[1];
	my $mut = $_[2];
	my $align = $_[3];
	my $comp;
	my @c;
	if($strain eq "reference") {
		$comp_seq = "";
		@c = split('', $ref);
		foreach my $char (@c) {
			$comp_seq = $comp_seq . $comp{$char};
		}
		$ref = reverse($comp_seq);
		$mut = reverse($comp_seq);
	} else {
		@c = split('', $mut);
		$comp_seq = "";
		foreach my $char (@c) {
			$comp_seq = $comp_seq . $comp{$char};
		}
		$mut = reverse($comp_seq);
		@c = split('', $ref);
		$comp_seq = "";
		foreach my $char (@c) {
			$comp_seq = $comp_seq . $comp{$char};
		}
		$ref = reverse($comp_seq);
		$align = reverse($align);
	}
	return($ref, $mut, $align);
}

sub get_genomic_seq_interspecies {
	my $start_genome = substr($_[0], 0, length($_[0]) - 2);
	my $stop_genome = substr($_[1], 0, length($_[1]) - 2);
	my $start = $_[0];
	my $stop  = $_[1];
	my $chr = $_[2];
	my $strand = $_[3]; 
	my $org2 = $_[4];
	my $strains2 = $_[5];
	my $print_gene = $_[6];
	my $align_gene = $_[7];
	my $print_protein = $_[8];
	my $align_protein = $_[9];
	my $first = 0;
	my $alleles;
	my $org1 = $genome;
	my ($seq1, $seq2, $command, $lift_db, $chr2, $chr1, $start2, $stop2, $shift_db, $result, $diff);

	$_ = () for my ($lfit_db, $inversion, $num);
	#Right order of genomes
	$sth = $dbh->prepare("SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org1 . "_" . $org2 . "\'") || die $DBH::errstr;
	$sth->execute();
	$num = $sth->fetchrow_hashref();
	$sth->finish();
	if($num->{'c'} == 0) {
#		print "SELECT COUNT(*) AS c FROM pg_tables WHERE tablename=\'lift_" . $org2 . "_" . $org1 . "\'\n";
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
		$sth = $dbh->prepare("SELECT count(*) AS c from lift_" . $org1 . "_" . $org2 . " WHERE chr2 = \'". $chr . "\' AND start2 <= " . $start . " AND end2 >= " . $stop);
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == 0) {
			print STDERR "This sequence has no corresponding sequence in $org2!\n";
			exit;
		}
		$command = "SELECT * from lift_" . $org1 . "_" . $org2 . " WHERE chr2 = \'" . $chr . "\' AND start2 <= " . $start . " AND end2 >= " . $stop;
	} else {
		$command = "SELECT * from lift_" . $org1 . "_" . $org2 . " WHERE chr1 = \'" . $chr . "\' AND start1 <= " . $start . " AND end1 >= " . $stop;
		$sth = $dbh->prepare("SELECT count(*) AS c from lift_" . $org1 . "_" . $org2 . " WHERE chr1 = \'". $chr . "\' AND start1 <= " . $start . " AND end1 >= " . $stop);
		$sth->execute();
		$num = $sth->fetchrow_hashref();
		$sth->finish();
		if($num->{'c'} == 0) {
			print STDERR "This sequence has no corresponding sequence in $org2!\n";
			exit;
		}

	}
#	print $command . "\n";
	$lift_db = $dbh->prepare($command);
	$lift_db->execute();
	my ($shift, $ref_org1, $ref_org2, $ref_align);
	while($lift = $lift_db->fetchrow_hashref()) {
		$inversion = $lift->{'inversion'};
		$chr2 = $lift->{'chr2'};
		$chr1 = $lift->{'chr1'};
		$start2 = $lift->{'start2'};
		$stop2 = $lift->{'end2'};
		if($inversion == 1) {
			print "INVERSION: " . $inversion . "\n";
			print "NEXT\n";
			next;
		}
		#Get shift vector
		$command = "SELECT * from offset_" . $org1 . "_" . $org2 . " WHERE chr1 = \'" . $chr1 . "\' AND chr2 = \'" . $chr2 . "\' AND start1 < " . $start . " AND start2 > " . $start2 . " AND start2 < " . $stop2 . " ORDER BY start1 DESC LIMIT 1";
		$shift_db = $dbh->prepare($command);
		$shift_db->execute();
		if(defined ($result = $shift_db->fetchrow_hashref())) {
			$shift = $result->{'shift'};
		} else {
			$shift = 0;
		}
		$diff = $lift->{'start2'} - $lift->{'start1'} + $shift;
		$start2 = $start + $diff;
		$stop2 = $stop + $diff;
		$seq1 = &get_seq(substr($start, 0, length($start) - 2), substr($stop, 0, length($stop) - 2), $chr, $start, $stop, 1, 2);
		$seq2 = &get_seq(substr($start2, 0, length($start2) - 2), substr($stop2, 0, length($stop2) - 2), $chr2, $start2, $stop2, 1, 2, $org2);
		my ($ref_org1, $ref_org2, $ref_align) = &align_seqs($org1, $seq1, "reference", 1, $start, $stop, $chr1, $start2, $stop2, $org2, $seq2, "reference", $chr2, "");
		if($html == 1) {
			print OUT "<b>Org1 chromosome " . $chr1 . " (" . $start . " - " . $stop . ") TO org2 chromosome " . $chr2 . " ($start2 - $stop2)</b><br><br><br>";
		} else {
			print "Org1 chromosome " . $chr1 . " (" . $start . " - " . $stop . ") TO org2 chromosome " . $chr2 . " ($start2 - $stop2)\n";

		}
		for(my $index = 0; $index < @{$strains}; $index++) {
			for(my $index2 = 0; $index2 < @{$strains2}; $index2++) {
				for(my $a = 0; $a < $allele; $a++) {
					if($strains->[$index] eq "reference" && $strains2->[$index2] eq "reference") {
						&print_routine_interspecies($index, $a, $print_gene, $align_gene, $print_protein, $align_protein, $ref_org1, $ref_org2, $ref_align, $org1, $org2, $strains2->[$index2]);
					}
					my ($refs, $muts, $align) = &align_seqs($org1, $ref_org1, $strains->[$index], ($a + 1), $start, $stop, $chr1, $start2, $stop2, $org2, $ref_org2, $strains2->[$index2], $chr2, $ref_align);
					&print_routine_interspecies($index, $a, $print_gene, $align_gene, $print_protein, $align_protein, $refs, $muts, $align, $org1, $org2, $strains2->[$index2]);
				}
			}
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

sub print_seqs {
	my $run = 0;
	my $seq = $_[0];
	my $strain = $_[1];
	my $a = $_[2];
	print $strain . " allele " . ($a+1) . "\n";
	while($run < length($seq) - 100) {
		print substr($seq, $run, 100) . "\n";
		$run = $run + 100;
	}
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
	$_ = "" for my ($output_seq, $strand, $seq, $mut, $chr, $align_string);	
	$_ = () for my (@all_seqs, $ref_seq, $diff, $end, $ass, $mut_seq, $last_exon, $general_shift, @w1, @w2);
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
			#Get the alignment seq
			$mut = $dbh->prepare("SELECT * FROM mut_" . $org1 . "_" . $org2 . " WHERE chr1 = \'$chr1\' AND chr2 = \'$chr2\' AND pos1 >= " . $mut_start. " AND pos1 < " . $mut_stop . " AND pos2 >= " . $ref_start . " AND pos2 < " . $ref_stop);
			my $mut_pos = 0;
			my $work_pos = 0;
			my $org1_pos = 0;
			my $org2_pos = 0;
			my $last = $mut_start;
			$mut->execute();
			while(my $m = $mut->fetchrow_hashref()) {
				$mut_pos = $m->{'pos1'} - $last;
				$last = $m->{'pos1'} + length($m->{'mut1'});
				$ref_seq .= substr($seq1, $org1_pos, $mut_pos) . $m->{'mut1'};
				$mut_seq .= substr($seq1, $org2_pos, $mut_pos) . $m->{'mut2'};
				for(my $i = 0; $i < $mut_pos; $i++) {
					$align_string .= "|";
					$org1_pos++;
					$org2_pos++;
				}
				if(substr($m->{'mut1'}, 0, 1) eq substr($m->{'mut2'}, 0, 1)) {
					$org1_pos++;
					$org2_pos++;
					$align_string .= "|";
				} else {
					$org1_pos++;
					$org2_pos++;
					$align_string .= ".";
				}
				my $ref_l = 1;
				my $str_l = 1;
				while(length($m->{'mut1'}) > $ref_l) {
					$mut_seq .= "-";
					$align_string .= "-";
					$org1_pos++;
					$ref_l++;
				}
				while(length($m->{'mut2'}) > $str_l) {
					$ref_seq .= "-";
					$org2_pos++;
					$align_string .= "-";
					$str_l++;
				}
			}
			$ref_seq .= substr($seq1, $org1_pos);
			for(my $i = $org1_pos; $i < $mut_stop - $mut_start; $i++) {
				$align_string .= "|";
			} 
			$mut_seq .= substr($seq1, $org1_pos);
			$mut->finish();
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
				if($t->{'pos'} + length($t->{'reference'}) - 1 > $mut_start) {
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

sub get_multiple_alignment{
	my $start_genome = substr($_[0], 0, length($_[0]) - 2);
	my $stop_genome = substr($_[1], 0, length($_[1]) - 2);
	my $start = $_[0];
	my $stop  = $_[1];
	my $chr = $_[2];
	my $strand = $_[3]; 
	my $print_gene = $_[4];
	my $align_gene = $_[5];
	my $print_protein = $_[6];
	my $align_protein = $_[7];
	my $alleles;
	my $return;
	my $i = 0;
	my $save_strains;
	my @split_ref;
	my @split_strain;
	my @ref_seqs;
	my @strain_seqs;
	@{$save_strains} = @{$strains};
	my $length = 0;
	for(my $i = 0; $i < @{$save_strains}; $i++) {
		if($save_strains->[$i] eq "reference") { next; }
		my @temp_strains = ();
		my $t_strains = \@temp_strains;
		push ($t_strains, $save_strains->[$i]);
		my $seq = get_genomic_seq($start, $stop, $chr, "+", 0, 1, 0, 0, $t_strains);
		my @seq_strain = split('\n', $seq);
		$ref_seqs[$i] = $seq_strain[1];
		$strain_seqs[$i] = $seq_strain[3];
	}
	my @final_ref;
	my @final_strain;

	for(my $i = 0; $i < @ref_seqs; $i++) {
		@{$split_ref[$i]} = split('', $ref_seqs[$i]);
		@{$split_strain[$i]} = split('', $strain_seqs[$i]);
		$final_ref[$i] = "";
		$final_strain[$i] = "";
	}
	my %gap;
	while($split_ref[-1]->[0] ne "") {
		%gap = ();
		for(my $j = 0; $j < @ref_seqs; $j++) {
			if($split_ref[$j]->[0] eq "-") {
				$gap{$j} = 1;
			}
		}
		if(keys %gap > 0) {
			for(my $j = 0; $j < @ref_seqs; $j++) {
				if(exists $gap{$j}) {
					$final_ref[$j] .= $split_ref[$j]->[0];
					$final_strain[$j] .= $split_strain[$j]->[0];
					shift $split_ref[$j];
					shift $split_strain[$j];
				} else {
					$final_ref[$j] .= "-";
					$final_strain[$j] .= "-";
				}
			}
		} else {
			for(my $j = 0; $j < @ref_seqs; $j++) {
				$final_ref[$j] .= $split_ref[$j]->[0];
				$final_strain[$j] .= $split_strain[$j]->[0];
				shift $split_ref[$j];
				shift $split_strain[$j];
			}
		}
	}
	my @return;
	for(my $i = 0; $i < @final_ref; $i++) {
		if($save_strains->[$i] eq "reference") { next; }
		push(@return, $final_ref[$i]);
		last;
	}
	for(my $i = 0; $i < @final_ref; $i++) {
		if($save_strains->[$i] eq "reference") { next; }
		push(@return, $final_strain[$i]);
	}
	return \@return;
}

sub get_seq{
	my $start = $_[0];
	my $stop = $_[1];
	my $chr = $_[2];
	my $mut_start = $_[3];
	my $mut_stop = $_[4];
	my $exon = $_[5];
	my $last_exon = $_[6];
	my $tmp_genome = $genome;
	if(@_ > 7) {
		$tmp_genome = $_[7];
	}
	my $seq = $dbh->prepare("SELECT * FROM " . $tmp_genome . "_ref_genome WHERE pos >= " . $start . " AND pos <= " . $stop . " AND chr=\'" . $chr . "\'");
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

#sub reverse{
#	my $comp_seq = "";
#	my @c = split('', $_[0]);
#	foreach my $char (@c) {
#		$comp_seq = $comp_seq . $comp{$char};
#	}
#	my $seq = reverse($comp_seq);
#	return $seq;
#}

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
