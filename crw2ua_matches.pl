#!/usr/bin/perl
#
# crw2ua_matches.pl
# Eric Nawrocki
# EPN, Thu May  7 09:56:57 2009
#
# Usage: perl crw2ua_matches.pl
#             <concatenated bpseq file from Comparative RNA Website (CRW)>
#             <CRW alignment (aligned fasta) with matching seqs to bpseq data>
#             <root name for output files>
#             
# Synopsis:
# Creates an unambiguous set of matches from the CRW '.bpseq' data and alignment data. 
# For each sequence in the concatenated '.bpseq' file check for a matching accession 
# in the alignment. If it exists, verify they are either (i) identical sequences or
# (ii) the aligned sequence is a subsequence of the '.bpseq' sequence. If so, this is
# an unambiguous match. 
#
# A new stocholm formatted alignment file using individual secondary structures 
# is created for the unambigous matches. This file is called <$out_root>.n<$nmatches>.indi.stk
# A 'summary' file is also created which gives information on which .bpseq sequences
# were skipped and why. 

$usage = "Usage: perl crw2ua_matches.pl\n\t<concatenated bpseq file from Comparative RNA Website (CRW)>\n\t<CRW alignment (aligned fasta) with matching seqs to bpseq data>\n\t<root name for output files>\n\nt";

if(@ARGV != 3)
{
    print $usage;
    exit();
}

($bpseq_file, $aln_file, $out_root) = @ARGV;

@aln2printA = ();
push(@aln2printA, "# STOCKHOLM 1.0\n\n");

$sum_file = $out_root . ".sum";
open(SUM, ">" . $sum_file);

# read aligned fasta file
%aname_H = ();
%aseq_H = ();
open(IN,  $aln_file);
$accn = "";
while($line = <IN>) {
    chomp $line;
    if($line =~ m/^>/) { 
	@elA = split("::", $line);
	if(scalar(@elA) != 3) { die "ERROR, expected three tokens divided by \"::\" in seq name, read $line.\n"; }
	($num, $name, $accn) = @elA;
	if(($accn !~ m/DIVIDER/) && ($accn !~ m/REFERENCE/)) { 
	    $anameline_H{$accn} = $line;
	    $anameline_H{$accn} =~ s/^\>//;
	    $aname_H{$accn} = $name;
	    $aseq_H{$accn}  = "";
	    $a_total++;
	}
    }
    elsif($line =~ m/\S/) { #seq line
	if($accn eq "")     { die "ERROR, first line should be a sequence name line beginning with \>\n"; }
	if($line =~ m/\s+/) { die "ERROR, there is whitespace in a sequence line: $line\n"; }
	if(exists($aseq_H{$accn})) { $aseq_H{$accn} .= $line; }
    }
}
#create unaligned seqs, and capitalize everything
foreach $accn (keys (%aseq_H)) { 
    $uaseq_H{$accn} = $aseq_H{$accn};
    $uaseq_H{$accn} =~ s/\W+//g;
    $uaseq_H{$accn} =~ tr/a-z/A-Z/;
    $aseq_H{$accn} =~ tr/a-z/A-Z/;
    #printf("$accn aln length: %d\n", length($aseq_H{$accn})); 
}

%bpseq_H = ();
# read the bpseq file
$nseq = 0;
@bpseq_accnA = ();
open(IN,  $bpseq_file);
while($line = <IN>) { 
    chomp $line;
    if($line =~ s/^Filename\:\s+//) {
	$filename = $line;
    }
    elsif($line =~ s/^Organism\:\s+//) {
	$organism = $line;
    }
    elsif($line =~ s/^Accession Number\:\s*//) {
	$baccn = $line;
	$bpseq_H{$baccn} = "";
	if($baccn !~ m/\w/) { 
	    $b_nonmatches++; 
	    printf SUM ("%2d. Ignoring .bpseq file: %-30s (the accession field is blank).\n", $b_nonmatches, $filename);
	    printf("%2d. Ignoring .bpseq file: %-30s (the accession field is blank).\n", $b_nonmatches, $filename);
	}
	push(@bpseq_accnA, $baccn);
	@{$ua_ctHA{$baccn}} = ();
	$ua_ctHA{$baccn}[0] = 0;
	$filename_H{$baccn} = $filename;
    }
    elsif($line =~ s/^Citation.+//) {
	;
    }
    elsif($line =~ m/^(\d+)\s+(\S+)\s+(\d+)/) { #only care about the sequence
	$i = $1;
	$a = $2;
	$j = $3;
	$bpseq_H{$baccn} .= $a;
	$ua_ctHA{$baccn}[$i] = $j;
    }
}

#create unaligned seqs, and capitalize everything
foreach $baccn (keys (%bpseq_H)) { 
    $bpseq_H{$baccn} =~ tr/a-z/A-Z/;
}

#For each bpseq accession, search for matching aligned sequence accessions
# There can be multiple accessions for a bpseq sequence, each is separated by whitespace
# There can be multiple accessions for a aligned sequence, but there's no good rule to delimit
#  them, so we ask if any of the bpseq accn are a subseq of the full aligned sequence accession
#  (which may contain > 1 accn concatenated together).
#  
foreach $baccn (@bpseq_accnA) { 
    if($baccn !~ m/\w/) { next; } #skip any data read from files with blank accessions 
    #printf("BACCN $baccn\n");
    $b_total++;
    @baA = split(/\s+/, $baccn);
    $found_match = 0;
    foreach $ba (@baA) { 
	foreach $aa (keys (%aseq_H)) { 
	    #skip seqs that have 0 As, Cs, Gs, Us (that are not seqs)
	    if($aseq_H{$aa} !~ m/[ACGU]/) { next; }

	    if($aa =~ m/$ba/) { 
		# test seq match, bpseq should be superseq of aligned seq
		$uaseq = $uaseq_H{$aa};

		if($bpseq_H{$baccn} =~ m/$uaseq/) { 
		    $first_hit = 1;
		    if(exists($ua_b2aH{$baccn})) { 
			#two alignment accession match a single bpseq accession
			#make sure the aligned sequences are 100% identical and same length
			if($aseq_H{$aa} ne $aseq_H{$ua_b2aH{$baccn}}) {
			    die "ERROR, two alignment accession match the same bpseq accession: $ua_b2aH{$baccn} and $aa\nbut their aligned seqs are not 100\% identical\n"; }
			#if they are identical, move on to next sequence
			$first_hit = 0;
			$a_dups++;
		    }		    
		    if(!($first_hit)) { next; }

		    $ua_b2aH{$baccn} = $aa;
		    $found_match = 1;

		    #create map 
		    @offsetA = ();
		    $offsetA[0] = 0;
		    @aseqA = split("", $aseq_H{$aa});
		    $alen = scalar(@aseqA);
		    $ngaps_seen = 0;
		    $first_matching_posn = index($bpseq_H{$baccn}, $uaseq) + 1;

		    if($first_matching_posn != 1) { 
			;#printf("WHOA first_matching_posn: $first_matching_posn for $aa $baccn\n");
		    }

		    $uapos = 1;
		    $ualen = length($uaseq);
		    $apos = 1;
		    while($uapos <= $ualen) { 
			if($aseqA[($apos-1)] =~ m/\w/) { 
			    #printf("uapos: $uapos offset: $ngaps_seen\n");
			    $offsetA[$uapos] = $ngaps_seen;
			    $uapos++; 
			    $apos++;
			}
			else { 
			    $ngaps_seen++; 
			    $apos++;
			}
		    }
		    #use offset map to update the ct array for this bpseq sequence
		    @$a_ctHA = ();
		    for($apos = 0; $apos <= $alen; $apos++) { 
			$a_ctHA[$apos] = 0;
		    }
		    #printf("fmp: $first_matching_posn\n");
		    for($uapos = 1; $uapos <= $ualen; $uapos++) { 
			$orig_uapos = $uapos + $first_matching_posn - 1;
			$apos = $uapos + $offsetA[$uapos];
			#$apos = $orig_uapos + $offsetA[$uapos];
			#printf("uapos: $uapos apos: $apos\n");
			$orig_j = $ua_ctHA{$baccn}[$orig_uapos];
			if(($orig_j == 0) || 
			   (($orig_j > ($ualen + $first_matching_posn - 1)) || #if the original bp has right half occurs after (3') of the unaligned subsequence that matches the bpseq, ignore it
			   (($orig_j < $first_matching_posn)))) {              #if the original bp has left half occurs before (5') the end of the unaligned subsequence that matches the bpseq, ignore it
			    $j = 0; 
			    $a_ctHA[$apos] = $j;
			} 
			else {
			    $j = $orig_j - $first_matching_posn + 1;
			    $a_ctHA[$apos] = $j + $offsetA[$j];
			}
			#printf("a_ctHA[$apos] = $a_ctHA[$apos]\t\t$orig_uapos $ua_ctHA{$baccn}[$orig_uapos]\t\tOJ: $orig_j ualen: $ualen\n");
		    }		    
		    $ss = "";
		    for($apos = 1; $apos <= $alen; $apos++) { 
			$j = $a_ctHA[$apos];
			#printf("$apos $j\n");
			if($j == 0)       { $ss .= "."; }
			elsif($j > $apos) { $ss .= "<"; }
			else              { $ss .= ">"; }
		    }
		    $name2print = $anameline_H{$aa};
		    #$name2print =~ s/\s+/\:\:/; 
		    $seq2print  = sprintf("%-85s $aseq_H{$aa}\n", $name2print);
		    $ss2print  = sprintf("#=GR %-80s $ss\n", $name2print . " SS");
		    push(@aln2printA, $seq2print);
		    push(@aln2printA, $ss2print);
		    $b_matches++;
		}
	    }
	}
	if(!($found_match)) { 
	    $b_nonmatches++; 
	    #printf("No alignment accession match for $filename_H{$baccn} (accn: $baccn). It will be ignored.\n");
	    printf ("%2d. Ignoring .bpseq file: %-30s (no alignment accession match for accn: $baccn).\n", $b_nonmatches, $filename_H{$baccn});
	    printf SUM ("%2d. Ignoring .bpseq file: %-30s (no alignment accession match for accn: $baccn).\n", $b_nonmatches, $filename_H{$baccn});
	    
	}
    }
}
printf("\n");
#if($b_nonmatches > 0) { die "ERROR, at least 1 bpseq sequence does not match to any alignment accessions."; }

    push(@aln2printA, "\/\/\n");

#print alignment file
$aln_file = $out_root . ".indi.stk";
open(OUT, ">" . $aln_file);
foreach $line (@aln2printA) { 
    printf OUT ("$line"); 
}
close(OUT);

printf SUM ("Matched %4d/%4d bpseq sequences to %4d/%4d aligned sequences (%4d duplicates)\n", $b_matches, ($b_total), $b_matches, ($a_total), $a_dups);
printf SUM ("New %4d sequence stockholm alignment file with individual secondary structure\nannotation printed to $aln_file\n", $b_matches);
printf ("Matched %4d/%4d bpseq sequences to %4d/%4d aligned sequences (%4d duplicates)\n", $b_matches, ($b_total), $b_matches, ($a_total), $a_dups);
printf ("New %4d sequence stockholm alignment file with individual secondary structure\nannotation printed to $aln_file\n", $b_matches);
close(SUM);
printf ("The information printed above was saved to file %s\n", $sum_file);

