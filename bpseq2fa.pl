#!/usr/bin/perl
#
# bpseq2fa.pl
# Eric Nawrocki
# EPN, Mon Jan 23 16:22:13 2017
#
# Usage: perl bpseq2fa.pl
#             <concatenated bpseq file from Comparative RNA Website (CRW)>
#             
# Synopsis:
# Converts the bpseq file into a FASTA sequence file. DOES NOT 
# validate bpseq file.

use strict;

my $usage = "Usage: perl bpseq2fa.pl\n\t<concatenated bpseq file from Comparative RNA Website (CRW)>\n\n";

if(scalar(@ARGV) != 1)
{
  print $usage;
  exit();
}

my ($bpseq_file) = (@ARGV);

# read the bpseq file and output fasta as we go
open(IN,  $bpseq_file) || die "ERROR unable to open bpseq file $bpseq_file";
my $filename = "";
my $organism = "";
my $accn = "";
my $seq  = "";
while(my $line = <IN>) { 
  chomp $line;
  if($line =~ s/^Filename\:\s+//) {
    if($seq ne "") { 
      printf(">%s:%s:%s\n%s\n", $filename, $organism, $accn, $seq);
      $organism = "";
      $accn = "";
      $seq = "";
    }
    $filename = $line;
    $filename =~ s/\s+/\-/g;
  }
  elsif($line =~ s/^Organism\:\s+//) {
    $organism = $line;
    $organism =~ s/\s+/\-/g;
  }
  elsif($line =~ s/^Accession Number\:\s*//) {
    $accn = $line;
    $accn =~ s/\s+/\-/g;
  }
  elsif($line =~ s/^Citation.+//) {
    ;
  }
  elsif($line =~ m/^(\d+)\s+(\S+)\s+(\d+)/) { #only care about the sequence
    my $i = $1;
    my $a = $2;
    my $j = $3;
    $seq .= $a;
  }
  else { 
    die "ERROR unable to parse line: $line";
  }
}
close(IN);

#output final seq
if($seq ne "") { 
  printf(">%s:%s:%s\n%s\n", $filename, $organism, $accn, $seq);
  $organism = "";
  $accn = "";
  $seq = "";
}

