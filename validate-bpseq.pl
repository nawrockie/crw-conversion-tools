$ctr = 0;
while($line = <>) { 
  $ctr++;
  # First 4 lines should look like this:
  # Filename: d.23.e.S.cerevisiae.bpseq
  # Organism: Saccharomyces cerevisiae
  # Accession Number: U53879
  # Citation and related information available at http://www.rna.icmb.utexas.edu
  if($ctr == 1) { 
    if($line !~ m/^Filename\:/) { die "ERROR line 1 unexpected format; $line"; }
  }
  if($ctr == 2) { 
    if($line !~ m/^Organism\:/) { die "ERROR line 2 unexpected format; $line"; }
  }
  if($ctr == 3) { 
    if($line !~ m/^Accession Number\:/) { die "ERROR line 3 unexpected format; $line"; }
  }
  if($ctr == 4) { 
    if($line !~ m/^Citation and related information/) { die "ERROR line 4 unexpected format; $line"; }
  }

  if($line =~ /^(\d+)\s+\S\s+(\d+)/) { 
    ($i, $j) = ($1, $2);
    if($j > 0) { 
      if($i > $j) { 
        if($bp_iH{$j} != $i) { die "ERROR inconsistent case 1 $line\n"; }
      }
      else { 
        $bp_iH{$i} = $j;
        $bp_jH{$j} = $i;
      }
    } 
    else { # $j is 0
      if(exists($bp_jH{$i})) { die "ERROR inconsistent case 2 $line\n"; }
    }
  }
}
