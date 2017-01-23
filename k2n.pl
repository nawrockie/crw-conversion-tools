while($line = <>) { 
    chomp $line;
    $path = $line;
    $file = $line;
    $file =~ s/^.+\///;
    $nfile = "nested." . $file;
    print("Processing $file ... ");
    system("python ~/src/k2n_standalone/knotted2nested.py -m OSP -f bpseq $path > $nfile\n");
    print("done.\n");
}

