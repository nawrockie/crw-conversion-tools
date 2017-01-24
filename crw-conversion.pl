#!/usr/bin/perl
#
# crw-conversion.pl
# Eric Nawrocki
# EPN, Wed May 20 07:15:02 2009
#
# Usage: perl crw-conversion.pl
#             <family name>
#             <output dir for new files>
#             <parameter file minimally defining:
#                $cmalign = path to and options for 'cmalign' infernal app
#                $cmbuild = path to and options for 'cmbuild' infernal app
#                $esl_construct = path to 'esl-construct' easel miniapp
#                $esl_alimanip  = path to 'esl-alimanip' easel miniapp
#                $esl_alistat   = path to 'esl-alistat' easel miniapp
#                $esl_reformat  = path to 'esl-reformat' easel miniapp
#                $crw2ua_matches = path to 'crw2ua_matches' perl script>
#             
# Synopsis:
# Converts Comparative RNA Website structural RNA data (.bpseq) and 
# alignments (.alnfasta) to consensus structure annotated stockholm alignments
# and iteratively refines them using cmalign. 
#
# Prerequisites:
# '$fam.nested.bpseq' must exist in current directory. 
#
use Getopt::Std;
use Cwd 'abs_path';
#use strict;

my $usage = 
    "Usage: crw-conversion.pl [-options]\n\t<family name>\n\t<CRW aligned fasta file>\n\t<NESTED .bpseq file from CRW>\n\t<output dir for new files>\n\t<parameters file>\n\n";
#    "<family name>\n\t" . 
#    "<CRW aligned fasta file>\n\t" . 
#    "<output dir for new files>\n\t" . 
#    "<parameter file minimally defining:\n\t" . 
#    "           \$cmalign = path to and options for 'cmalign' infernal app\n\t"
#    "           \$cmbuild = path to and options for 'cmbuild' infernal app\n\t"
#    "           \$esl_construct = path to 'esl-construct' easel miniapp\n\t"
#    "           \$esl_alimanip  = path to 'esl-alimanip' easel miniapp\n\t"
#    "           \$esl_reformat  = path to 'esl-reformat' easel miniapp\n\t"
#    "           \$crw2ua_matches = path to 'crw2ua_matches' perl script>\n\n"
$options_usage  = "where options are:\n";
$options_usage .= "  -h     : show brief help on version and usage\n";
$options_usage .= "  -F     : allow file clobbering\n";
$options_usage .= "  -i <n> : set max number of cmalign iterations to <n>\n";
$options_usage .= "  -n <x> : set min fraction of median length allowed [df:0.8] (0.0 for no limit)\n";
$options_usage .= "  -x <x> : set max fraction of median length allowed [df:1.2] (10000.0 for no limit)\n";
$options_usage .= "  -a <n> : set max number of ambiguous residues allowed to <n>\n";
$options_usage .= "  -y <x> : use cmbuild --symfrac <x> [df:none]\n";
$options_usage .= "  -w     : use default cmbuild entropy weighting [df:use --enone]\n";
$options_usage .= "  -e <x> : use cmbuild --ere <x> [df:use --enone]\n";
$options_usage .= "  -s     : use cmalign --sub [df:do not use it]\n";
$options_usage .= "  -c <n> : max \# of basepair conflicts to allow during iterative refinement [df:15]\n";
$options_usage .= "  -o     : add --v1p0 to cmbuild call, to build using Infernal 1.0 parameterization\n\n";

###################
# Process options #
###################
getopts('hFi:a:n:x:y:we:sc:o');
my $do_clobber = 0;
my $max_iter = 10;
my $xambig = 5;
my $lnfract = 0.8;
my $lxfract = 1.2;
my $cmb_sf  = "";
#my $cmb_ent = "--enone";
my $cma_opts = "-g";
my $maxconflicts = 15;
my $cmb_1p0 = "";

if (defined $opt_h) { print $usage . "\n"; print $options_usage; exit(1); } 
if (defined $opt_F) { $do_clobber = 1; } 
if (defined $opt_i) { $max_iter = $opt_i; }
if (defined $opt_a) { $xambig   = $opt_a; }
if (defined $opt_n) { $lnfract  = $opt_n; }
if (defined $opt_x) { $lxfract  = $opt_x; }
if (defined $opt_y) { $cmb_sf   = "--symfrac $opt_y"; }
if (defined $opt_w) { $cmb_ent  = ""; }
if (defined $opt_e) { $cmb_ent  = "--ere $opt_e"; }
if (defined $opt_s) { $cma_opts = "--sub --notrunc -g"; }
if (defined $opt_c) { $maxconflicts  = $opt_c; }
if (defined $opt_o) { $cmb_1p0  = "--v1p0"; }

# Check for incompatible option combinations and non-sensical option parameters
if(scalar(@ARGV) != 5) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print $options_usage;
  exit(1);
}
my ($fam, $crwalnfa_file, $bpseq_file, $out_dir, $params_file) = @ARGV;
my $out_stk = "";
my $cmbuild_options = "$cmb_sf $cmb_ent $cmb_1p0";
my $cmalign_options = "$cma_opts";
my $out_dir_root = $out_dir . "/" . $fam;
my $log_file = $out_dir_root . ".log";

########################################################
# Preliminary steps, before cmbuild/cmalign iterations #
########################################################
printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
printf("Preliminary steps\n\n");
printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
validate_input_files(\$cmalign, \$cmbuild, \$esl_construct, \$esl_alimanip, \$esl_reformat, \$crw2ua_matches);
create_output_dir($do_clobber, $out_dir); 

open(LOG, ">" . $log_file);
printf LOG ("                                                      prefilter step   \n");
printf LOG ("                                                 ----------------------\n");
printf LOG ("round     nseq  fminx  ncbp    nbp   ncf    fcf   nrm  namb  nsrt  nlng\n");
printf LOG ("------    ----  -----  ----  -----  ----  -----  ----  ----  ----  ----\n");

# Find matching seqs in CRW structures (.bpseq) and alignment (.alnfasta)
match_step($crwalnfa_file, "$fam.0");
$ntotal_init = determine_num_seqs_in_stk("$out_dir_root.0.indi.stk");

# Determine consensus structure of initial alignment, and --fmin <x>
printf LOG ("input   ");
consensusize_step("$out_dir_root.0.indi.stk", "$out_dir_root.input.cons.stk", "$out_dir_root.input.construct");

# Remove some sequences based on length and number of ambiguities
prefilter_step("$out_dir_root.0.indi.stk", "$out_dir_root.0.fil.indi.stk", \$nrm_ambig, \$nrm_short, \$nrm_long, $nrm_total);

printf("\nPrefilter removed %3d/%3d sequences\n\t(%2d had > %2d ambiguous residues)\n\t(%2d were < %4.2f fraction median length)\n\t(%2d were > %4.2f fraction median length)\n\n", 
       $nrm_total, $ntotal_init, $nrm_ambig, $maxambig, $nrm_short, $lnfract, $nrm_long, $lxfract);
printf LOG ("  %4d  %4d  %4d  %4d\n", $nrm_total, $nrm_ambig, $nrm_short, $nrm_long);

# Create a consensus structure based on the individual structure annotations
printf LOG ("iter 0  ");
consensusize_step("$out_dir_root.0.fil.indi.stk", "$out_dir_root.0.cons.stk", "$out_dir_root.0.construct");
printf LOG ("\n");

############################################################
# Iterative alignment refinement using cmbuild and cmalign #
############################################################
$iter = 1;
$stabilized = 0;
while(($iter <= $max_iter) && (!($stabilized))) { 
  $nseq_left = determine_num_seqs_in_stk("$out_dir_root." . ($iter-1) . ".cons.stk");
  printf("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  printf("Iteration %2d of cmbuild -> cmalign alignment refinement (%3d/%3d seqs remain)\n\n", $iter, $nseq_left, $ntotal_init);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  cmbuild_step(     "$out_dir_root." .($iter-1) . ".cons.stk", 
                    "$out_dir_root." . $iter    . ".cm", 
                    "$out_dir_root." . $iter    . ".cmbuild");

  fastaize_step(    "$out_dir_root." .($iter-1) . ".cons.stk", 
                    "$out_dir_root." . $iter . ".fa");

  cmalign_step(     "$out_dir_root." . $iter . ".cm", 
                    "$out_dir_root." . $iter . ".fa", 
                    "$out_dir_root." . $iter . ".cmalign.stk",
                    "$out_dir_root." . $iter . ".cmalign");

  alnfastaize_step ("$out_dir_root." . $iter . ".cmalign.stk",
                    "$out_dir_root." . $iter . ".cmalign.alnfasta");		      

  match_step       ("$out_dir_root." . $iter . ".cmalign.alnfasta",
                    "$fam.$iter");

  printf LOG ("iter $iter  ");
  consensusize_step("$out_dir_root." . $iter . ".indi.stk", 
                    "$out_dir_root." . $iter . ".nonfil.cons.stk", 
                    "$out_dir_root." . $iter . ".nonfil.construct");
  printf LOG ("\n");
  
  $nremoved = rmconflicts_step("$out_dir_root." . $iter . ".nonfil.cons.stk", 
                               "$out_dir_root." . $iter . ".cons.rm" . $maxconflicts. ".list", 
                               "$out_dir_root." . $iter . ".cons.stk");

  printf("iter: %2d.  %2d seqs with > %2d conflicts removed.\n", $iter, $nremoved, $maxconflicts);
  if($nremoved == 0) { 
    $stabilized = 1; 
    printf("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    run_command("cp $out_dir_root.$iter.cons.stk $out_dir_root.final.stk");
    printf("\n$out_dir_root.final.stk is your seed alignment.\n\n");
    printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  }
  $iter++;
}
printf LOG ("\n");
close(LOG);

# success, exit normally
exit(0); 

###############
# subroutines #
###############

#####################################################################
# subroutine: match_step
# incept:     EPN, Wed May 20 11:52:31 2009
# 
# purpose:    Match CRW sequences we have bpseq data for to their 
#             aligned versions in the aligned fasta file using
#             crw2ua_matches.pl (see that script for more details).
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub match_step { 
  $narg_expected = 2;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, match_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($alnfa_file, $root) = @_;

  $system_call = "perl $crw2ua_matches $bpseq_file $alnfa_file $root > /dev/null";
  $system_call2print = $system_call;
  $system_call2print =~ s/^perl \S+\//perl /;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");
  run_command("mv $root.indi.stk $out_dir");
  run_command("mv $root.sum $out_dir");
  return;
}


#####################################################################
# subroutine: prefilter_step
# incept:     EPN, Wed May 20 11:52:31 2009
# 
# purpose:    Remove sequences that fail our prefilter quality 
#             requirements using esl-alimanip.
#             -  sequences must be >= <x1> fraction the median length 
#                of the sequences in the alignment
#             - sequences must be <= <x2> fraction the median length 
#                of the sequences in the alignment
#             - sequences must have less than or equal to <x3> 
#               ambiguous bases. 
#
#               default values: 
#               <x1> = 0.8 (changed to <x> with -n <x>)
#               <x2> = 1.2 (changed to <x> with -x <x>)
#               <x3> = 5   (changed to <n> with -a <n>)
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub prefilter_step { 
  $narg_expected = 6;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, prefilter_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_stk, $out_stk, $nrm_ambig_ref, $nrm_short_ref, $nrm_long_ref, $nrm_total_ref) = @_;

  if($lxfract > 999) { 
    $system_call = "$esl_alimanip -o $out_stk --xambig $xambig --lnfract $lnfract $in_stk";
  }
  else { 
    $system_call = "$esl_alimanip -o $out_stk --xambig $xambig --lnfract $lnfract --lxfract $lxfract $in_stk";
  }
  $system_call2print = $system_call;
  $system_call2print =~ s/^\S+\///;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");

  #rerun esl-alimanip many times wastefully, just to figure out how many were removed at each step
  #determine number of sequences initially 
  $nseqs_init = determine_num_seqs_in_stk($in_stk);

  $nrm_total = $nseqs_init - determine_num_seqs_in_stk($out_stk);

  $system_call = "$esl_alimanip -o $fam.tmp.stk --xambig $xambig $in_stk";
  run_command("$system_call");
  $nrm_ambig = $nseqs_init - determine_num_seqs_in_stk("$fam.tmp.stk");

  $system_call = "$esl_alimanip -o $fam.tmp.stk --lnfract $lnfract $in_stk";
  run_command("$system_call");
  $nrm_short = $nseqs_init - determine_num_seqs_in_stk("$fam.tmp.stk");

  if($lxfract < 999) { 
    $system_call = "$esl_alimanip -o $fam.tmp.stk --lxfract $lxfract $in_stk";
    run_command("$system_call");
    $nrm_long = $nseqs_init - determine_num_seqs_in_stk("$fam.tmp.stk");
    run_command("rm $fam.tmp.stk");
  }
  else { 
    $nrm_long = 0;
  }

  $$nrm_ambig_ref = $nrm_ambig;
  $$nrm_short_ref = $nrm_short;
  $$nrm_long_ref  = $nrm_long;
  $$nrm_total_ref = $nrm_total;
  return;
}


#####################################################################
# subroutine: consensusize_step
# incept:     EPN, Wed May 20 12:59:55 2009
# 
# purpose:    Create a consensus structure annotated stockholm
#             alignment.  The consensus structure is defined as the
#             set of basepairs between alignment columns i:j that are
#             basepaired in > x% of the individual structures.  'x' is
#             defined as the minimum possible value that creates a
#             'consistent' well-nested consensus
#             structure. 'consistent' is defined as a structure for
#             which no two basepairs i:j and k:l exist such that
#             either: (i==k and j!=l) or (i!=k and j==l). A
#             well-nested structure is one that has 0 pseudoknots, a
#             pseudoknot exists if two basepairs i:j and k:l if (i < k
#             < j < l).  purpose: Match CRW sequences we have bpseq
#             data for to their aligned versions in the aligned fasta
#             file using crw2ua_matches.pl (see that script for more
#             details).
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub consensusize_step { 
  $narg_expected = 3;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, match_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($indi_stk, $cons_stk, $construct_output) = @_;

  $system_call = "$esl_construct --fmin -o $cons_stk $indi_stk > $construct_output";
  $system_call2print = $system_call;
  $system_call2print =~ s/^\S+\///;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");

  # parse esl-construct output to print --fmin value
  open(TMP, $construct_output); 
  printf("construct_output: $construct_output\n");
  while($line = <TMP>) { 
    chomp $line;
    if($line =~ /^\s*(0\.\d+)\s+\d+\s+(\d+)/) { 
      $x = $1; $nbps = $2;
    }
  }
  close(TMP);
  
  # one final step, call esl-construct once more to get stats on the new cons struct,
  # and print them to the log file
  $output = `$esl_construct $cons_stk`;
  @lines = split(/\n/, $output);
  $nseq = 0;
  foreach $line (@lines) { 
    if($line =~ /^\s*SS\_cons \(consensus\)\s+(\d+)/) { $ncons = $1; }
    elsif($line =~ /\s*(\d+)\/\s*(\d+)\s+\((\S+)\)\s+conflict/) { 
      $nconflict = $1;
      $ntotal = $2;
      $fconflict = $3;
    }
    elsif($line =~ m/^\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+/) { 
      if($line !~ m/SS\_cons\(consensus\)/) { 
        $nseq++;
      }
    }
  }
  printf LOG ("  %4d  %.3f  %4d  %5d  %4d  %.3f", $nseq, $x, $nbps, $ntotal, $nconflict, $fconflict);
  return;
}

#####################################################################
# subroutine: fastaize_step
# incept:     EPN, Wed May 20 12:59:55 2009
# 
# purpose:    Create an unaligned fasta sequence file from a stockholm
#             alignment file.
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub fastaize_step { 
  $narg_expected = 2;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, fastaize_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_stk, $out_fa) = @_;

  $system_call = "$esl_reformat -o $out_fa fasta $in_stk";
  $system_call2print = $system_call;
  $system_call2print =~ s/^\S+\///;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");
  return;
}

#####################################################################
# subroutine: alnfastaize_step
# incept:     EPN, Wed May 20 14:17:50 2009
# 
# purpose:    Create an aligned fasta sequence file 
#             from a stockholm alignment file.
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub alnfastaize_step { 
  $narg_expected = 2;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, alnfastaize_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_stk, $out_fa) = @_;

  $system_call = "$esl_reformat -o $out_fa afa $in_stk";
  run_command("$system_call");

  return;
}

#####################################################################
# subroutine: cmbuild_step
# incept:     EPN, Wed May 20 13:03:15 2009
# 
# purpose:    Build a CM from a stockholm file.
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub cmbuild_step { 
  $narg_expected = 3;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, cmbuild_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_stk, $out_cm, $out_cmbuild) = @_;

  $system_call = "$cmbuild -F $cmbuild_options $out_cm $in_stk > $out_cmbuild";
  $system_call2print = $system_call;
  $system_call2print =~ s/^\S+\///;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");
  return;
}
#####################################################################
# subroutine: cmalign_step
# incept:     EPN, Wed May 20 13:04:43 2009
# 
# purpose:    Align sequences to a CM.
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub cmalign_step { 
  $narg_expected = 4;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, cmalign_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_cm, $in_fa, $out_stk, $out_cmalign) = @_;
  
  #align seqs to the CM
  $system_call = "$cmalign $cmalign_options -o $out_stk $in_cm $in_fa > $out_cmalign";
  $system_call2print = $system_call;
  $system_call2print =~ s/^\S+\///;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");
  return;
}


#####################################################################
# subroutine: rmconflicts_step
# incept:     EPN, Wed May 20 13:15:28 2009
# 
# purpose:    Remove sequences that fail our 'number of conflicts' 
#             quality requirements using esl-alimanip.
#             -  sequences must have <= <n> basepair conflicts between
#                their CRW individual structure annotation and the
#                implicit individual structure defined by imposing the
#                current consensus structure.
#
#             The default value of <n> is 15. It can be changed to
#             <n> with -c <n>.
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub rmconflicts_step { 
  $narg_expected = 3;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, prefilter_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_stk, $out_list, $out_stk) = @_;

  $system_call = "$esl_construct --lmax $maxconflicts -l $out_list $in_stk > /dev/null";
  $system_call2print = $system_call;
  $system_call2print =~ s/^\S+\///;
  printf("Running '$system_call2print'\n");
  run_command("$system_call");

  if(-s $out_list) { # only call esl-alimanip if the out list has >= 1 sequences in it
    $system_call = "$esl_alimanip -o $out_stk --seq-r $out_list $in_stk";
    $system_call2print = $system_call;
    $system_call2print =~ s/^\S+\///;
    printf("Running '$system_call2print'\n");
    run_command("$system_call");
    
    #determine number of seqs removed 
    $nremoved = 0;
    open(IN, $out_list); 
    while($line = <IN>) { 
      if($line =~ m/\w/) { $nremoved++; }
    }
    close(IN);
  }
  else { 
    $system_call = "cp $in_stk $out_stk";
    $system_call2print = $system_call;
    $system_call2print =~ s/^\S+\///;
    printf("Running '$system_call2print'\n");
    run_command("$system_call");
    $nremoved = 0;
  }
  return $nremoved;
}

#####################################################################
# subroutine: validate_input_files()
# incept:     EPN, Wed May 20 09:10:24 2009
# 
# purpose:    Validate that required files exist. Read the parameters
#             parameters file <$params_file> and verify that the
#             required variables defining locations of programs we
#             need to be able to run are defined, and that those
#             programs exist.
# 
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub validate_input_files { 
  my $narg_expected = 6;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, validate_input_files() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my($cmalign_ref, $cmbuild_ref, $esl_construct_ref, $esl_alimanip_ref, $esl_reformat_ref, $crw2ua_matches_ref) = @_;

  #make sure ARGV files exist
  if(!(-e $bpseq_file))    { printf STDERR ("ERROR, bpseq file $bpseq_file does not exist.\n"); exit(1); }
  if(!(-e $crwalnfa_file)) { printf STDERR ("ERROR, aligned fasta file $crwalnfa_file does not exist.\n"); exit(1); }
  if(!(-e $params_file))   { printf STDERR ("ERROR, parameter file $params_file does not exist.\n"); exit(1); }
  if($crwalnfa_file !~ m/\.alnfasta$/) { printf STDERR ("ERROR, CRW aligned fasta file name must end with \".alnfasta\", $crwalnfa_file does not.\n"); exit(1); }

  #Read parameter file
  $cmalign = "";
  $cmbuild = "";
  $esl_construct = "";
  $esl_alimanip = "";
  $esl_reformat = "";
  $crw2au_matches = "";
  require $params_file;
  # check to make sure required variables were set and are valid
  if($cmalign eq "")       { printf STDERR ("ERROR, \$cmalign not defined in parameter file.\n"); exit(1); }
  if($cmbuild  eq "")      { printf STDERR ("ERROR, \$cmbuild not defined in parameter file.\n"); exit(1); }
  if($esl_construct eq "") { printf STDERR ("ERROR, \$esl_construct not defined in parameter file.\n"); exit(1); }
  if($esl_alimanip eq "")  { printf STDERR ("ERROR, \$esl_alimanip not defined in parameter file.\n"); exit(1); }
  if($esl_reformat eq "")  { printf STDERR ("ERROR, \$esl_reformat not defined in parameter file.\n"); exit(1); }
  if($crw2ua_matches eq "") { printf STDERR ("ERROR, \$crw2ua_matches not defined in parameter file.\n"); exit(1); }

  # required variable paths can't have ~ prefix indicating home dir
  if($cmalign  =~ /^~/)     { printf STDERR ("ERROR, \$cmalign definition must be absolute path, \~ cannot be used as shortcut for home directory.\n"); exit(1); }
  if($cmbuild  =~ /^~/)     { printf STDERR ("ERROR, \$cmbuild definition must be absolute path, \~ cannot be used as shortcut for home directory.\n"); exit(1); }
  if($esl_construct =~ /^~/)    { printf STDERR ("ERROR, \$esl_construct definition must be absolute path, \~ cannot be used as shortcut for home directory.\n"); exit(1); }
  if($esl_alimanip  =~ /^~/)    { printf STDERR ("ERROR, \$esl_alimanip definition must be absolute path, \~ cannot be used as shortcut for home directory.\n"); exit(1); }
  if($esl_reformat  =~ /^~/)    { printf STDERR ("ERROR, \$esl_reformat definition must be absolute path, \~ cannot be used as shortcut for home directory.\n"); exit(1); }
  if($crw2ua_matches =~ /^~/)    { printf STDERR ("ERROR, \$crw2ua_matches definition must be absolute path, \~ cannot be used as shortcut for home directory.\n"); exit(1); }

  # required variable paths can't have ~ prefix indicating home dir
  if($cmalign  =~ /^\.\./)     { printf STDERR ("ERROR, \$cmalign definition must be absolute path, \.\. cannot be used as shortcut.\n"); exit(1); }
  if($cmbuild  =~ /^\.\./)     { printf STDERR ("ERROR, \$cmbuild definition must be absolute path, \.\. cannot be used as shortcut.\n"); exit(1); }
  if($esl_construct =~ /^\.\./)    { printf STDERR ("ERROR, \$esl_construct definition must be absolute path, \.\. cannot be used as shortcut.\n"); exit(1); }
  if($esl_alimanip  =~ /^\.\./)    { printf STDERR ("ERROR, \$esl_alimanip definition must be absolute path, \.\. cannot be used as shortcut.\n"); exit(1); }
  if($esl_reformat  =~ /^\.\./)    { printf STDERR ("ERROR, \$esl_reformat definition must be absolute path, \.\. cannot be used as shortcut.\n"); exit(1); }
  if($crw2ua_matches =~ /^\.\./)    { printf STDERR ("ERROR, \$crw2ua_matches definition must be absolute path, \.\. cannot be used as shortcut.\n"); exit(1); }

  # required variable paths can't have ~ prefix indicating home dir
  if($cmalign  =~ /^\./)     { printf STDERR ("ERROR, \$cmalign definition must be absolute path, \. cannot be used as shortcut.\n"); exit(1); }
  if($cmbuild  =~ /^\./)     { printf STDERR ("ERROR, \$cmbuild definition must be absolute path, \. cannot be used as shortcut.\n"); exit(1); }
  if($esl_construct =~ /^\./)    { printf STDERR ("ERROR, \$esl_construct definition must be absolute path, \. cannot be used as shortcut.\n"); exit(1); }
  if($esl_alimanip  =~ /^\./)    { printf STDERR ("ERROR, \$esl_alimanip definition must be absolute path, \. cannot be used as shortcut.\n"); exit(1); }
  if($esl_reformat  =~ /^\./)    { printf STDERR ("ERROR, \$esl_reformat definition must be absolute path, \. cannot be used as shortcut.\n"); exit(1); }
  if($crw2ua_matches =~ /^\./)    { printf STDERR ("ERROR, \$crw2ua_matches definition must be absolute path, \. cannot be used as shortcut.\n"); exit(1); }

  # required variable paths must exist
  if(!(-e $cmalign))        { printf STDERR ("ERROR, cmalign executable $cmsearch does not exist.\n"); exit(1); }
  if(!(-e $cmbuild))        { printf STDERR ("ERROR, cmalign executable $cmsearch does not exist.\n"); exit(1); }
  if(!(-e $esl_construct))       { printf STDERR ("ERROR, esl-construct executable $esl_construct does not exist.\n"); exit(1); }
  if(!(-e $esl_alimanip))        { printf STDERR ("ERROR, esl-alimanip executable $esl_alimanip does not exist.\n"); exit(1); }
  if(!(-e $esl_reformat))        { printf STDERR ("ERROR, esl-reformat executable $esl_reformat does not exist.\n"); exit(1); }
  if(!(-e $crw2ua_matches))      { printf STDERR ("ERROR, crw2ua_matches.pl script $crw2ua_matches does not exist.\n"); exit(1); }

  $$cmalign_ref = $cmalign;
  $$cmbuild_ref = $cmbuild;
  $$esl_construct_ref   = $esl_construct;
  $$esl_alimanip_ref   = $esl_alimanip;
  $$esl_reformat_ref   = $esl_reformat;
  $$crw2ua_matches_ref   = $crw2ua_matches;

  return;
}


#####################################################################
# subroutine: create_output_dir()
# incept:     EPN, Mon Nov  3 15:22:26 2008
# 
# purpose:    Create a new directory <$outdir> for the output files. 
#             If the directory already exists and <$do_clobber>, 
#             delete all of the files within it.
# 
# returns:    Nothing.
#
# exits:      If the directory already exists and <$do_clobber> is
#             FALSE. If a system call unexpectedly fails and returns 
#             a non-zero status code. 
# 
####################################################################
sub create_output_dir { 
  my $narg_expected = 2;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, create_output_dir() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my($do_clobber, $out_dir) = @_;

  # create output directory
  if(-d $out_dir) { 
    if(!($do_clobber)) { 
      printf STDERR ("ERROR, output directory $out_dir already exists. Delete it or use -F to overwrite it.\n"); exit(1); 
    }
    else { # dir exists, but -F enabled, so we remove it
      if($out_dir eq "") { printf STDERR ("ERROR, trying to create directory named \"\"\n"); exit(1); }
      run_command("rm -rf $out_dir/*");
      run_command("rmdir $out_dir");
      if(-d $out_dir) { printf STDERR ("ERROR, output directory $out_dir still exists after a system\(\"rmdir $out_dir\"\) call.\n"); exit(1); }
    }
  }
  # if we get here, either $out_dir does not yet exist, or it does but -F was set on command line
  run_command("mkdir $out_dir");
  #run_command("cp $bpseq_file $out_dir\/");
  return;
}


#####################################################################
# subroutine: determine_num_seqs_in_stk
# incept:     EPN, Wed May 20 13:28:13 2009
# 
# purpose:    Determine number of sequences in a stk alignment.
#
# returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub determine_num_seqs_in_stk {
  $narg_expected = 1;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, match_step() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($in_stk) = $_[0];

  $system_call = "$esl_alistat --rna --list $fam.tmp.list $in_stk > /dev/null";
  run_command("$system_call");
  my $nseqs = 0; 
  open(IN, "$fam.tmp.list"); 
  while($line = <IN>) { 
    $nseqs++; 
  } 
  close(IN);
  run_command("rm $fam.tmp.list");
  return $nseqs;
}


#####################################################################
# subroutine: run_command
# incept:     EPN, Tue Jan 24 09:51:51 2017
# 
# purpose:    Run a command using 'system' and die if it fails.
#
# returns:    Nothing
# 
####################################################################
sub run_command {
  $narg_expected = 1;
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, run_command() entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my ($cmd) = $_[0];
  system($cmd);

  if($? != 0) {
    die "ERROR in run_command, the following command failed:\n$cmd\n";
  }

  return;
}
