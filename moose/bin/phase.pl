#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

my $pedigree_test_filename = shift; # file with:   acc_id  bad_count  mat_id  pat_id  ...

my $n_sites_to_print_on_line = shift // 0;

my $max_bad = shift // 6;

open my $fhp, "<", "$pedigree_test_filename";
my %parentalidpair_offspringidstr = ();
while (my $line = <$fhp>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($accid, $matid, $patid) = @cols[0,2,3];
  my $parentalidpair = ($matid gt $patid)? "$matid $patid" : "$patid $matid";
  $parentalidpair_offspringidstr{$parentalidpair} .= "$accid ";
}

my @marker_ids;
my %id_genotypes = ();
while (my $line = <>) {	# out_genotypes file; each line like:  acc_id  0010120001201310...
  next if($line =~ /^\s*#/);
  if ($line =~ /^MARKER/) {
    @marker_ids = split(" ", $line);
    shift @marker_ids;		# get rid of initial 'MARKER'
  } elsif ($line =~ /^\s*(\S+)\s+(\S+)\s*$/) {
    my $id = $1;
    my $genotypes = $2;
    $id_genotypes{$id} = $genotypes;
  }
}

while (my($paridpair, $progstr) = each %parentalidpair_offspringidstr) {
  my ($p1id, $p2id) = split(" ", $paridpair);
  my $p1gts = $id_genotypes{$p1id};
  my $p2gts = $id_genotypes{$p2id};

  my @offspringids = split(" ", $progstr);
  my @prog_0counts = ();
  my @prog_1counts = ();
  my @prog_2counts = ();
  my @prog_mdcounts = ();
  for my $offid (@offspringids) {
    my $offgts = $id_genotypes{$offid};
    for (my $i = 0; $i < length $offgts; $i++) {
      my $offgt = substr($offgts, $i, 1);
      if ($offgt eq '0') {
	$prog_0counts[$i]++;
      } elsif ($offgt eq '1') {
	$prog_1counts[$i]++;
      } elsif ($offgt eq '2') {
	$prog_2counts[$i]++;
      } elsif ($offgt eq '3') {	# missing data
	$prog_mdcounts[$i]++;
      } else {
	die "gt is $offgt (should be 0,1,2, or 3)\n";
      }
    }
  }

  my @p1heterozyg_indices = ();
  my @p2subset_gts = ();
  for (my $i=0; $i<length $p1gts; $i++) {
    if (substr($p1gts, $i, 1) eq '1') {
      my $p2gt = substr($p2gts, $i, 1);
      if (
	  ($p2gt eq '0'  and  ($prog_2counts[$i] // 0) < $max_bad)
	  or
	  ($p2gt eq '2'  and  ($prog_0counts[$i] // 0) < $max_bad)
	 ) {
	push @p1heterozyg_indices, $i;
	push @p2subset_gts, $p2gt;
      }
    }
  }

 
  my @progsubset_0counts = ();
  my @progsubset_1counts = ();
  my @progsubset_2counts = ();
  my @progsubset_mdcounts = ();
  my @cis_counts = ();
  my @trans_counts = ();
  my @md_counts = ();
  my @impossible_counts = ();

  my @cis2_counts = ();
  my @trans2_counts = ();
  my @md2_counts = ();
  my @impossible2_counts = ();

  

  my $n_p1het_p2hom = scalar @p1heterozyg_indices;

  for my $offid (@offspringids) {
    my $offgts = $id_genotypes{$offid};
    for (my $i=0; $i < scalar @p1heterozyg_indices; $i++) {

      my $j1 = $p1heterozyg_indices[$i];
      my $p2gt1 = $p2subset_gts[$i];
      my $offgt1 = substr($offgts, $j1, 1);

      if ($offgt1 eq '0') {
	$progsubset_0counts[$i]++;
      } elsif ($offgt1 eq '1') {
	$progsubset_1counts[$i]++;
      } elsif ($offgt1 eq '2') {
	$progsubset_2counts[$i]++;
      } else {			# missing data
	$progsubset_mdcounts[$i]++;
      }
      if ($i < ($n_p1het_p2hom - 1)) {
	my $j2 = $p1heterozyg_indices[$i+1];
	my $p2gt2 = $p2subset_gts[$i+1];
        my $offgt2 = substr($offgts, $j2, 1);
	increment_cis_trans_etc($i, $p2gt1, $p2gt2, $offgt1, $offgt2, \@cis_counts, \@trans_counts, \@impossible_counts, \@md_counts);
	if ($i < ($n_p1het_p2hom - 2)) {
	  my $j3 =  $p1heterozyg_indices[$i+2];
	  my $p2gt3 = $p2subset_gts[$i+2];
	  my $offgt3 = substr($offgts, $j3, 1);
   
	  increment_cis_trans_etc($i+1, $p2gt1, $p2gt3, $offgt1, $offgt3, \@cis2_counts, \@trans2_counts, \@impossible2_counts, \@md2_counts);
	}
      }
    }
  }

  if ($n_sites_to_print_on_line > 0) {
    print("idxs:  ");
    for my $i (@p1heterozyg_indices[0..$n_sites_to_print_on_line]) {
      printf("%3i ", $i);
    }
    print("\n");
    print("p2gts: ");
    for my $p2gt (@p2subset_gts[0..$n_sites_to_print_on_line]) {
      printf("%3i ", $p2gt);
    }
    print("\n");

    print("prog0  ");
    for my $prog0count (@progsubset_0counts[0..$n_sites_to_print_on_line]) {
      printf("%3i ", $prog0count // 0);
    }
    print("\n");
    print("prog1  ");
    for my $prog1count (@progsubset_1counts[0..$n_sites_to_print_on_line]) {
      printf("%3i ", $prog1count // 0);
    }
    print("\n");
    print("prog2  ");
    for my $prog2count (@progsubset_2counts[0..$n_sites_to_print_on_line]) {
      printf("%3i ", $prog2count // 0);
    }
    print("\n");
    print("progmd ");
    for my $progmdcount (@progsubset_mdcounts[0..$n_sites_to_print_on_line]) {
      printf("%3i ", $progmdcount // 0);
    }
    print("\n");
  

    print("cis      ");
    for my $cis_count (@cis_counts[0..($n_sites_to_print_on_line-1)]) {
      printf("%3i ", $cis_count // 0);
    }
    print("\n");

    print("trans    ");
    for my $trans_count (@trans_counts[0..($n_sites_to_print_on_line-1)]) {
      printf("%3i ", $trans_count // 0);
    }
    print("\n");

    print("imp      ");
    for my $impossible_count (@impossible_counts[0..($n_sites_to_print_on_line-1)]) {
      printf("%3i ", $impossible_count // 0);
    }
    print("\n");

    print("md       ");
    for my $md_count (@md_counts[0..($n_sites_to_print_on_line-1)]) {
      printf("%3i ", $md_count // 0);
    }
    print("\n");

  } else {			# 1 line for each marker used
    my $prev_marker_position_number = 0;
    while (my ($i, $j) = each @p1heterozyg_indices) {
      my $marker_id = $marker_ids[$j];
      my $marker_position_number = ($marker_id =~ /[:](\d+)/)? $1 : 0;
      
      printf("%-16s %8ld %4i %3i    ", $marker_id,
	     $marker_position_number - $prev_marker_position_number, $j, $p2subset_gts[$i]);
      if ($p2subset_gts[$i] == '0') {
	printf("%3i %3i %2i %2i   ",
	       $progsubset_0counts[$i] // 0, $progsubset_1counts[$i] // 0,
	       $progsubset_2counts[$i] // 0, $progsubset_mdcounts[$i] // 0);
      } else {
	printf("%3i %3i %2i %2i   ",
	       $progsubset_2counts[$i] // 0, $progsubset_1counts[$i] // 0,
	       $progsubset_0counts[$i] // 0, $progsubset_mdcounts[$i] // 0);
      }

      my ($cis1_count, $trans1_count) = ($cis_counts[$i] // 0, $trans_counts[$i] // 0);
      my ($cis2_count, $trans2_count) = ($cis2_counts[$i] // 0, $trans2_counts[$i] // 0);
      printf("%3i %3i %2i %2i %3i   %3i %3i %2i %2i %3i  \n",
	     $cis1_count, $trans1_count, 
	     $impossible_counts[$i] // 0, $md_counts[$i] // 0,
	     min($cis1_count, $trans1_count),
	     
	     $cis2_count, $trans2_count, 
	     $impossible2_counts[$i] // 0, $md2_counts[$i] // 0,
	       min($cis2_count, $trans2_count) );
      $prev_marker_position_number = $marker_position_number;
    }
    print("\n");

  }
}

sub increment_cis_trans_etc{
  my $i = shift;
  my $p2gt1 = shift;
  my $p2gt2 = shift;
  my $offgt1 = shift;
  my $offgt2 = shift;
  my $cis_counts = shift;
  my $trans_counts = shift;
  my $impossible_counts = shift;
  my $md_counts = shift;

  my $p2gtpair = "$p2gt1 $p2gt2";
  my $offgtpair = "$offgt1 $offgt2";
  if ($p2gtpair eq "2 0") {
    $p2gtpair = "0 2";
    $offgtpair = "$offgt2 $offgt1";
  }

  if ($offgt1 eq '3'  or  $offgt2 eq '3') {
    $md_counts->[$i]++;
  } elsif ($p2gtpair eq "0 0") {
    if ($offgtpair eq "0 0"  or $offgtpair eq "1 1") {
      $cis_counts->[$i]++;
    } elsif ($offgtpair eq "0 1"  or  $offgtpair eq "1 0") {
      $trans_counts->[$i]++;
    } else {
      $impossible_counts->[$i]++;
    }
  } elsif ($p2gtpair eq "0 2") {
    if ($offgtpair eq "0 1"  or  $offgtpair eq "1 2") {
      $cis_counts->[$i]++;
    } elsif ($offgtpair eq "0 2"  or  $offgtpair eq "1 1") {
      $trans_counts->[$i]++;
    } else {
      $impossible_counts->[$i]++;
    }
  } elsif ($p2gtpair eq "2 2") {
    if ($offgtpair eq "1 1"  or  $offgtpair eq "2 2") {
      $cis_counts->[$i]++;
    } elsif ($offgtpair eq "1 2"  or  $offgtpair eq "2 1") {
      $trans_counts->[$i]++;
    } else {
      $impossible_counts->[$i]++;
    }
  } else {
    die "p2gtpair: $p2gtpair (shouldn't occur)\n";
  }
}
