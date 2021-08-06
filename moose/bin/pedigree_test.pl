#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
# use Graphics::GnuplotIF qw(GnuplotIF);
# use Math::GSL::SF  qw( :all );

use File::Basename 'dirname';
use Time::HiRes qw( gettimeofday );
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;
print STDERR "# libdir: $libdir \n";


use GenotypesSet;
use PedigreesDAG;
use CheckPedigrees;
use Cluster1d;

# *****  run C program to get agmr, hgmr, r for parent-offspring pairs according to pedigree file,
# *****  the further analysis with perl.

my $gtfilename = undef;
my $pedigree_table_filename = undef;
my $delta = 0.05; # real-number genotypes are rounded to nearest integer (0,1,2) if within +- $delta
my $max_bad_gt_fraction = 1.0;
# my $n_random_parents = 0;
my $base_output_filename = 'out';
my $output_pedigrees = 0;
my $output_genotype_matrix = 0;
my $checkpedigrees_progress_interval = 500;
my $genotypeset_progress_interval = 500;
my $pow = 'log';
my $find_alternatives = 0;

my $c_output_filename = 'c_out.tmp';
GetOptions(
	   'gtsfile|gtfile|genotypesfile=s' => \$gtfilename,
	   'pedigreein|pedigreefile|pedtable=s' => \$pedigree_table_filename,
	   'delta=f' => \$delta,
	   'max_bad_gt_fraction=f' => \$max_bad_gt_fraction,
	   'baseout|basename=s' => \$base_output_filename,
	   'outpedigrees!' => \$output_pedigrees,
	   'outmatrix!' => \$output_genotype_matrix,
	   'pow=s' => \$pow,
	   'alternatives=i' => \$find_alternatives, # 0: none, 1: only for 'bad' pedigrees, >=2: do for all pedigrees.
	  );

die "No genotypes matrix filename provided.\n" if(!defined $gtfilename);

# Read in the pedigree table:

my $t_start = gettimeofday();
print STDERR "# Creating PedigreesDAG object from file: $pedigree_table_filename\n";
my $pedigrees = PedigreesDAG->new({pedigree_filename => $pedigree_table_filename});
print STDERR "# PedigreesDAG object created.\n\n";
if ($output_pedigrees) {
  my $pedigree_output_filename = $base_output_filename . '_pedigrees';
  open my $fhout, ">", "$pedigree_output_filename";
  print $fhout $pedigrees->as_string();
  close $fhout;
}
print STDERR "# Time to read pedigree file, create PedigreesDAG obj: ", gettimeofday() - $t_start, "\n";

# Read in the genotype matrix file
my $wd = `pwd`;
print "wd: $wd\n";
my $command = "~/simgtsrch/src/pedigree_test -g $gtfilename -p $pedigree_table_filename  -w $delta  -x $max_bad_gt_fraction  -o $c_output_filename ";
print STDERR "command: $command \n";
system "$command";

# exit;
open my $fhin, "<", "$c_output_filename" or die "Couldn't open $c_output_filename for reading.\n";

my @lines = <$fhin>;

print scalar @lines, "\n";

my @matpat_agmrs = ();
#my @array_of_lines_as_cols = ();
while (my ($j, $line) = each @lines) {
  my @cols = split(" ", $line);
  #   push @array_of_lines_as_cols, \@cols;
  push @matpat_agmrs, $cols[5];
}
my $cluster1d_obj = Cluster1d->new({label => 'agmr between parents', xs => \@matpat_agmrs, pow => $pow});
my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
#print "clustering of hgmr (female parent): $n_pts k-means: $km_n_L below $km_n_R above $km_h_opt;  kde: $kde_n_L below $kde_n_R above $kde_h_opt.\n";
printf("clustering of agmr between parents: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);

my $agmr_h = $km_h_opt;

my @hgmr_denoms = ();
my @hgmrs = ();
my @r_denoms = ();
my @rs = ();
my @d1_denoms = ();
my @d1s = ();



# my @self_hgmrs = ();		# self (according to pedigree)
# my @bip_hgmrs = ();		# biparental (according to pedigree)
# my @self_rs = ();
# my @bip_rs = ();
# my @self_d1s = ();
# my @bip_d1s = ();


# #####  clustering  ##########################
while (my ($j, $line) = each @lines) {
  my @cols = split(" ", $line);

  my ($accid, $bad_gt_count) = @cols[0,1];
  my ($mat_id, $pat_id) = @cols[2,3];
  
  my $matpat_agmr = $cols[5];
  my ($mat_hgmr, $pat_hgmr) = @cols[7,11];
  my ($mat_r, $pat_r) = @cols[9,13];
  my ($d1, $d2) = @cols[15,17];
  
  push @hgmr_denoms, @cols[6,10];
  push @hgmrs, ($mat_hgmr, $pat_hgmr);
  push @r_denoms, @cols[8,12];
  push @rs, ($mat_r, $pat_r);
  push @d1_denoms, $cols[14];
  push @d1s, $d1;
}
my $median_matpat_agmr_denom = $matpat_agmrs[int(scalar @matpat_agmrs / 2)];
my $median_hgmr_denom = $hgmr_denoms[int(scalar @hgmr_denoms / 2)];
my $median_r_denom = $r_denoms[int(scalar @r_denoms / 2)];
my $median_d1_denom = $d1_denoms[int(scalar @d1_denoms / 2)];

$cluster1d_obj = Cluster1d->new({label => 'hgmr', xs => \@hgmrs, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("clustering of hgmr: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_hgmr = $km_h_opt;

$cluster1d_obj = Cluster1d->new({label => 'r', xs => \@rs, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("clustering of    r: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_self_r = $km_h_opt;

$cluster1d_obj = Cluster1d->new({label => 'd1', xs => \@d1s, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("clustering of   d1: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_d1 = $km_h_opt;



# if($hgmr_denoms_ok  and  $r_denoms_ok  and  $d1_denoms_ok){
#   if($Fhgmr <= $max_ok_hgmr  and  $Mhgmr <= $max_ok_hgmr  and  $d1 <= $max_ok_d1){
#   }
# }else{ # high missing data

# }

if ($find_alternatives > 0) {

  $command = "~/simgtsrch/src/pedigree_test -g $gtfilename -p $pedigree_table_filename  -w $delta  -x $max_bad_gt_fraction -A $find_alternatives -o pedigrees_w_alternatives ";
  $command .= " -a $agmr_h  -h $max_ok_hgmr -r $max_self_r -D $max_ok_d1";
  print STDERR "command: $command \n";
  system "$command";

}

open my $fh_alt, "<", "pedigrees_w_alternatives";
@lines = <$fh_alt>;

my $factor = 0.33;
my %category_counts = ();
my $bad_denoms_count = 0;
while (my ($j, $line) = each @lines) {
  my @cols = split(" ", $line);
  my ($acc_id, $acc_md) = @cols[0,1];
  my $ped_d1 = $cols[15];
  my $category_string = '';
  my $denoms_ok = are_denoms_ok(\@cols, 4, $median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_d1_denom, $factor);

  if($denoms_ok){
    $category_string = category(\@cols, 4, $agmr_h, $max_ok_hgmr, $max_self_r, $max_ok_d1);
    $category_counts{$category_string}++;
  }else{
    $category_string = 'bad_denoms';
    $bad_denoms_count++;
  }

    # alternative pedigrees:
  my $n_alternatives = $cols[18] // 0; #
  my $good_alternative_count = 0;
  my %alt_d1 = ();
  for (my $i = 0; $i < $n_alternatives; $i++) {
    my $first = 21 + $i*16;
    $denoms_ok = are_denoms_ok(\@cols, $first, $median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_d1_denom, 0.33);
    if($denoms_ok){
      my $alt_id_pair = $cols[$first-2] . ' ' . $cols[$first-1];
      my $alt_d1 = $cols[$first+11];
      my $alt_category_string = category(\@cols, $first, $agmr_h, $max_ok_hgmr, $max_self_r, $max_ok_d1);
      if($alt_category_string eq '0 00 00 0'  or $alt_category_string eq '1 01 01 0'){
	$good_alternative_count++;
	$alt_d1{$alt_id_pair . "  " . $alt_category_string} = $alt_d1;
      }
    }
  }
  print $cols[0];
  if($category_string eq '0 00 00 0'  or  $category_string eq '1 01 01 0'){
    print "  ped  ", $cols[2], "  ", $cols[3], "   $category_string  $ped_d1  $good_alternative_count\n";
  }elsif($good_alternative_count > 0){
    my @sorted_alts = sort {$alt_d1{$a} <=> $alt_d1{$b} } keys %alt_d1;
    my $best_alt = $sorted_alts[0];
    print "  alt  ", $best_alt, "  ", $alt_d1{$best_alt}, "  $good_alternative_count \n";
  }else{
    print "  none \n";
  }
}
print "bad denoms count: $bad_denoms_count\n";
while (my($cat, $count) = each %category_counts) {
  print "$cat $count\n";
}



exit;


sub xxx{
  my @cols = @{my $cls = shift};
  my $Fparent_id = shift @cols;
  my $Mparent_id = shift @cols;
  my $agmr_denom = shift @cols;
  my $agmr = shift @cols;
  my $Fhgmr_denom = shift @cols;
  my $Fhgmr = shift @cols;
  my $Fr_denom = shift @cols;
  my $Fr = shift @cols;
  my $Mr_denom = shift @cols;
  my $Mr = shift @cols;
  my $d1_denom = shift @cols;
  my $d1 = shift @cols;
  my $d2_denom = shift @cols;
  my $d2 = shift @cols;
 
}

sub are_denoms_ok{
  my @cols = @{my $cls = shift};
  my $first = shift;
  my ($median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_d1_denom) = @_;
  my $last = $first + 11;
  my ($FMagmr_denom, $FMagmr, $Fhgmr_denom, $Fhgmr, $Fr_denom, $Fr, $Mhgmr_denom, $Mhgmr, $Mr_denom, $Mr, $d1_denom, $d1) = @cols[$first..$last];
  my $agmr_denom_ok = ($FMagmr_denom >= $factor*$median_matpat_agmr_denom);
  my $hgmr_denoms_ok = ($Fhgmr_denom >= $factor*$median_hgmr_denom  and  $Mhgmr_denom >= $factor*$median_hgmr_denom);
  my $r_denoms_ok = ($Fr_denom >= $factor*$median_r_denom  and  $Mr_denom >= $factor*$median_r_denom);
  my $d1_denom_ok = ($d1_denom >= $factor*$median_d1_denom);
  my $denoms_ok = ($agmr_denom_ok  and  $hgmr_denoms_ok  and  $r_denoms_ok  and  $d1_denom_ok);
  return $denoms_ok;
}

sub category{
  my @cols = @{my $cls = shift};
  my $first = shift;
  my ($agmr_h, $max_ok_hgmr, $max_self_r, $max_ok_d1) = @_;
  my $last = $first + 11;
  my ($FMagmr_denom, $FMagmr, $Fhgmr_denom, $Fhgmr, $Fr_denom, $Fr, $Mhgmr_denom, $Mhgmr, $Mr_denom, $Mr, $d1_denom, $d1) = @cols[$first..$last];
  # my $agmr_denom_ok = ($FMagmr_denom >= $factor*$median_matpat_agmr_denom);
  # my $hgmr_denoms_ok = ($Fhgmr_denom >= $factor*$median_hgmr_denom  and  $Mhgmr_denom >= $factor*$median_hgmr_denom);
  # my $r_denoms_ok = ($Fr_denom >= $factor*$median_r_denom  and  $Mr_denom >= $factor*$median_r_denom);
  # my $d1_denom_ok = ($d1_denom >= $factor*$median_d1_denom);
  # my $denoms_ok = ($agmr_denom_ok  and  $hgmr_denoms_ok  and  $r_denoms_ok  and  $d1_denom_ok);
  my $category_string = '';
  $category_string .= ($FMagmr <= $agmr_h)? '0' : '1';
  $category_string .= ' ';
  $category_string .= ($Fhgmr <= $max_ok_hgmr)? '0' : '1';
  $category_string .= ($Fr <= $max_self_r)? '0' : '1';
  $category_string .= ' ';
  $category_string .= ($Mhgmr <= $max_ok_hgmr)? '0' : '1';
  $category_string .= ($Mr <= $max_self_r)? '0' : '1';
  $category_string .= ' ';
  $category_string .= ($d1 <= $max_ok_d1)? '0' : '1';
  return $category_string;
}
