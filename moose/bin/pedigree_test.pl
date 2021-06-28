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
my $n_random_parents = 0;
my $base_output_filename = 'out';
my $output_pedigrees = 0;
my $output_genotype_matrix = 0;
my $output_27triple_counts = 0;
my $output_14_counts = 0;
my $checkpedigrees_progress_interval = 500;
my $genotypeset_progress_interval = 500;

my $c_output_filename = 'c_out.tmp';
 GetOptions(
	    'gtsfile|gtfile|genotypesfile=s' => \$gtfilename,
	    'pedigreein|pedigreefile|pedtable=s' => \$pedigree_table_filename,
	    'delta=f' => \$delta,
	    'max_bad_gt_fraction=f' => \$max_bad_gt_fraction,
	    'n_random_parents=i' => \$n_random_parents,
	    'baseout|basename=s' => \$base_output_filename,
	    'outpedigrees!' => \$output_pedigrees,
	    'outmatrix!' => \$output_genotype_matrix,
	    'out27!' => \$output_27triple_counts,
	    'out14!' => \$output_14_counts,

	   );

die "No genotypes matrix filename provided.\n" if(!defined $gtfilename);

# Read in the pedigree table:

my $t_start = gettimeofday();
print STDERR "# Creating PedigreesDAG object from file: $pedigree_table_filename\n";
my $pedigrees = PedigreesDAG->new({pedigree_filename => $pedigree_table_filename});
print STDERR "# PedigreesDAG object created.\n\n";
if($output_pedigrees){
  my $pedigree_output_filename = $base_output_filename . '_pedigrees';
  open my $fhout, ">", "$pedigree_output_filename";
  print $fhout $pedigrees->as_string();
  close $fhout;
}
print STDERR "# Time to read pedigree file, create PedigreesDAG obj: ", gettimeofday() - $t_start, "\n";

# Read in the genotype matrix file

system "../../pedigree_test -g $gtfilename -p $pedigree_table_filename  -d $delta  -x $max_bad_gt_fraction  -o $c_output_filename ";

open my $fhin, "<", "$c_output_filename" or die "Couldn't open $c_output_filename for reading.\n";

my @mat_hgmrs = ();
my @mat_rs = ();
my @pat_hgmrs = ();
my @pat_rs = ();
my @matpat_agmrs = ();
while(my $line = <$fhin>){
  my @cols = split(" ", $line);
  my ($accid, $bad_gt_count, $mat_id, $pat_id, $mat_hgmr, $mat_r, $pat_hgmr, $pat_r, $matpat_agmr) =
    @cols[0,1,2,3,7,9,13,15,17];
  push @mat_hgmrs, $mat_hgmr;
  push @mat_rs, $mat_r;
  push @pat_hgmrs, $pat_hgmr;
  push @pat_rs, $pat_r;
  push @matpat_agmrs, $matpat_agmr;
}


 my $cluster1d_obj = Cluster1d->new({label => 'hgmr(female parent)', xs => \@mat_hgmrs});
  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
print "clustering of hgmr (female parent): $n_pts k-means: $km_n_L below, and $km_n_R above $km_h_opt; kde: $kde_n_L below and $kde_n_R above $kde_h_opt.\n";

$cluster1d_obj = Cluster1d->new({label => 'hgmr(male parent)', xs => \@pat_hgmrs});
  ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
print "clustering of hgmr (male parent): $n_pts k-means: $km_n_L below, and $km_n_R above $km_h_opt; kde: $kde_n_L below and $kde_n_R above $kde_h_opt.\n";

$cluster1d_obj = Cluster1d->new({label => 'r(female parent)', xs => \@mat_rs});
  ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
print "clustering of r (female parent): $n_pts k-means: $km_n_L below, and $km_n_R above $km_h_opt; kde: $kde_n_L below and $kde_n_R above $kde_h_opt.\n";

$cluster1d_obj = Cluster1d->new({label => 'r(male parent)', xs => \@pat_rs});
  ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
print "clustering of r (male parent): $n_pts k-means: $km_n_L below, and $km_n_R above $km_h_opt; kde: $kde_n_L below and $kde_n_R above $kde_h_opt.\n";

$cluster1d_obj = Cluster1d->new({label => 'agmr(between parents)', xs => \@matpat_agmrs});
  ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
 print "clustering of agmr (between parents): $n_pts k-means: $km_n_L below, and $km_n_R above $km_h_opt; kde: $kde_n_L below and $kde_n_R above $kde_h_opt.\n";



exit;

$t_start = gettimeofday();
print STDERR "# Creating GenotypesSet object.\n# Reading in the gts matrix file.\n";
my $gtset = GenotypesSet->new({gt_matrix_filename => $gtfilename, delta => $delta, max_bad_gt_fraction => $max_bad_gt_fraction, progress_report_interval => $genotypeset_progress_interval});
print STDERR "# GenotypesSet object created.\n\n";
print STDERR "# Time to read in genotypes file, create GenotypesSet obj: ", gettimeofday() - $t_start, "\n";
# while(my ($k, $v) = each %{$gtset->accid_genotypes()}){
# print "xxx: [", length $v->genotypes(), "]\n";
# }

# if( ($max_bad_gt_fraction < 1.0) or ($min_hw_qual_param > 0.0) ){
#   print STDERR "# removing bad markers. max bad fraction: $max_bad_gt_fraction   min hw qual: $min_hw_qual_param \n";
#   $gtset->clean_marker_set($max_bad_gt_fraction, $min_hw_qual_param);
# }

if($output_genotype_matrix){
  print STDERR "Writing output genotype matrix.\n";
  my $gtmatrix_output_filename = $base_output_filename . '_genotypeset';
  open my $fhout, ">", "$gtmatrix_output_filename";
  print $fhout  $gtset->as_string();
  close $fhout;
}

$t_start = gettimeofday();
print STDERR "# Create CheckPedigrees object.\n";
my $summary_output_filename = $base_output_filename . "_summary";
my $output_filename_27triple_counts = ($output_27triple_counts)? $base_output_filename . "_27" : undef;
my $output_filename_14counts = ($output_14_counts)? $base_output_filename . "_14" : undef;
my $check_peds = CheckPedigrees->new({
				      genotypes_set => $gtset,
				      pedigrees => $pedigrees,
				      n_random_parents => $n_random_parents,
				      summary_filename => $summary_output_filename,
				      triplecounts27_filename => $output_filename_27triple_counts,
				      counts14_filename => $output_filename_14counts,
				      progress_report_interval => $checkpedigrees_progress_interval,
				     });
print STDERR "# PedigreeChecks object created.\n";
print STDERR "# Time to create CheckPedigrees obj: ", gettimeofday() - $t_start, "\n";

$t_start = gettimeofday();
printf("# Female parent - offspring hgmr clusters: n pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $check_peds->m_hgmr_cluster());
printf("#   Male parent - offspring hgmr clusters: n_pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $check_peds->p_hgmr_cluster());

$check_peds->categorize();

print STDERR "# Time of cluster and categorize: ", gettimeofday() - $t_start, "\n";
