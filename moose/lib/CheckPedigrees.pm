package CheckPedigrees;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );
use Genotypes;
use GenotypesSet;
#use Pedigrees;
use PedigreeCheck;

# accession id and string of genotypes

has genotypes_set => (
		      isa => 'Object',
		      is => 'ro',
		      required => 1,
		     );

has pedigrees => (
		  isa =>  'Object',
		  is => 'ro',
		  required => 1,
		 );

has n_random_parents => (
			 isa => 'Int',
			 is => 'ro',
			 default => 0,
			);

has summary_filename => (
			 isa => 'Str',
			 is => 'ro',
			 default => 'summary',
			);

has triplecounts27_filename => (
				isa => 'Maybe[Str]',
				is => 'ro',
				default => undef,
			       );
has counts14_filename => (
			  isa => 'Maybe[Str]',
			  is => 'ro',
			  default => undef,
			 );

has pedigree_checks => ( # keys: ids, values: pedigree_check objects.
			isa => 'HashRef',
			is => 'rw',
			default => sub { {} },
		       );

has matrand_checks => (
		       isa => 'HashRef',
		       is => 'rw',
		       default => sub { {} },
		      );

has patrand_checks => (
		       isa => 'HashRef',
		       is => 'rw',
		       default => sub { {} },
		      );

has randrand_checks => (
			isa => 'HashRef',
			is => 'rw',
			default => sub { {} },
		       );

has progress_report_interval => (
				 isa => 'Int',
				 is => 'ro',
				 default => 500,
				);

has Fhgmr_h => ( # dividing value between the two clusters
		isa => 'Maybe[Num]',
		is => 'rw',
		default => undef,
	       );

has Mhgmr_h => ( # dividing value between the two clusters
		isa => 'Maybe[Num]',
		is => 'rw',
		default => undef,
	       );

has FMagmr_h => ( # dividing value between the two clusters
		isa => 'Maybe[Num]',
		is => 'rw',
		default => undef,
	       );

has rx01_h => ( # dividing value between the two clusters
		isa => 'Maybe[Num]',
		is => 'rw',
		default => undef,
	       );

has r0x1_h => ( # dividing value between the two clusters
		isa => 'Maybe[Num]',
		is => 'rw',
		default => undef,
	      );

has id_category => ( # keys: ids, values: category strings (e.g. '01010')
			  isa => 'HashRef[Str]',
			  is => 'rw',
			  default => sub{ {} },
			 );

has category_ids => ( # keys: category strings (e.g. '01101'), values: space separated strings of ids.
		     isa => 'HashRef[Str]',
		     is => 'rw',
		     default => sub{ {} },
		     );

sub BUILD{
  my $self = shift;

  my $gtsetobj = $self->genotypes_set();
  my $accid_gts = $gtsetobj->accid_genotypes(); # values are Genotypes objects
  my $peds = $self->pedigrees();
  #my $accid_parents = $peds->accid_parents();

  my $ostr = $gtsetobj->summary_info_string();

  print STDERR $ostr;
  open my $fhsummary, ">", $self->summary_filename();
  print $fhsummary $ostr;
  my $fh27;
  if (defined $self->triplecounts27_filename()) {
    open $fh27, ">", $self->triplecounts27_filename();
    print $fh27 $ostr;
  }
  my $fh14;
  if (defined $self->counts14_filename()) {
    open $fh14, ">", $self->counts14_filename();
    print $fh14 $ostr;
  }

  my $n_pedigree_checks = 0;
  my @accession_ids = sort keys %$accid_gts;
#  print STDERR "Number of accessions in CheckPedigrees BUILD: ", scalar @accession_ids, "\n";
  for my $accid (@accession_ids) {
    my $accgtobj = $accid_gts->{$accid};
    my ($matid, $patid) = $peds->parents($accid); # [matid, patid] if defined
  #  print STDERR $matid // 'undef', "  ", $patid // 'undef', "\n";
    next if(!defined $matid  or  !defined $patid);

    my $matgtobj = $accid_gts->{$matid} // undef;
 #   print STDERR "$matid  $patid   ", (defined $matgtobj)? $matgtobj->id() : 'matgtobj undefined' , "\n";
    next if (!defined $matgtobj);
    my $patgtobj = $accid_gts->{$patid} // undef;
    next if (!defined $patgtobj);
    my $pedcheck = PedigreeCheck->new({ mat_gtsobj => $matgtobj, pat_gtsobj => $patgtobj, acc_gtsobj => $accgtobj, n_random_parents => $self->n_random_parents() } );
    #print STDERR "n_random_parents:  ", $self->n_random_parents(), "\n";
    #  print STDERR "Storing pedcheck for $accid in CheckPedigrees obj.\n";
     
    $self->pedigree_checks()->{$accid} = $pedcheck;
    $pedcheck->compare_to_random_parents($accid_gts);
     my $progeny_string = $self->pedigrees()->id_node()->{$accid}->offspring();
   my @progeny = split(" ", $progeny_string);
    print $fhsummary $pedcheck->as_string(), "  ", scalar @progeny, "\n";
    print $fh27 $pedcheck->as_string_27(), "\n" if(defined $self->triplecounts27_filename());
    print $fh14 "$accid  ", $accgtobj->quality_counts(), "  $matid $patid  ",
      $pedcheck->as_string_14(), "\n" if(defined $self->counts14_filename());
    $n_pedigree_checks++;
    print STDERR "# pedigree checks done: $n_pedigree_checks \n" if(($n_pedigree_checks % $self->progress_report_interval()) == 0);
  }
  print $fhsummary $self->as_string();
  printf($fhsummary "# Female parent - offspring hgmr clusters: n pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $self->m_hgmr_cluster());
  printf($fhsummary "#   Male parent - offspring hgmr clusters: n_pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $self->p_hgmr_cluster());
  printf($fhsummary "# agmr of parents                clusters: n_pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $self->mp_agmr_cluster());
  printf($fhsummary "# rx01                           clusters: n_pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $self->rx01_cluster());
  printf($fhsummary "# r0x1                           clusters: n_pts: %4i  k-means: %4i %4i %8.5f  kde: %4i %4i %8.5f \n", $self->r0x1_cluster());
}

sub as_string{
  my $self = shift;
  my $gtsetobj = $self->genotypes_set();
  my $the_string = '';
  $the_string .= "# n_accessions: " . $gtsetobj->n_accessions() . "\n";
  $the_string .= "# n_markers: " . scalar @{$gtsetobj->marker_ids()} . "\n";
  $the_string .= "# delta: " . $gtsetobj->delta() . "\n";

  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
    #   my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
    #    $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
    $the_string .= $pedchk->as_string();
    $the_string .= "\n";
  }
  return $the_string;
}

sub m_hgmr_cluster {
  my $self = shift;

  my @hgmrs = ();
  print STDERR "In m_hgmr_cluster. number of pedigree_checks: ", scalar keys %{$self->pedigree_checks()}, "\n";
  while (my($aid, $pedcheck) = each %{$self->pedigree_checks()}) {
    # print STDERR "aid: $aid  m_hgmr: ", $pedcheck->m_hgmr(), "\n";
    push @hgmrs, $pedcheck->m_hgmr();
  }
  print STDERR "size of hgmrs array: ", scalar @hgmrs, "\n";
  my $cluster1d_obj = Cluster1d->new({label => 'hgmr(female parent)', xs => \@hgmrs});
  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
 $self->Fhgmr_h($kde_h_opt);
  return ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt); # $cluster1d_obj->one_d_2cluster();
}

  sub p_hgmr_cluster {
    my $self = shift;

    my @hgmrs = ();
    while (my($aid, $pedcheck) = each %{$self->pedigree_checks()}) {
      push @hgmrs, $pedcheck->p_hgmr();
    }
    my $cluster1d_obj = Cluster1d->new({label => 'hgmr(male parent)', xs => \@hgmrs});
    my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
    $self->Mhgmr_h($kde_h_opt);
    return ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt); # $cluster1d_obj->one_d_2cluster();
  }

sub mp_agmr_cluster {
  my $self = shift;

  my @agmrs = ();
  while (my($aid, $pedcheck) = each %{$self->pedigree_checks()}) {
    push @agmrs, $pedcheck->mp_agmr();
  }
  my $cluster1d_obj = Cluster1d->new({label => 'agmr(between parents)', xs => \@agmrs});
  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
  $self->FMagmr_h($kde_h_opt);
  return ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt); # $cluster1d_obj->one_d_2cluster();
}

sub rx01_cluster {
  my $self = shift;

  my @rs = ();
  while (my($aid, $pedcheck) = each %{$self->pedigree_checks()}) {
    push @rs, $pedcheck->rx01();
  }
  my $cluster1d_obj = Cluster1d->new({label => 'r0x1', xs => \@rs});
  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
  $self->rx01_h($kde_h_opt);
  return ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt); # $cluster1d_obj->one_d_2cluster();
}

sub r0x1_cluster {
  my $self = shift;

  my @rs = ();
  while (my($aid, $pedcheck) = each %{$self->pedigree_checks()}) {
    push @rs, $pedcheck->r0x1();
  }
  my $cluster1d_obj = Cluster1d->new({label => 'rx01', xs => \@rs});
  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
  $self->r0x1_h($kde_h_opt);
  return ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $kde_n_L, $kde_n_R, $kde_h_opt); # $cluster1d_obj->one_d_2cluster();
}

sub categorize{
  my $self = shift;
  my %accid_catstr = ();
  my %catstr_ids = ();
  while (my($accid, $pedcheck) = each %{$self->pedigree_checks()}){
    my $category_string = _x_to_012($pedcheck->mp_agmr(), $self->FMagmr_h());

    $category_string .= _x_to_012($pedcheck->m_hgmr(), $self->Fhgmr_h());
    $category_string .= _x_to_012($pedcheck->p_hgmr(), $self->Mhgmr_h());

    $category_string .= _x_to_012($pedcheck->rx01(), $self->rx01_h());
    $category_string .= _x_to_012($pedcheck->r0x1(), $self->r0x1_h());
    $accid_catstr{$accid} = $category_string;
    if(! exists $catstr_ids{$category_string}){
      $catstr_ids{$category_string} = '';
    }
    $catstr_ids{$category_string} .= "$accid ";
  }
  $self->category_ids(\%catstr_ids);
  $self->id_category(\%accid_catstr);
  # while(my($id, $cs) = each %accid_catstr){
  #   print "$id  $cs \n";
  # }
  my @sortedcs = sort keys %catstr_ids;
  my $n_in_all_categories = 0;
  #while(my($cs, $ids) = each %catstr_ids){
    for my $cs (@sortedcs){
      my $ids = $catstr_ids{$cs};
      my $n_in_category = scalar split(" ", $ids);
      $n_in_all_categories += $n_in_category;
    print join('   ', split('', $cs)), "    $n_in_category\n"; # "  $ids \n";
    }
  print "all categories:   $n_in_all_categories \n";
}


### non methods

sub _get_a_pedcheck{
  my $aid_gts = shift;
  my ($matid, $patid, $accid) = @_;
  my $accgtobj = $aid_gts->{$accid} // return undef;
  my $matgtobj = $aid_gts->{$matid} // return undef;
  my $patgtobj = $aid_gts->{$patid} // return undef;
  return  PedigreeCheck->new( { mat_gtsobj => $matgtobj, pat_gtsobj => $patgtobj, acc_gtsobj => $accgtobj } );
}

sub _x_to_012{
  my $x = shift;
  my $h = shift;
  my $result;
  if($x > $h){ # above threshold
    $result = 1;
  }elsif($x >= 0){ # below threshold
    $result = 0;
  }else{ # negative; invalid
    $result = 2;
  }
  return $result;
}


# unused

# sub print_summary{
#   my $self = shift;
#   my $output_filename = shift;
#   open my $fh, ">", "$output_filename";
#   my $the_string = "#n_accessions: " . $gtsetobj->n_accessions() . "\n";
#   $the_string .= "#n_markers: " . $gtsetobj->n_markers() . "\n";
#   $the_string .= "#delta: " . $gtsetobj->delta() . "\n";
# print $fh "$the_string";
#   while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
#     #   my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
#     #    $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
#     print $fh  $pedchk->as_string(), "\n";
#   }
#   close $fh;
#   #return $the_string;
# }


# sub as_string_ns{
#   my $self = shift;
#   my $gtsetobj = $self->genotypes_set();
#   my $the_string = '';
#   $the_string .= "#n_accessions: " . $gtsetobj->n_accessions() . "\n";
#   $the_string .= "#n_markers: " . $gtsetobj->n_markers() . "\n";
#   $the_string .= "#delta: " . $gtsetobj->delta() . "\n";

#   while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
#     my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
#     $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
#     $the_string .= $pedchk->as_string_ns() . "\n";
#   }
#   return $the_string;
# }

# sub as_string_Ns{
#   my $self = shift;
#   my $gtsetobj = $self->genotypes_set();
#   my $the_string = '';
#   $the_string .= "# n_accessions: " . $gtsetobj->n_accessions() . "\n";
#   $the_string .= "# n_markers: " . $gtsetobj->n_markers() . "\n";
#   $the_string .= "# delta: " . $gtsetobj->delta() . "\n";

#   #  while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
#   for my $accid (@{$gtsetobj->accession_ids()}) {
#     my $gtsobj = $gtsetobj->accid_genotypes()->{$accid};
#     my $pedchk = $self->pedigree_checks()->{$accid} // undef;
#     $the_string .= $accid . "  " . $gtsobj->quality_counts() . "  ";
#     $the_string .= (defined $pedchk)? $pedchk->as_string_Ns() . "\n" : "  No pedigree \n";
#   }
#   return $the_string;
# }

# sub as_string_xs{
#   my $self = shift;

#   my $the_string = '';
#   while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
#     $the_string .= $pedchk->as_string_xs() . "   ";
#     my $mrchk = $self->matrand_checks()->{$accid};
#     $the_string .= $mrchk->as_string_xs() . "   ";
#     my $prchk = $self->patrand_checks()->{$accid};
#     $the_string .= $prchk->as_string_xs() . "   ";
#     $the_string .= $self->randrand_checks()->{$accid}->as_string_xs() . "   ";
#     $the_string .= "\n";
#   }
#   return $the_string;
# }

# sub as_string_zs{
#   my $self = shift;

#   my $the_string = '';
#   while (my ($accid, $pedchk) = each %{$self->pedigree_checks()}) {
#     $the_string .= $pedchk->as_string_zs() . "   ";
#     my $mrchk = $self->matrand_checks()->{$accid};
#     $the_string .= $mrchk->as_string_zs() . "   ";
#     my $prchk = $self->patrand_checks()->{$accid};
#     $the_string .= $prchk->as_string_zs() . "   ";
#     $the_string .= $self->randrand_checks()->{$accid}->as_string_zs() . "   ";
#     $the_string .= "\n";
#   }
#   return $the_string;
# }

# sub _nxyzs{	     # given mat, pat, child gts, get n000, n001, etc.
#   my $mat_gts = shift;
#   my $pat_gts = shift;
#   my $child_gts = shift;

#   my %n_c = ('000' => 0, '001' => 0, '002' => 0,   '010' => 0, '011' => 0, '012' => 0,
# 	     '020' => 0, '021' => 0, '022' => 0,   '110' => 0, '111' => 0, '112' => 0,
# 	     '120' => 0, '121' => 0, '122' => 0,   '220' => 0, '221' => 0, '222' => 0,
# 	     'X' => 0);
 
#   my @m_gts = split("", $mat_gts);
#   my @p_gts = split("", $pat_gts);
#   my @ch_gts = split("", $child_gts);
#   die "gt string length prob. \n" if (scalar @m_gts != scalar @p_gts  or  scalar @m_gts != scalar @ch_gts);
#   #  die "gt string length prob. " . length $mat_gts . " " . length $pat_gts . " " . length $child_gts . "\n"
#   #    if (length $mat_gts != length $pat_gts  or  length $mat_gts != length $child_gts);

#   #  for (my $i = 0; $i < length $mat_gts; $i++) {
#   while (my($i, $c) = each @ch_gts) {
#     my $m = $m_gts[$i];
#     my $p = $p_gts[$i];
#     my $tr; # e.g. '001'
#     #    print STDERR "$c  $m $p  $tr \n";
#     #   $tr .= $c;
#     if (uc $c eq 'X'  or  uc $m eq 'X'  or  uc $p eq 'X') { # missing data case
#       $tr = 'X';
#     } else {			# all 3 gts present. 
#       $tr = ($m < $p)? "$m$p$c" : "$p$m$c";
#     }
#     $n_c{$tr}++;
#   }
#   my @nxyzs = ($n_c{'000'}, $n_c{'001'}, $n_c{'002'},
# 	       $n_c{'010'}, $n_c{'011'}, $n_c{'012'},
# 	       $n_c{'020'}, $n_c{'021'}, $n_c{'022'},
# 	       $n_c{'110'}, $n_c{'111'}, $n_c{'112'},
# 	       $n_c{'120'}, $n_c{'121'}, $n_c{'122'},
# 	       $n_c{'220'}, $n_c{'221'}, $n_c{'222'}, $n_c{'X'}); #, scalar @ch_gts);
#   die if(scalar @ch_gts != sum(@nxyzs));
#   return \@nxyzs;
# }

1;
