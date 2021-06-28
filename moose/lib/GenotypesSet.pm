package GenotypesSet;
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
# use Marker;

# accession id and string of genotypes

has gt_matrix_filename => (
			   isa => 'Str',
			   is => 'ro',
			   required => 1,
			  );

has delta => (
	      isa => 'Num',
	      is => 'ro',
	      default => 0.5,
	     );

has n_accessions => (
		     isa => 'Num',
		     is => 'rw',
		     default => -1,
		    );

has accession_ids => (isa => 'ArrayRef[Str]',
		      is => 'rw',
		      required => 0,
		     );


has marker_ids => ( isa => 'ArrayRef[Str]', 
		    is => 'rw',
		    required => 0,
		    # default =>  sub { [] },
		  );

has accid_genotypes => (		  # Can point to either allgenotypes or to goodgenotypes
			isa => 'HashRef', # values are Genotypes objects
			is => 'rw',
			default => sub { {} },
		       );

has allmarker_ids => ( isa => 'ArrayRef[Str]',
		    is => 'rw',
		    required => 0,
		    # default =>  sub { [] },
		  );
has accid_allgenotypes => (		  # 
			isa => 'HashRef', # values are Genotypes objects
			is => 'rw',
			default => sub { {} },
		       );

# 'good' genotype info

has marker_gt_counts  => ( isa => 'ArrayRef[ArrayRef[Int]]', 
			   is => 'rw',
			   required => 0,
			   # default =>  sub { [] },
			 );
has max_bad_gt_fraction => (
		    isa => 'Num',
		    is => 'rw',
		    default => 1,
			   );
has min_hw_qual_param => (
			  isa => 'Num',
			  is => 'rw',
			  default => 0,
			  );
has goodmarker_ids => ( isa => 'ArrayRef[Str]', 
		    is => 'rw',
		    required => 0,
		    # default =>  sub { [] },
		  );

has accid_goodgenotypes => (		  # 
			isa => 'HashRef', # values are Genotypes objects
			is => 'rw',
			default => sub { {} },
		       );

has gt_set_as_01234string => ( # 01234
			 isa => 'Str',
			 is => 'rw',
			 default => '',
			     );

has progress_report_interval => (
				 isa => 'Int',
				 is => 'ro',
				 default => 500,
				 );

sub BUILD {
  my $self = shift;
  my $filename = $self->gt_matrix_filename();

  open my $fh, "<", "$filename" or die "Couldn't open $filename for reading.\n";
  my @lines = <$fh>;

  # my $line1;
  # while(($line1 = <$fh>) =~ /^\s*#/){
  # };

  while($lines[0] =~ /^\s*#/){
    shift @lines; # just skipping comments
  }
  my @acc_ids = ();
  my $marker_id_line = shift @lines;
  my @mrkr_ids = split(" ", $marker_id_line); # yes, input is whitespace separated.
  my $M = shift @mrkr_ids;
  die if($M ne 'MARKER');
  $self->allmarker_ids(\@mrkr_ids);
  $self->marker_ids($self->allmarker_ids());
  $self->goodmarker_ids($self->allmarker_ids());
  my @mrkr_gtcounts = ();
  for (@mrkr_ids) {
    push @mrkr_gtcounts, [0, 0, 0, 0, 0];
  }
  $self->marker_gt_counts(\@mrkr_gtcounts);

  my $n_markers = scalar @mrkr_ids; print STDERR "# n markers: $n_markers \n";
  my %id_gts = ();
  my $accessions_read = 0;
 #  while (my $s = <$fh>) {
    for my $s (@lines){
    next if($s =~ /^\s*($|#)/);	# skip whitespace-only lines and comments
    my @gts = split(" ", $s);
    my $acc_id = shift @gts;
    push @acc_ids, $acc_id;
    #   print STDERR "acc id: [$acc_id] \n";
    next if (scalar @gts == 0);
    my ($s01234, $q_counts_str) = $self->resolve_to_01234(\@gts);
    my ($n0, $n1, $n2, $n3, $n4) = split(" ", $q_counts_str);
    my $acc_bad_fraction = ($n3+$n4)/($n0+$n1+$n2+$n3+$n4);
   # print STDERR "$acc_id  $q_counts_str  $acc_bad_fraction \n";

    die "length: ", length $s01234, "  n markers: $n_markers" if(length $s01234 != $n_markers);
    $id_gts{$acc_id} = Genotypes->new({id => $acc_id, genotypes => $s01234, quality_counts => $q_counts_str});
    $accessions_read++;
    print STDERR "# accessions read: $accessions_read \n" if(($accessions_read % $self->progress_report_interval()) == 0);
  }
  $self->n_accessions(scalar keys %id_gts);
  $self->accession_ids(\@acc_ids);
#  $self->n_markers($n_markers);
  $self->accid_allgenotypes(\%id_gts);
  $self->accid_genotypes($self->accid_allgenotypes());

  if( ($self->max_bad_gt_fraction() < 1.0) or ($self->min_hw_qual_param() > 0.0) ){
    print STDERR "# removing bad markers. max bad fraction: ", $self->max_bad_gt_fraction(), "   min hw qual: ", $self->min_hw_qual_param(), "\n";
    $self->clean_marker_set(); #$max_bad_gt_fraction, $min_hw_qual_param);
  }
  #$self->gt_set_as_string($out_string_012Xx);
}

sub resolve_to_01234{ # take gts which are non-integers in range [0,2], and round to 0, 1, 2, (or X if not withing delta of 0, 1, or 2)
  my $self = shift;
  my $gts = shift;		# array ref 
  my $delta = $self->delta();

  my $gtstr01234;
  #print "genotypes:    ", $gts->[0], "\n";
  if (scalar @$gts == 1) {	# already resolved to 01234
    $gtstr01234 = $gts->[0];
  } else {			# float values in range [0,2]
    my @gts01234 = ();
    for my $agt (@$gts) {
      if ($agt <= $delta) {	# round to 0
	push @gts01234, 0;
      } elsif (abs($agt - 1) <= $delta) { 
	push @gts01234, 1;
      } elsif ($agt >= 2-$delta) { # round to 2
	push @gts01234, 2;
      } else {
	push @gts01234, 3; 
	# if ($agt < 1) {
	#   push @gts01234, 3;
	# } else {
	#   push @gts01234, 4;
	# }
      }
    }

    $self->increment_marker_gt_counts(\@gts01234);

    $gtstr01234 = join('', @gts01234);
  }
  #   $self->increment_marker_gt_counts_alt($gtstr01234);
  my $n0 = ($gtstr01234 =~ tr/0//);
  my $n1 = ($gtstr01234 =~ tr/1//);
  my $n2 = ($gtstr01234 =~ tr/2//);
  my $nX = ($gtstr01234 =~ tr/3//);
  my $nx = 0 ; # ($gtstr01234 =~ tr/4//);
  return ($gtstr01234, "$n0 $n1 $n2 $nX $nx"); # e.g. 001012010201 etc.
}

sub increment_marker_gt_counts_alt{
  my $self = shift;
  my $igt_string = shift;
  my $mrkr_gt_counts = $self->marker_gt_counts();
  for my $i (0 .. (length $igt_string) -1) {
    my $igt = substr($igt_string, $i, 1);
    $mrkr_gt_counts->[$i]->[$igt]++;
  }
}

sub increment_marker_gt_counts{
  my $self = shift;
  my $igts = shift;		# array ref w 01234
  my $mrkr_gt_counts = $self->marker_gt_counts();
  while (my($i, $igt) = each @$igts) {
    $mrkr_gt_counts->[$i]->[$igt]++;
  }
}

sub summary_info_string{
  my $self = shift;
  my $str =  "# n_accessions: " . $self->n_accessions() . "\n";
  $str .= "# delta: " . $self->delta() . "\n";
  $str .=  "# n_markers. all: " .  scalar @{$self->allmarker_ids()} .
    "  max bad fraction: " . $self->max_bad_gt_fraction() .
    "  n markers used:  " .  scalar @{$self->goodmarker_ids()} .
    "\n";
  return $str;
}

sub as_string{
  my $self = shift;
  my $X = shift // 0; # default is 01234; call with $X = 1 to get 012Xx
  my @marker_ids =  @{$self->marker_ids()};
#  my $s = "# delta: " . $self->delta() . "\n";
#  $s .= "# n_markers. all: " .  scalar @{$self->allmarker_ids()} . "  good (max bad fraction: " . $self->max_bad_gt_fraction() . "):  " .  scalar @{$self->marker_ids()} . "\n";
  my $s = $self->summary_info_string();
  $s .= "#FULL_MARKER_SET " . join(" ", @marker_ids) . "\n";
  my @count_labels = ('n0', 'n1', 'n2', 'nX', 'nx');
  for my $j (0..4) {
    $s .= "# " . $count_labels[$j] . " ";
    #for my $mid (@marker_ids) {
    while (my($i, $mid) = each @marker_ids) {
      #    $s .= sprintf(" %1i", $self->markerid_marker()->{$mid}->counts()->[$j] // 0);
      $s .= sprintf(" %1i", $self->marker_gt_counts()->[$i]->[$j] // 0);
    }
    $s .= "\n";
  }

  #while (my ($accid, $gtsobj) = each %{$self->accid_genotypes()}) {
  # my $agtsobj = $self->accid_genotypes()->{$self-accession_ids()->[0]};
  $s .= "MARKER " . join(" ", @{$self->goodmarker_ids()}) . "\n";
  for my $acc_id (@{$self->accession_ids()}) {
    my $gtsobj = $self->accid_genotypes()->{$acc_id};
    my $gtstr = $gtsobj->as_string();
  #   if($X){
  #   $gtstr =~ s/3/X/g;
  #   $gtstr =~ s/4/x/g;
  # }#
    $s .= $gtstr . "\n";
  }
  return $s;
}

sub marker_qual_string{
  my $self = shift;
  my @marker_ids =  @{$self->marker_ids()};
  my $s = "# delta: " . $self->delta() . "\n";
  $s .= "# marker_id  n0 n1 n2  nX nx  hw\n";
  while (my($i, $mid) = each @marker_ids) {
    $s .= "$mid  ";
    my @counts =  @{$self->marker_gt_counts()->[$i]};
    my $n012 = sum(@counts[0..2]);
    die if($self->n_accessions() != sum(@counts));
    for my $j (0..4) {
      my $n = $counts[$j] // 0;
      my $denom = ($j <= 2)? $n012 : sum(@counts);
      $s .= sprintf(" %1i %7.5f", $n, ($denom>0)? $n/$denom : -1);
    }
    $s .= "  " . hardy_weinberg(\@counts) . "\n";
  }
  return $s;
}

sub clean_marker_set{
  my $self = shift;
  my $max_bad_gt_fraction = shift // $self->max_bad_gt_fraction();
  my $min_hw_qual_param = shift // $self->min_hw_qual_param();

  my @good_marker_ids = ();
  my @mrkr_ids = @{$self->marker_ids()};
  my @marker_thumbsupdown = (1) x scalar @mrkr_ids; # 1: good, 0: bad
  my $mrkr_gt_counts = $self->marker_gt_counts();
  while (my ($i, $markerid) = each @mrkr_ids) {
    my $marker_counts = $self->marker_gt_counts()->[$i];
    my $n_markers = scalar @{$self->marker_ids()};
    my $n_accessions = $self->n_accessions();
    my $hw_qual_param = hardy_weinberg($marker_counts);
    my $bad_gt_fraction = ($marker_counts->[3] + $marker_counts->[4])/$n_accessions // 1;
    if ( ($bad_gt_fraction > $max_bad_gt_fraction) or ($hw_qual_param < $min_hw_qual_param) ) { # bad
      $marker_thumbsupdown[$i] = 0; # mark this marker as bad
    } else {			    # OK
      push @good_marker_ids, $markerid;
    }
  }
  print "In clean_marker_set size of marker_thumbsupdown array: ", scalar @marker_thumbsupdown, "\n";
  #$self->n_markers(scalar @good_marker_ids);
  $self->goodmarker_ids(\@good_marker_ids);

  my @accids =  keys %{$self->accid_genotypes()};
  $self->max_bad_gt_fraction($max_bad_gt_fraction);
  for my $accid (@accids){
    my $gtsobj = $self->accid_genotypes()->{$accid};
    # print "length gtsobj->genotypes(): ", length $gtsobj->genotypes(), "\n";
    my $cleaned_gtsobj = $gtsobj->construct_cleaned_genotypes_object(\@marker_thumbsupdown);

     # print "length cleaned_gtsobj->genotypes(): ", length $cleaned_gtsobj->genotypes(), "\n";
  #  $gtsobj->remove_bad_markers(\@marker_thumbsupdown);
    $self->accid_goodgenotypes()->{$accid} = $cleaned_gtsobj;
  }
  $self->accid_genotypes($self->accid_goodgenotypes());
  ## $self->marker_ids($self->goodmarker_ids());
  
}

######### non-methods ##################

sub hardy_weinberg{ # if hardy-weinberg eq. frequencies, return value close to 1.
  my $counts = shift;		#
  my $n0 = $counts->[0] // 0;
  my $n1 = $counts->[1] // 0;
  my $n2 = $counts->[2] // 0;
  my $ntotal = $n0 + $n1 + $n2;
  my $hw_one = ($ntotal > 0)? ($n0**0.5 + $n2**0.5)/$ntotal**0.5 : -1;
  return $hw_one;
}

1;
