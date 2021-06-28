package PedigreeCheck;		# checking a single pedigree.
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

has acc_gtsobj => (
		   isa => 'Object',
		   is => 'ro',
		   required => 1,
		  );
has mat_gtsobj => (
		   isa => 'Object',
		   is => 'ro',
		   required => 1,
		  );
has pat_gtsobj => (
		   isa => 'Object',
		   is => 'ro',
		   required => 1,
		  );

has am_distances => (
		     isa => 'ArrayRef', # [agmr, hgmr]
		     is => 'rw',
		     default => sub { [-1, -1] }, # agmr hgmr
		    );

has ap_distances => (
		     isa => 'ArrayRef',
		     is => 'rw',
		     default => sub { [-1, -1] },
		    );

has mp_distances => (
		     isa => 'ArrayRef',
		     is => 'rw',
		     default => sub { [-1, -1] },
		    );

has nxyzs => (
	      isa => 'ArrayRef[Int]',
	      is => 'rw',
	      required => 0,
	     );

has mxyzs => (		# m020 = N020 + N202, m012 = N012 + N210, etc.
	      isa => 'ArrayRef[Int]',
	      is => 'rw',
	      required => 0,
	     );

has n00x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n01x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n02x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n11x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n12x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has n22x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has nX => (
	   isa => 'Int',
	   is => 'rw',
	   default => 0,
	  );

has Nxyzs => (
	      isa => 'ArrayRef[Int]',
	      is => 'rw',
	      required => 0,
	     );

has m_hgmr_n_d => ( 
		   isa => 'ArrayRef[Int]', # [numerator, denominator]
		   is => 'rw',
		   required => 0,
		  );

has p_hgmr_n_d => ( 
		   isa => 'ArrayRef[Int]', # [numerator, denominator]
		   is => 'rw',
		   required => 0,
		  );

has mp_agmr_n_d => (  # numerator and denom of agmr between mat & pat.
		    isa => 'ArrayRef[Int]', # [numerator, denominator]
		    is => 'rw',
		    required => 0,
		   );

has rx01_n_d => (	      # 
	      isa => 'ArrayRef[Int]', # [numerator, denominator]
	      is => 'rw',
	      required => 0,
	     );

has r0x1_n_d => (	      #
	      isa => 'ArrayRef[Int]', # [numerator, denominator]
	      is => 'rw',
	      required => 0,
	     );

has n_random_parents => (
		 isa => 'Int',
		 is => 'rw',
		 default => 0,
		);

has mhgmr_rank => ( # replace maternal parent with random accession from data set.
		    isa => 'Int', # [n_less, n_total]
		    is => 'rw',
		    default => 0,
		  );

has phgmr_rank => ( # replace maternal parent with random accession from data set.
		    isa => 'Int', # [n_less, n_total]
		    is => 'rw',
		    default => 0,
		   );

has N00x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N01x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N02x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N10x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N11x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N12x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N20x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N21x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has N22x => (
	     isa => 'ArrayRef[Int]',
	     is => 'rw',
	     default => sub{ [0, 0, 0, 0] },
	    );

has NX => (
	   isa => 'Int',
	   is => 'rw',
	   default => 0,
	  );

sub BUILD{
  my $self = shift;

  $self->triple_counts_27();
  my @Nxyzs = @{$self->Nxyzs()};

  my @N00xs = (@Nxyzs[0..2], sum(@Nxyzs[0..2]));
  $self->N00x(\@N00xs);

  my @N01xs = (@Nxyzs[3..5], sum(@Nxyzs[3..5]));
  $self->N01x(\@N01xs);

  my @N02xs = (@Nxyzs[6..8], sum(@Nxyzs[6..8]));
  $self->N02x(\@N02xs);


  my @N10xs = (@Nxyzs[9..11], sum(@Nxyzs[9..11]));
  $self->N10x(\@N10xs);

  my @N11xs = (@Nxyzs[12..14], sum(@Nxyzs[12..14]));
  $self->N11x(\@N11xs);

  my @N12xs = (@Nxyzs[15..17], sum(@Nxyzs[15..17]));
  $self->N12x(\@N12xs);


  my @N20xs = (@Nxyzs[18..20], sum(@Nxyzs[18..20]));
  $self->N20x(\@N20xs);

  my @N21xs = (@Nxyzs[21..23], sum(@Nxyzs[21..23]));
  $self->N21x(\@N21xs);

  my @N22xs = (@Nxyzs[24..26], sum(@Nxyzs[24..26]));
  $self->N22x(\@N22xs);


  $self->NX($Nxyzs[27]);

  $self->nxyzs(	    # not distinguishing maternal and paternal parents
	       [ $Nxyzs[0],  $Nxyzs[1],  $Nxyzs[2], # n000 = N000
		 $Nxyzs[3] + $Nxyzs[9],  $Nxyzs[4] + $Nxyzs[10],  $Nxyzs[5] + $Nxyzs[11], # N010 + N100, 
		 $Nxyzs[6] + $Nxyzs[18],  $Nxyzs[7] + $Nxyzs[19],  $Nxyzs[8] + $Nxyzs[20],
		 $Nxyzs[12],  $Nxyzs[13],  $Nxyzs[14],
		 $Nxyzs[15] + $Nxyzs[21],  $Nxyzs[16] + $Nxyzs[22],  $Nxyzs[17] + $Nxyzs[23],
		 $Nxyzs[24], $Nxyzs[25], $Nxyzs[26],
		 $Nxyzs[27] ]
	      );

  my @nxyzs = @{$self->nxyzs()}; # e.g. @nxyzs[0..2] = (n000, n001, n002)

  my @n00xs = (@nxyzs[0..2], sum(@nxyzs[0..2]));
  $self->n00x(\@n00xs);

  my @n01xs = (@nxyzs[3..5], sum(@nxyzs[3..5]));
  $self->n01x(\@n01xs);

  my @n02xs = (@nxyzs[6..8], sum(@nxyzs[6..8]));
  $self->n02x(\@n02xs);

  my @n11xs = (@nxyzs[9..11], sum(@nxyzs[9..11]));
  $self->n11x(\@n11xs);

  my @n12xs = (@nxyzs[12..14], sum(@nxyzs[12..14]));
  $self->n12x(\@n12xs);

  my @n22xs = (@nxyzs[15..17], sum(@nxyzs[15..17]));
  $self->n22x(\@n22xs);

  $self->nX($Nxyzs[27]);


  $self->mxyzs(
	       [ $Nxyzs[0] + $Nxyzs[26], # m000 = N000 + N222
		 $Nxyzs[1] + $Nxyzs[25], # m001 = N001 + N221
		 $Nxyzs[2] + $Nxyzs[24],
		 $Nxyzs[3] + $Nxyzs[23],
		 $Nxyzs[4] + $Nxyzs[22], # m010 = N010 + N212
		 $Nxyzs[5] + $Nxyzs[21],
		 $Nxyzs[6] + $Nxyzs[20],
		 $Nxyzs[7] + $Nxyzs[19], # N021 + N201
		 $Nxyzs[8] + $Nxyzs[18],
		 $Nxyzs[9] + $Nxyzs[17],
		 $Nxyzs[10] + $Nxyzs[16],
		 $Nxyzs[11] + $Nxyzs[15],
		 $Nxyzs[12] + $Nxyzs[14],
		 $Nxyzs[13],
		 $Nxyzs[27] ]
	      );

  
  my $rx01_num = sum(@Nxyzs[1,10,19,7,16,25]); # n001 + n101 + n201 + n021 + n121 + n221 i.e. n0x1 + n2x1
  my $rx01_denom = sum(@Nxyzs[0,9,20,6,17,26]) + $rx01_num; # n000 + n100 + n020 + n022 + n122 + n222  +  $r_num  i.e. n0x0 + n2x2 + $r_num
  $self->rx01_n_d([$rx01_num, $rx01_denom]);

   my $r0x1_num = sum(@Nxyzs[1,4,7,19,22,25]); # n001 + n011 + n021 + n201 + n211 + n221 i.e. n0x1 + n2x1
  my $r0x1_denom = sum(@Nxyzs[0,3,6,20,23,26]) + $r0x1_num; # n000 + n010 + n020 + n202 + n212 + n222  +  $r_num  i.e. n0x0 + n2x2 + $r_num
  $self->r0x1_n_d([$r0x1_num, $r0x1_denom]);


  
  #  $self->distances();
  
}

sub as_string_ns{
  my $self = shift;
  my $the_string = sprintf("%s %s  ", $self->mat_gtsobj->id(), $self->pat_gtsobj->id());
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n00x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n01x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n02x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n11x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n12x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->n22x()});
  $the_string .= sprintf("%2i", $self->nX());

  return $the_string;
}

sub as_string_27{
  my $self = shift;
  #my $the_string = sprintf("%s %s  ", $self->mat_gtsobj->id(), $self->pat_gtsobj->id());
   my $the_string = sprintf("%s %s %s %s  ", $self->acc_gtsobj()->id(),  $self->acc_gtsobj()->quality_counts(), $self->mat_gtsobj()->id(), $self->pat_gtsobj->id());
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N00x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N01x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N02x()});

  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N10x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N11x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N12x()});

  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N20x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N21x()});
  $the_string .= sprintf("%2i %2i %2i %2i   ", @{$self->N22x()});

  $the_string .= sprintf("%2i", $self->NX());

  # $the_string .= sprintf("  %7.5f %7.5f", @{$self->am_distances()});
  # $the_string .= sprintf("  %7.5f %7.5f", @{$self->ap_distances()});
  # $the_string .= sprintf("  %7.5f %7.5f", @{$self->mp_distances()});

  $the_string .= sprintf("   %5i %5i  %7.5f   %5i %5i  %7.5f", @{$self->m_hgmr_n_d()}, $self->m_hgmr(), @{$self->p_hgmr_n_d()}, $self->p_hgmr());
  $the_string .= sprintf("   %5i %5i  %7.5f", @{$self->mp_agmr_n_d()}, $self->mp_agmr());
  
  return $the_string;
}

sub as_string_14{
  my $self = shift;
  #my $the_string = sprintf("%s %s  ", $self->mat_gtsobj->id(), $self->pat_gtsobj->id());
  my $the_string = '';
  #sprintf("%s %s %s %s  ", $self->acc_gtsobj()->id(),  $self->acc_gtsobj()->quality_counts(), $self->mat_gtsobj()->id(), $self->pat_gtsobj->id());

  my @ms = @{$self->mxyzs()};
  # $the_string .= sprintf("%2i %2i %2i   ", @ms[0,1,2]);
  # $the_string .= sprintf("%2i %2i %2i   ", @ms[3,4,5]);
  # $the_string .= sprintf("%2i %2i %2i   ", @ms[6,7,8]);
  # $the_string .= sprintf("%2i %2i %2i   ", @ms[9,10,11]);
  # $the_string .= sprintf("%2i %2i %2i   ", @ms[12,13,14]);

    $the_string .= sprintf("%1i %1i %1i  ", @ms[0,1,2]);
  $the_string .= sprintf("%1i %1i %1i  ", @ms[3,4,5]);
  $the_string .= sprintf("%1i %1i %1i  ", @ms[6,7,8]);
  $the_string .= sprintf("%1i %1i %1i  ", @ms[9,10,11]);
  $the_string .= sprintf("%1i %1i  %1i  ", @ms[12,13,14]);

#  $the_string .= sprintf("%7.5f   %7.5f   %7.5f", $self->m_hgmr(), $self->p_hgmr(), $self->mp_agmr());
  return $the_string;
}


sub as_string{
  my $self = shift;
  my $the_string = sprintf("%s  %s  %s  %s   ", $self->acc_gtsobj()->id(),  $self->acc_gtsobj()->quality_counts(), $self->mat_gtsobj()->id(), $self->pat_gtsobj->id());
  $the_string .= sprintf("%7.5f %7.5f %7.5f %7.5f %7.5f   ", $self->m_hgmr(), $self->p_hgmr(), $self->mp_agmr(), $self->rx01(), $self->r0x1());
  $the_string .= $self->as_string_14() . " ";
  $the_string .= sprintf("%1i %1i %1i", $self->mhgmr_rank(), $self->phgmr_rank(), $self->n_random_parents());
  
  return $the_string;
}

sub as_string_xs{
  my $self = shift;
  my @n_xyzs = @{$self->nxyzs()};
  my $the_string = sprintf("%s %s %s  ", $self->accession_id(), $self->maternal_id(), $self->paternal_id());
  for my $i (0..5) {
    my $j = 3*$i;
    my $ntot = sum(@n_xyzs[$j..$j+2]);

    $the_string .= ($ntot > 0)? sprintf("%5.4f %5.4f %5.4f %2i  ", map($_/$ntot, @n_xyzs[$j..$j+2]), $ntot ) : sprintf("-1 -1 -1 0  ");
  }
  return $the_string;
}

sub as_string_z{		# lump together n001 n221, etc.
  my $self = shift;
  my $the_string = sprintf("%s %s %s  ", $self->acc_gtsobj->id(), $self->mat_gtsobj->id(), $self->pat_gtsobj->id());

  my $n00n22 = $self->n00x()->[3] + $self->n22x()->[3];
  my $z001221 = ($n00n22 > 0)? ($self->n00x->[1] + $self->n22x->[1])/$n00n22 : -1;
  my $z002220 = ($n00n22 > 0)? ($self->n00x->[2] + $self->n22x->[0])/$n00n22 : -1;
  
  my $n01n12 = $self->n01x()->[3] + $self->n12x()->[3];
  my $z012120 = ($n01n12 > 0)? ($self->n01x->[2] + $self->n12x->[0])/$n01n12 : -1;

  my $n02 = $self->n02x()->[3];
  my $z020022 = ($n02 > 0)? ($self->n02x->[0] + $self->n02x->[2])/$n02 : -1;
  $the_string .= sprintf("%5.4f %3i %5.4f %3i %5.4f %3i %5.4f %3i",
			 $z001221, $n00n22, $z002220, $n00n22, $z012120, $n01n12, $z020022, $n02);
  return $the_string;
}

# sub double_counts_6{		# just
#   my $self = shift;

#   my $mat_gtstr = $self->mat_gtsobj()->genotypes();
#   my $pat_gtstr = $self->pat_gtsobj()->genotypes();
#   my $acc_gtstr = $self->acc_gtsobj()->genotypes();

#   for my $i (0 .. (length $acc_gtstr) - 1){
#     my $acc_gt = substr($acc_gtstr, $i, 1);
#     next if($acc_gt == 3);
    

#   }
# }

sub triple_counts_27{ # get N000, N001, etc. for the pedigree.
  my $self = shift;

  my $mat_gtstr = $self->mat_gtsobj()->genotypes();
  my $pat_gtstr = $self->pat_gtsobj()->genotypes();
  my $acc_gtstr = $self->acc_gtsobj()->genotypes();
  $mat_gtstr =~ s/4/3/g;
  $pat_gtstr =~ s/4/3/g;
  $acc_gtstr =~ s/4/3/g;
  if (length $mat_gtstr != length $acc_gtstr  or  length $pat_gtstr != length $acc_gtstr) {
    print STDERR "genotype string lengths (should all be equal): ",
      length $mat_gtstr, "  ", length $pat_gtstr, "  ", length $acc_gtstr, "\n";
    die "in PedigreeCheck gt string length prob. \n";
  }

  my @tcs = (0) x 64;
  for my $i (0 .. (length $acc_gtstr) - 1) {
    my $m = substr($mat_gtstr, $i, 1);
    my $p = substr($pat_gtstr, $i, 1);
    my $c = substr($acc_gtstr, $i, 1);
    $tcs[16*$m + 4*$p + $c]++;
  }	     # end loop over gts in $acc_gtstr, $mat_gtstr, $pat_gtstr.

  my ($mp_agmr_numerator, $mp_agmr_denominator) = (0, 0);
  my ($m_hgmr_numerator, $m_hgmr_denominator, $n_amok) = (0, 0, 0);
  my ($p_hgmr_numerator, $p_hgmr_denominator, $n_pmok) = (0, 0, 0);
  my @triplecounts_27 = (0) x 28;
  for my $i (0..3) { # maternal gt 
    my $i16 = 16*$i;
    for my $j (0..3) {
      my $j4 = 4*$j;
      for my $k (0..3) {
	my $tcount = $tcs[$i16 + $j4 + $k];
	if ($i <= 2  and $j <= 2  and  $k <= 2) {
	  $triplecounts_27[9*$i + 3*$j + $k] = $tcount;
	}
	if($k <= 2){
	if ($i <= 2) { # mat and acc ok
	  $n_amok++;
	  if ($i != 1  and  $k != 1) {
	    $m_hgmr_denominator += $tcount;
	    $m_hgmr_numerator += $tcount if($i != $k); # 02 or 20
	  }
	}
	if ($j <= 2) { # pat and acc ok
	  $n_pmok++;
	  if ($j != 1  and  $k != 1) {
	    $p_hgmr_denominator += $tcount;
	    $p_hgmr_numerator += $tcount if($j != $k); # 02 or 20
	  }
	}
      }
	if ($i <= 2 and $j <= 2) { # mat and pat ok
	  $mp_agmr_denominator += $tcount;
	  $mp_agmr_numerator += $tcount if($i != $j);
	}
      }
    }
  }

  my $nX = sum(@tcs) - sum(@triplecounts_27);
  $triplecounts_27[27] = $nX;
  $self->Nxyzs( \@triplecounts_27 );

  if ( length $acc_gtstr != sum(@{$self->Nxyzs()}) ) {
    print STDERR $self->acc_gtsobj()->id(), "  ", length $acc_gtstr, "  ";
    while (my ($trip, $c) = each @{$self->Nxyzs()}) {
      printf( STDERR  "%3i %3i  ", $trip, $c);
    }
    print STDERR "\n";
    die;
  }

  $self->m_hgmr_n_d( [$m_hgmr_numerator, $m_hgmr_denominator, $n_amok] );
  $self->p_hgmr_n_d( [$p_hgmr_numerator, $p_hgmr_denominator, $n_pmok] );
  $self->mp_agmr_n_d( [$mp_agmr_numerator, $mp_agmr_denominator] );
}


sub m_hgmr{
  my $self = shift;
  my $nandd = $self->m_hgmr_n_d();
return ($nandd->[1] > 0)? $nandd->[0]/$nandd->[1] : -1;
}

sub p_hgmr{
    my $self = shift;
    my $nandd = $self->p_hgmr_n_d();
    return ($nandd->[1] > 0)? $nandd->[0]/$nandd->[1] : -1;
  }

sub mp_agmr{
  my $self = shift;
  my $nandd = $self->mp_agmr_n_d();
  return ($nandd->[1] > 0)? $nandd->[0]/$nandd->[1] : -1;
}

sub rx01{
  my $self = shift;
  my $nandd = $self->rx01_n_d();
  return  ($nandd->[1] > 0)? $nandd->[0]/$nandd->[1] : -1;
}

sub r0x1{
  my $self = shift;
  my $nandd = $self->r0x1_n_d();
  return  ($nandd->[1] > 0)? $nandd->[0]/$nandd->[1] : -1;
}

sub distances{
  my $self = shift;
  my $mat_gts = $self->mat_gtsobj();
  my $pat_gts = $self->pat_gtsobj();
  my $child_gts = $self->acc_gtsobj();
  $self->am_distances($child_gts->agmr_hgmr($mat_gts));
  $self->ap_distances($child_gts->agmr_hgmr($pat_gts));
  $self->mp_distances($mat_gts->agmr_hgmr($pat_gts));
}

sub compare_to_random_parents{ 
  my $self = shift;		# PedigreeCheck obj.
  my $id_gtsobj = shift;	# keys: ids, values: genotype string
  my $n_random = $self->n_random_parents();
#  print "in compare ... n_random: ", $n_random, "\n";
  return if($n_random == 0);

  my $offspring_gtsobj = $self->acc_gtsobj();

  my $m_hgmr = $self->m_hgmr();
  my $p_hgmr = $self->p_hgmr();
  my $m_lt_count = 0;
  my $p_lt_count = 0;
  my @ids = keys %$id_gtsobj;
  for (1..$n_random) {
    my $random_parent_gtsobj = $id_gtsobj->{$ids[int(rand(scalar @ids))]};
    my $random_parent_hgmr = $offspring_gtsobj->agmr_hgmr($random_parent_gtsobj)->[1];
    $m_lt_count++ if($random_parent_hgmr < $m_hgmr);
    $p_lt_count++ if($random_parent_hgmr < $p_hgmr);
  }
  $self->mhgmr_rank($m_lt_count);
  $self->phgmr_rank($p_lt_count);
}

# sub compare_to_random_parents{ 
#   my $self = shift;		 # PedigreeCheck obj.
#   my $id_genotypestring = shift; # keys: ids, values: genotype string
#   #my $ped_7tfs = shift;	 # seven triple frequency for pedigree parents
#   #  my $offspring_gts = shift;	# 'offspring' genotype string
#   my $parent1_gts = shift; # use this as one 'parent' together with a randomly chosen one, or 'undef' to get both parents randomly chosen.
#   my $n_random = shift;

#   my $offspring_gts = $self->acc_gtsobj()->genotypes();

#   my $i5pct = ($n_random/20)-1;
#   $i5pct = 0 if($i5pct < 0);
#   my @ids = keys %$id_genotypestring;
#   my @stfars = ([], [], [], [], [], [], []); # seven triple frequency array refs
#   for (1..$n_random) {
#     $parent1_gts = $id_genotypestring->{$ids[int(rand(scalar @ids))]} if(!defined $parent1_gts);
#     my $parent2_gts = $id_genotypestring->{$ids[int(rand(scalar @ids))]};
#     my $xyz_counts =  nxxxs($parent1_gts, $parent2_gts, $offspring_gts);
#     my @stfsr = seven_triple_frequencies(@$xyz_counts);
#     my $xX = pop @stfsr;	# 
#     while (my($i, $f) = each @stfsr) {
#       push @{$stfars[$i]}, $f if($f >= 0); #
#     }
#   }
#   my @nlesses = ();
#   while ( my($i, $stfar) = each @stfars) {
#     push @nlesses, (n_less_than($ped_7tfs->[$i], @$stfar), scalar @$stfar);
#   }
#   my $nlessstr = sprintf("%3i %3i  %3i %3i  %3i %3i  %3i %3i  %3i %3i  %3i %3i  %3i %3i", @nlesses);

#   my $fivepctstr = '';
#   for my $stfar (@stfars) {
#     my $nok = scalar @$stfar;
#     $i5pct = int(0.05*$nok - 0.01);
#     my @sorted_stfs = sort {$a <=> $b} @$stfar;
#     $fivepctstr .= '  '. ($nok > 0)? sprintf("%7.5f ", $sorted_stfs[$i5pct]) : "-1 ";
#     #  $fivepctstr .= '  '. sprintf("[%7.5f %7.5f %7.5f] ", $sorted_stfs[0], $sorted_stfs[$n_random-1],  $sorted_stfs[$i5pct]);
#   }

#   return ($nlessstr, $fivepctstr);

# }


# sub seven_triple_freqs{ # x001, x002, x012, x02x, x120, x220, x221. x001 = n001/n00, etc.
#   # i.e. the 7 combinations that should, ideally, not occur for two parents and offspring
#   # in reality we expect these to be small for correct pedigrees, but often non-zero.
#   # e.g. 001 means supposed parents both have genotype 0,
#   # but 'offspring' has genotype 1, which should not occur if pedigree is correct.
#   # x02x = (n020 + n022)/
  
#   my ($n000, $n001, $n002,
#       $n010, $n011, $n012,
#       $n020, $n021, $n022,
#       $n110, $n111, $n112,
#       $n120, $n121, $n122,
#       $n220, $n221, $n222, $nX) = @_;

#   my $n00 = ($n000 + $n001 + $n002);
#   my $x001 = ($n00 > 0)? $n001/$n00 : -1;
#   my $x002 = ($n00 > 0)? $n002/$n00 : -1;

#   my $n01 = ($n010 + $n011 + $n012);
#   my $x012 = ($n01 > 0)? $n012/$n01 : -1;

#   my $n02 = ($n020 + $n021 + $n022);
#   my $x02x = ($n02 > 0)? ($n020 + $n022)/$n02 : -1;

#   my $n11 = ($n110 + $n111 + $n112);
  
#   my $n12 = ($n120 + $n121 + $n122);
#   my $x120 = ($n12 > 0)? $n120/$n12 : -1;

#   my $n22 = ($n220 + $n221 + $n222);
#   my $x220 = ($n22 > 0)? $n220/$n22 : -1;
#   my $x221 = ($n22 > 0)? $n221/$n22 : -1;

#   my $xX = $nX/($n00 + $n01 + $n02 + $n11 + $n12 + $n22 + $nX);
#   #  my $str = sprintf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f", $x001, $x002, $x012, $x02x, $x120, $x220, $x221);
#   #  return $str;
#   return ($x001, $x002, $x012, $x02x, $x120, $x220, $x221, $xX);
# }

# sub nxxxs{	     # given mat, pat, child gts, get n000, n001, etc.
#   my $mat_gts = shift;		# genotypes string
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
#     my $tr; 
#     #    print STDERR "$c  $m $p  $tr \n";
#     #   $tr .= $c;
#     if ($c eq 'X'  or  $m eq 'X'  or  $p eq 'X') { # missing data case
#       $tr = 'X';
#     } else {			# all 3 gts present. 
#       $tr = ($m < $p)? "$m$p$c" : "$p$m$c";
#     }
#     $n_c{$tr}++;
#   }
#   my @nxxxs = ($n_c{'000'}, $n_c{'001'}, $n_c{'002'},
# 	       $n_c{'010'}, $n_c{'011'}, $n_c{'012'},
# 	       $n_c{'020'}, $n_c{'021'}, $n_c{'022'},
# 	       $n_c{'110'}, $n_c{'111'}, $n_c{'112'},
# 	       $n_c{'120'}, $n_c{'121'}, $n_c{'122'},
# 	       $n_c{'220'}, $n_c{'221'}, $n_c{'222'}, $n_c{'X'}, scalar @ch_gts);
#   return \@nxxxs;
# }

1;
