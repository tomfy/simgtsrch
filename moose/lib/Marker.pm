package Marker;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );
#use Genotypes;

has id => (
	   isa => 'Str',
	   is => 'ro',
	   required => 1,
	  );

has counts => ( # e.g. [500, 11, 400, 13, 100, 1024] counts of 0, 1, 2, in between 0 and 1, and in between 1 and 2, and sum of all 5
	       isa => 'ArrayRef',
	       is => 'rw',
	       default => sub { [0, 0, 0, 0, 0, 0] },
	      );

sub increment_counts{
  my $self = shift;
  my $gt = shift;		# 0, X, 1, x, or 2
  if ($gt eq '0') {
    $self->counts()->[0]++;
  } elsif ($gt eq '1') {
    $self->counts()->[1]++;
  } elsif ($gt eq '2') {
    $self->counts()->[2]++;
  } elsif ($gt eq 'X') {
    $self->counts()->[3]++;
  } elsif($gt eq 'x'){
    $self->counts()->[4]++;
  }else{
    die "Unknown genotype: $gt \n";
  }
  $self->counts()->[5]++;
}

sub hardy_weinberg{ # if hardy-weinberg eq. frequencies, return value close to 1.
  my $self = shift;
  my $n0 = $self->counts()->[0] // 0;
  my $n1 = $self->counts()->[1] // 0;
  my $n2 = $self->counts()->[2] // 0;
  my $ntotal = $n0 + $n1 + $n2;
  my $hw_one = ($ntotal > 0)? ($n0**0.5 + $n2**0.5)/$ntotal**0.5 : -1;
  return $hw_one;
}

1;
