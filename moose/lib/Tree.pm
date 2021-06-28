package Tree;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

has Root => (
		isa => 'Maybe[Node]', # Node object
		is => 'rw',
		default => undef,
	       );
