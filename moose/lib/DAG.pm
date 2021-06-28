package DAG;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

has id_node => (
		isa => 'HashRef[Maybe[Str]]', # key: id (string), value: Node object
		is => 'rw',
		default => sub{ {} },
	       );

has childless_ids => (
		      isa => 'HashRef[Maybe[Str]]', #  key: id (string), value: Node Obj
		      is => 'rw',
		      default => sub{ {} },
		     );

has parentless_ids => (
		       isa => 'HashRef[Maybe[Str]]', #  key: id (string), value: Node Obj
		       is => 'rw',
		       default => sub{ {} },
		      );

sub add_node{
  my $self = shift;
  my $node = shift;
  my $node_id = $node->id();
  $self->id_node()->{$node_id} = $node;
  if ($node->offspring() eq '') {
    $self->childless_ids()->{$node_id} = $node;
  }
  if (!defined $node->female_parent()  and  !defined $node->male_parent()) {
    $self->parentless_ids()->{$node_id} = $node;
  }
}

sub is_it_acyclic{
  my $self = shift;
  for my $childless_node ( map($self->childless_ids()->{$_}, sort keys %{$self->childless_ids()} ) ){
    my $res = $childless_node->ancestors_acyclic( {} );
    print STDERR "node ", $childless_node->id(), " and ancestors acyclic?  $res \n\n";
    return 0 if($res == 0);
  }
  return 1;
}

# sub is_it_acyclic_down{
#   my $self = shift;
#   for my $parentless_node ( map($self->parentless_ids()->{$_}, sort keys %{$self->parentless_ids()} ) ){
#     my $res = $parentless_node->descendants_acyclic( {} );
#     return 0 if($res == 0);
#   }
#   return 1;
# }

1;
