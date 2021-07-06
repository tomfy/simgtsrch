#!/usr/bin/perl -w
use strict;

my $marker_id = shift;

my %markerid_index = ();
my $n_markers = undef;
my @md_counts = ();
my $the_md_count = 0;
my $accession_count = 0;
while(<>){
  next if(/^\s*#/);
  if(/^MARKER\s/){
    my @marker_ids = split(" ", $_);
    shift @marker_ids;
    $n_markers = scalar @marker_ids;
    while(my($i, $mrkrid) = each @marker_ids){
      $markerid_index{$mrkrid} = $i;
    }
  }else{
    my ($accid, $genotypes) = split(" ", $_);
    die if(length $genotypes  !=  $n_markers);
   # for(my $i=0; $i < length $genotypes; $i++){
   #   if(substr($genotypes, $i, 1) == 3){#
   #   $md_counts[$i]++;
    #   }
    $accession_count++;
    $the_md_count++ if(substr($genotypes, $markerid_index{$marker_id}, 1) == 3);
  }

}

print "marker  $marker_id  has $the_md_count  missing_data out of $accession_count.\n";
