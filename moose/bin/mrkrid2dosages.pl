#!/usr/bin/perl -w
use strict;

my $marker_id = shift;

my %markerid_index = ();
my $n_markers = undef;
my @md_counts = ();
my $the_md_count = 0;
my $accession_count = 0;
my $the_marker_index;
while(<>){
  next if(/^\s*#/);
  if(/^MARKER\s/){
    my @marker_ids = split(" ", $_);
    shift @marker_ids;
    $n_markers = scalar @marker_ids;
    while(my($i, $mrkrid) = each @marker_ids){
      $markerid_index{$mrkrid} = $i;
    }
    $the_marker_index = $markerid_index{$marker_id};
print "# marker id: $marker_id  $the_marker_index \n"; 
exit;
  }else{
	  my @dosages = split(" ", $_);
	  my $accid = shift @dosages;
    die print scalar @dosages, "  $n_markers" if(scalar @dosages  !=  $n_markers);
   # for(my $i=0; $i < length $genotypes; $i++){
   #   if(substr($genotypes, $i, 1) == 3){#
   #   $md_counts[$i]++;
    #   }
    $accession_count++;
    #$the_md_count++ if(substr($genotypes, $markerid_index{$marker_id}, 1) == 3);
	print "$accid  ", $dosages[$the_marker_index], "\n";  
}

}

# print "marker  $marker_id  has $the_md_count  missing_data out of $accession_count.\n";
