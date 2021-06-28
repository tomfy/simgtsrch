#!/usr/bin/perl -w
use strict;

my %parentidpair_offspringids = (); # key: ordered id pair of parents, value string of offspring ids

while (<>) {			# read pedigree_test output file
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($accid, $matid, $patid) = @cols[0,2,3];
  if ($matid gt $patid) {
    $parentidpair_offspringids{"$matid $patid"} .= "$accid ";
  } else {
    $parentidpair_offspringids{"$patid $matid"} .= "$accid ";
  }
}

while (my($k, $v) = each %parentidpair_offspringids) {
  my @accids = split(" ", $v);
  print "$k  ", scalar @accids, "  $v\n";
}
