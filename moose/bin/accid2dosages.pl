#!/usr/bin/perl -w
use strict;

my $accid = shift;

while(<>){
	next if(/^\s*#/);
	my @cols = split(" ", $_);
	my $aid = shift @cols;
	if($aid eq $accid){
		print "# $accid \n";
		print join("\n", @cols), "\n";
	exit;
	}
}

