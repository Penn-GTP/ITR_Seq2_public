#!/bin/env perl
# This script is used to get uniq IDs using UMI tags and alignment locations given reference alignment file
use strict;
use warnings;

my $usage = "Usage: $0 SAM-INFILE OUTFILE";
my $in = shift or die $usage;
my $out = shift or die $usage;

# open input
open(IN, $in) || die "Unable to open $in: $!";

# open outputs
open(OUT, ">$out") || die "Unable to write to $out: $!";

# Scan alignments get unique mapped loc
my %uid2mapQ; # uid2best mapQ
my %uid2qname; # uid2qname
while(my $line = <IN>) {
	chomp $line;
	my ($qname, $flag, $rname, $start, $mapQ, $cigar) = split(/\t/, $line);
	my ($UMI) = $qname =~ /(UMI:\w+)/;
	my $align_len = get_align_len_by_cigar($cigar);
	my $end = $start + $align_len - 1;
	my $uid = "$UMI:$rname:$start:$end";
	if(!exists $uid2mapQ{$uid} || $mapQ > $uid2mapQ{$uid}) {
		$uid2mapQ{$uid} = $mapQ;
		$uid2qname{$uid} = $qname;
	}
}

# output
while(my ($uid, $qname) = each %uid2qname) {
	print OUT "$qname\n";
}

sub get_align_len_by_cigar {
	my $cigar = shift;
	my $len = 0;
  while($cigar =~ /(\d+)([MIDNSHP=X])/g) {
    my ($oplen, $op) = ($1, $2);
    if($op eq 'M' || $op eq 'D' || $op eq 'N' || $op eq '=' || $op eq 'X') {
      $len += $oplen;
    }
  }
  return $len;
}

