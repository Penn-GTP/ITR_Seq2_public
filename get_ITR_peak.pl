#!/usr/bin/env perl
# This script is ued to get ITR insertion site (of given flanking length) from dedup-novec alignments
use strict;
use warnings;

my $usage = "Usage: $0 INFILE OUTFILE [STRAND]";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $keep_strand = shift || 0;

if(!($keep_strand >= 0)) {
	print STDERR "STRAND must be non-negative\n";
	exit;
}

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $names, $score, $strands) = split(/\t/, $line);
	my %strand2UMI;
	my @names = split(/,/, $names);
	my @strands = split(/,/, $strands);
	if(@names != @strands) {
		print STDERR "number of fields in names and strands don't match\n";
		exit;
	}
	for(my $i = 0; $i < @names; $i++) {
		if($names[$i] =~ /\/1$/) { # this is a forward read, where UMI exists
			my ($UMI) = $names[$i] =~ /(UMI:[A-Za-z]+)/;
			$strand2UMI{$strands[$i]}{$UMI}++;
		}
	}

	my $strand = 0;
  foreach my $str (sort keys %strand2UMI) {
		my $rstr = $str eq '+' ? 0x1 : $str eq '-' ? 0x2 : 0;
		$strand |= $rstr;
	}

# filter input
	if($keep_strand == 0 || ($strand & $keep_strand) == $keep_strand) {
		print OUT "$chr\t$start\t$end\t$names\t$score\t.\n"; # merged map has no strand
	}
}

close(IN);
close(OUT);
