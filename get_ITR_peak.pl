#!/usr/bin/env perl
# This script is ued to get ITR insertion site (of given flanking length) from dedup-novec alignments
use strict;
use warnings;

my $usage = "Usage: $0 INFILE OUTFILE";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $names, $score, $strands) = split(/\t/, $line);
	my @read_count = (0, 0, 0, 0); # R1 +,-|R2 +,-

	my @names = split(/,/, $names);
	my @strands = split(/,/, $strands);
	if(@names != @strands) {
		print STDERR "number of fields in names and strands don't match\n";
		exit;
	}
	for(my $i = 0; $i < @names; $i++) {
		my ($mate) = $names[$i] =~ /\/(\d)$/;
		my $mate_idx = $mate == 1 ? 0 : 2;
		my $strand_idx = $strands[$i] eq '+' ? 0 : 1;
		$read_count[$mate_idx + $strand_idx]++;
	}

	my ($fwd_plus, $fwd_minus, $rev_plus, $rev_minus) = @read_count;

# filter input
	if($fwd_plus > 0 && $fwd_minus > 0 || $rev_plus > 0 && $rev_minus > 0) {
	  print OUT "$chr\t$start\t$end\t$names\t$score\t$fwd_plus,$fwd_minus|$rev_plus,$rev_minus\n";
	}
}

close(IN);
close(OUT);
