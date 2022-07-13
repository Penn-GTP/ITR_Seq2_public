#!/usr/bin/env perl
# This script is ued to get ITR clones (unique clone sites) from merged peaks
use strict;
use warnings;
use Getopt::Long;

my $min_loc = 2;

my $usage = "Usage: $0 PEAK-INFILE ALIGN-INFILE OUTFILE [--min-loc INT ($min_loc)]";
my $peak_infile = shift or die $usage;
my $aln_infile = shift or die $usage;
my $outfile = shift or die $usage;

my $samtools = 'samtools';

GetOptions(
"min-loc=i" => \$min_loc)
or die "Error in command line arguments, usage: $usage";

if(!($min_loc > 0)) {
	print STDERR "--min-loc must be greater than 0\n";
	exit;
}

open(PEAK, "<$peak_infile") || die "Unable to open $peak_infile: $!";
if($aln_infile =~ /\.sam$/) {
	open(ALN, "<$aln_infile") || die "Unable to open $aln_infile: $!";
}
elsif($aln_infile =~ /\.bam$/) {
	open(ALN, "samtools view $aln_infile |") || die "Unable to open $aln_infile: $!";
}
else {
	print STDERR "Unable to open $aln_infile, unknown format\n";
	exit;
}
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# Read in R1 read loc
my %name2loc;
while(my $line = <ALN>) {
	next if($line =~ /^@/);
	chomp $line;
	my ($name, $flag, $chr, $loc) = split(/\t/, $line);
	$name2loc{$name} = "$chr:$loc";
}

# Scan peaks and output
while(my $line = <PEAK>) {
	chomp $line;
	my ($chr, $start, $end, $names, $score, $strands) = split(/\t/, $line);
	my @name_loc;
	my %loc2count;

	foreach my $name (split(/,/, $names)) {
		$name =~ s/\/(\d)$//;
		my $loc = $name2loc{$name};
		push(@name_loc, "$name:$loc");
		$loc2count{$name2loc{$name}}++;
	}
	#my $num_loc = scalar keys %loc2count;
	my $name_locs = join(",", @name_loc);
	my $clone_count = scalar keys %loc2count;
	my $clone_loci = join(",", map { "$_:$loc2count{$_}" } sort keys %loc2count);

# output
	if($clone_count >= $min_loc) {
		print OUT "$chr\t$start\t$end\t$name_locs\t$score\t$clone_loci\n";
	}
}

close(PEAK);
close(ALN);
close(OUT);
