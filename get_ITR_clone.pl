#!/usr/bin/env perl
# This script is ued to get ITR clones (unique clone locus) from merged peaks
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
	$name2loc{$name} = "$chr|$loc";
}

# Scan peaks and output
my $id = 0;
while(my $line = <PEAK>) {
	chomp $line;
	my ($chr, $start, $end, $rnames, $score, $strands) = split(/\t/, $line);

	my @name_loc;
	my %loc2count;

	foreach my $rname (split(/,/, $rnames)) {
		$rname =~ s/\/(\d)$//;
		my $loc = $name2loc{$rname};
		push(@name_loc, "$rname|$loc");
		$loc2count{$name2loc{$rname}}++;
	}
	my $num_umi = scalar @name_loc;
	my $num_loc = scalar keys %loc2count;
	my $loc_freq = join(",", map { "$_:$loc2count{$_}" } sort keys %loc2count);


# output
	if($num_loc >= $min_loc) {
		my $name = "loc" . (++$id);
		my $clone_name = qq(Name=$name;UMICount=$num_umi;LocCount=$num_loc;LocFreq=$loc_freq;);
		print OUT "$chr\t$start\t$end\t$clone_name\t$score\t.\n";
	}
}

close(PEAK);
close(ALN);
close(OUT);
