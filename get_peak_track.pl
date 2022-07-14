#!/usr/bin/env perl
# This script is ued to get UCSC BED file with track lines for being loaded to UCSC or IBV like browsers
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $sample_desc = "ITR-Seq ref mapped peak";
my $usage = "Usage: $0 INFILE OUTFILE [--name SAMPLE-NAME (INFILE)] [--desc SAMPLE-DESC ($sample_desc)]";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $sample_name = basename($infile);

GetOptions(
"name=s" => \$sample_name,
"desc=s" => \$sample_desc)
or die "Error in command line arguments, usage: $usage";

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
print OUT "#gffTags\n";
print OUT qq(track name="$sample_name-ITR-peak" description="$sample_desc"\n);

while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $rnames, $score, $peak_strand_counts) = split(/\t/, $line);
	my $name = "$chr:$start-$end"; # use loc as name
	my $id = "$sample_name:$name";
# get read count and UMI count
	my %name_count;
	foreach my $rname (split(/,/, $rnames)) {
		$rname =~ s/\/\d+$//; # remove trailing /1 or /2 suffix as R1 or R2 indications
		$name_count{$rname}++;
	}
	my $read_count = scalar keys %name_count;
	my $bedName = qq(ID=$id;Name=$name;ReadCount=$read_count;ReadStrandCounts=$peak_strand_counts;);
	print OUT "$chr\t$start\t$end\t$bedName\t$score\t.\n"; # merged map has no strand
}

close(IN);
close(OUT);
