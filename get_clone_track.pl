#!/usr/bin/env perl
# This script is ued to get UCSC BED file with track lines for being loaded to UCSC or IBV like browsers
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $track_desc = "ITR-Seq ref mapped clone";
my $usage = "Usage: $0 INFILE OUTFILE [--name TRACK-NAME (INFILE)] [--desc TRACK-DESC ($track_desc)]";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $track_name = basename($infile);

GetOptions(
"name=s" => \$track_name,
"desc=s" => \$track_desc)
or die "Error in command line arguments, usage: $usage";

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
print OUT "#gffTags\n";
print OUT qq(track name="$track_name" description="$track_desc"\n);

while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $rnames, $score, $loc_counts) = split(/\t/, $line);
	my $name = "$chr:$start-$end"; # use loc as name
	my $id = "$track_name:$name";
# get clone loc count
  my $read_count = scalar split(/,/, $rnames);
	my $loc_count = scalar split(/,/, $loc_counts);
	my $bedName = qq(ID=$id;Name=$name;ReadCount=$read_count;LocCount=$loc_count;);
	print OUT "$chr\t$start\t$end\t$bedName\t$score\t.\n"; # merged loc has no strand
}

close(IN);
close(OUT);
