#!/usr/bin/env perl
# This script is ued to get UCSC BED file with track lines for being loaded to UCSC or IBV like browsers
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $sample_desc = "ITR-Seq ref mapped clone";
my $usage = "Usage: $0 INFILE TRACK-OUTFILE INFO-OUTFILE [--name SAMPLE-NAME (INFILE)] [--desc SAMPLE-DESC ($sample_desc)]";
my $infile = shift or die $usage;
my $track_outfile = shift or die $usage;
my $info_outfile = shift or die $usage;
my $sample_name = basename($infile);

GetOptions(
"name=s" => \$sample_name,
"desc=s" => \$sample_desc)
or die "Error in command line arguments, usage: $usage";

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(TRACK, ">$track_outfile") || die "Unable to write to $track_outfile: $!";
open(INFO, ">$info_outfile") || die "Unable to write to $info_outfile: $!";

# read and output
print TRACK "#gffTags\n";
print TRACK qq(track name="$sample_name-ITR_clone" description="$sample_desc"\n);

print INFO "sample_name\tclone_id\tUMI_count\tloc_count\n";

while(my $line = <IN>) {
	chomp $line;
	print STDERR "$line\n";
	my ($chr, $start, $end, $rnames, $score, $loc_counts) = split(/\t/, $line);
	my $name = "$chr:$start-$end"; # use loc as name
	my $id = "$sample_name:$name";
# get clone loc count
  my $umi_count = scalar split(/,/, $rnames);
	my $loc_count = scalar split(/,/, $loc_counts);
	my $bedName = qq(ID=$id;Name=$name;UMICount=$umi_count;LocCount=$loc_count;);
	print TRACK "$chr\t$start\t$end\t$bedName\t$score\t.\n"; # merged loc has no strand
	print INFO "$sample_name\t$id\t$umi_count\t$loc_count\n";
}

close(IN);
close(TRACK);
close(INFO);
