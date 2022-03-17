#!/usr/bin/env perl
# This script is ued to get UCSC BED file with track lines for being loaded to UCSC or IBV like browsers
use strict;
use warnings;
use File::Basename;

my $usage = "Usage: $0 INFILE OUTFILE [TRACK-NAME]";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $track_name = shift || basename($infile);

my $track_desc = "ITR-Seq reference peak";

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
print OUT "#gffTags\n";
print OUT qq(track name="$track_name" description="$track_desc"\n);

while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $rnames, $score, $strand) = split(/\t/, $line);
	my $name = "$chr:$start-$end"; # use loc as name
	my $id = "$track_name:$name";
# get read count and UMI count
	my %UMI_count;
	my %name_count;
	foreach my $rname (split(/,/, $rnames)) {
		$rname =~ s/\/\d+$//; # remove trailing /1 or /2 suffix as R1 or R2 indications
		my ($UMI) = $rname =~ /UMI:([A-Za-z]+)/;
		$UMI_count{$UMI}++;
		$name_count{$rname}++;
	}
	my $read_count = scalar keys %name_count;
	my $UMI_count = scalar keys %UMI_count;
	my $UMIs = join(',', sort keys %UMI_count);
	my $bedName = qq(ID=$id;Name=$name;ReadCount=$read_count;UMICount=$UMI_count;UMIs=$UMIs;);
	print OUT "$chr\t$start\t$end\t$bedName\t$score\t.\n"; # merged map has no strand
}

close(IN);
close(OUT);
