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
	if($line =~ /^#/) {
		print OUT $line;
		next;
	}
	chomp $line;
	my @fields = split(/\t/, $line);
	my ($name) = $fields[3] =~ /Name=([^;]+)/;
	my $id = "$sample_name:$name";
	$fields[3] = "ID=$id;" . $fields[3];

	print OUT join("\t", @fields), "\n";
}

close(IN);
close(OUT);
