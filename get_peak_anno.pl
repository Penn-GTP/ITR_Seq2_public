#!/usr/bin/env perl
# This script is ued to get peak annotation details from anno BED files
use strict;
use warnings;

my $usage = "Usage: $0 GFF-FILE INFILE OUTFILE";
my $gff = shift or die $usage;
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

open(IN, "bedtools intersect -a $infile -b $gff -wao |") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
print OUT "chrom\tstart\tend\tname\tscore\tstrand\toverlap_details\n";

# read in peak anno
my @peak_names;
my %peak2info;
my %peak2detail;

while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $name, $score, $strand, $over_chr, $over_src, $over_type, $over_start, $over_end, $over_score, $over_strand, $over_frame, $over_attr, $over_len) = split(/\t/, $line);
	my $info = "$chr\t$start\t$end\t$name\t$score\t$strand";
  push(@peak_names, $name) if(!exists $peak2info{$name});
	$peak2info{$name} = $info;
	if($over_len > 0 && $over_attr !~ /Parent=/) { # overlap exists and is a top-level feature
		push(@{$peak2detail{$name}},"Source=$over_src;Type=$over_type;Locus=$over_chr:$over_start-$over_end:$over_strand;Overlap=$over_len;$over_attr"); # prepend additional information
	}
}

# output
foreach my $name (@peak_names) {
	my $info = $peak2info{$name};
	my $details = exists $peak2detail{$name} ? join(' // ', @{$peak2detail{$name}}) : '';
	print OUT "$info\t$details\n";
}

close(IN);
close(OUT);
