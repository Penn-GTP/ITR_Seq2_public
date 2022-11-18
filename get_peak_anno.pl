#!/usr/bin/env perl
# This script is ued to get peak overlap/flanking annotation from anno BED files
use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any none);

my @inc_types;
my @exc_types = qw(chromosome region biological_region mRNA exon pseudogenic_transcript);

my $options = qq([OPTIONS]
OPTIONS:
  --include-type [TYPE1[,TYPE2]]: include given types in the GFF file as the annotation source, multiple values allowed separated by comma or by given multiple times
  --exclude-type [TYPE3[,TYPE4]]: exclude given types in the GFF file as the annotation source, multiple values allowed separated by comma or by given multiple times; [default: ) . join(",", @exc_types) . qq(]
);

my $usage = "Usage: $0 GFFFILE INFILE OUTFILE $options";

my $gff = shift or die $usage;
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

GetOptions(
"include-type=s" => \@inc_types,
"exclude-type=s" => \@exc_types)
or die "Error in command line arguments, usage: $usage";

# expand potential comma-separated opts
$gff =~ s/,/ /g;
@inc_types = split(/,/, join(',', @inc_types));
@exc_types = split(/,/, join(',', @exc_types));

open(IN, "bedtools closest -a $infile -b $gff -bed -D ref -nonamecheck |") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
print OUT "chrom\tstart\tend\tname\tscore\tstrand\toverlap_details\n";

# read in peak anno
my @peak_names;
my %peak2info;
my %peak2detail;
my %peak2seen;

while(my $line = <IN>) {
	chomp $line;
	my @fields = split(/\t/, $line);
	if(@fields > 16) { # multiple -b options in action
		splice(@fields, 6, 1); # remove the fileID field
	}
	my ($chr, $start, $end, $name, $score, $strand, $over_chr, $over_src, $over_type, $over_start, $over_end, $over_score, $over_strand, $over_frame, $over_attr, $over_dist) = @fields;
	next unless( (any { $over_type eq $_ } @inc_types) || (none { $over_type eq $_ } @exc_types) );

	my $info = "$chr\t$start\t$end\t$name\t$score\t$strand";
  push(@peak_names, $name) if(!exists $peak2info{$name});
	$peak2info{$name} = $info;
	my $over_orient = $over_dist == 0 ? "overlap" : $over_dist < 0 ? "upstream" : "downstream";
	my $tss_dist = $over_attr =~ /Parent=/ ? 'NA' : get_tss_dist($chr, $start, $end, $strand, $over_chr, $over_start, $over_end, $over_strand);
	if(!exists $peak2seen{"$name:$over_type:$over_chr:$over_start-$over_end:$over_strand"}) {
		push(@{$peak2detail{$name}},"Source=$over_src;Type=$over_type;Locus=$over_chr:$over_start-$over_end:$over_strand;Orient=$over_orient;Dist=$over_dist;TssDist=$tss_dist;$over_attr"); # prepend additional information
	}
	$peak2seen{"$name:$over_type:$over_chr:$over_start-$over_end:$over_strand"}++;
}

# output
foreach my $name (@peak_names) {
	my $info = $peak2info{$name};
	my $details = join(' // ', @{$peak2detail{$name}});

	print OUT "$info\t$details\n";
}

close(IN);
close(OUT);

sub get_tss_dist {
	my ($chr1, $start1, $end1, $str1, $chr2, $start2, $end2, $str2) = @_;
	if($chr1 ne $chr2) {
		return 'Inf';
	}
	my $tss = $str2 eq '+' ? $start2: $end2;
	my $dist = 0;

	if($start1 < $tss && $tss <= $end1) {
		$dist = 0;
	}
	elsif($end1 <= $tss) {
		$dist = $tss - $end1;
	}
	else {
		$dist = $tss - $start1 - 1;
	}
	return $dist;
}
