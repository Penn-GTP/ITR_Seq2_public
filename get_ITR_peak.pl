#!/usr/bin/env perl
# This script is ued to get ITR peaks (frequent insertion sites) from merged ITR peaks
use strict;
use warnings;
use Getopt::Long;

my $keep_strand = 3;

my $usage = "Usage: $0 INFILE OUTFILE [--keep-strand 0|1|2|3 ($keep_strand)]";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

GetOptions(
"keep-strand=i" => \$keep_strand)
or die "Error in command line arguments, usage: $usage";

if(!(0 <= $keep_strand && $keep_strand <= 3)) {
	print STDERR "--keep-strand must be one of 0 (no-requirement), 1 (matches plus strand), 2 (matches minus strand) or 3 (matches both strands)\n";
	exit;
}

open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
my $id = 0;
while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $rnames, $score, $strands) = split(/\t/, $line);

  my %name_count;
  foreach my $rname (split(/,/, $rnames)) {
    $rname =~ s/\/\d+$//; # remove trailing /1 or /2 suffix as R1 or R2 indications
    $name_count{$rname}++;
  }
  my $read_count = scalar keys %name_count;

	my @read_count = (0, 0, 0, 0); # R1 +,-|R2 +,-
	my $strand_bit = 0; # strand bit indicating which strand have hits

	my @rnames = split(/,/, $rnames);
	my @strands = split(/,/, $strands);
	if(@rnames != @strands) {
		print STDERR "number of fields in names and strands don't match\n";
		exit;
	}
	for(my $i = 0; $i < @rnames; $i++) {
		my ($mate) = $rnames[$i] =~ /\/(\d)$/;
		my $mate_idx = $mate == 1 ? 0 : 2;
		my $strand_idx = $strands[$i] eq '+' ? 0 : 1;
		$read_count[$mate_idx + $strand_idx]++;
		$strand_bit |= $strand_idx + 1;
	}

	my ($fwd_plus, $fwd_minus, $rev_plus, $rev_minus) = @read_count;


# filter input
	if($keep_strand == 0 || ($strand_bit & $keep_strand) >= $keep_strand) {
		my $name = "peak" . (++$id);
		my $peak_name = qq(Name=$name;ReadCount=$read_count;ReadStrandCounts=$fwd_plus,$fwd_minus|$rev_plus,$rev_minus;);
	  print OUT "$chr\t$start\t$end\t$peak_name\t$score\t.\n";
	}
}

close(IN);
close(OUT);
