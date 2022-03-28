#!/usr/bin/env perl
# This script is ued to get ITR insertion site (of given flanking length) in BED format from dedup-novec alignments
use strict;
use warnings;

my $usage = "Usage: $0 INFILE OUTFILE [--insert-size 20] [--min-softclip 29] [--keep-dup]";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $insert_size = 20;
my $min_softclip = 29;
my $keep_dup = 0;

# parse options
for(my $i = 0; $i < @ARGV; $i++) {
	if($ARGV[$i] eq '--insert-size') {
		$insert_size = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '--min-softclip') {
		$min_softclip = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '--keep-dup') {
		$keep_dup = $ARGV[++$i];
	}
	else {
		print STDERR "Error: unknown option $ARGV[$i]\n";
		exit;
	}
}

if(!($insert_size > 0)) {
	print STDERR "--insert-size must be positive\n";
	exit;
}

if(!($min_softclip > 0)) {
	print STDERR "--min-softclip must be positive\n";
	exit;
}


# open BAM input
my $flags = $keep_dup ? "" : "-F 0x400";
open(IN, "samtools view $flags $infile | ") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read and output
while(my $line = <IN>) {
	chomp $line;
	my ($qname, $flag, $chr, $start, $mapQ, $cigar) = split(/\t/, $line);
	my $strand = ($flag & 0x10) ? '-' : '+';
	my $mate = ($flag & 0x40) ? 1 : 2;
	$start--; # use 0-based start
	my $clip_end;
	if($mate == 1) {
		$clip_end = $strand eq '+' ? 3 : 5;
	}
	else {
		$clip_end = $strand eq '-' ? 3 : 5;
	}

	my $align_len = get_align_len_from_cigar($cigar);
	my $clip_len = get_clip_len_from_cigar($cigar, $clip_end);
	my $end = $start + $align_len;
	if($clip_len >= $min_softclip) {
		my $insert_pos = $mate == 1 ? $end : $start;
		print OUT "$chr\t", ($insert_pos - $insert_size / 2), "\t", ($insert_pos + $insert_size / 2), "\t$qname/$mate\t$mapQ\t$strand\n";
	}
}

close(IN);
close(OUT);

sub get_align_len_from_cigar {
	my $cigar = shift;
	my $align_len = 0;
	while($cigar =~ /(\d+)([MIDNSHPX=])/g) {
    if($2 eq 'M' || $2 eq 'D' || $2 eq 'N' || $2 eq '=' || $2 eq 'X') {
			$align_len += $1;
		}
	}
	return $align_len;
}

sub get_clip_len_from_cigar {
	my ($cigar, $clip_end) = @_;
	my $clip_regex = $clip_end == 3 ? qr/(\d+)S$/ : qr/^(\d+)S/;
	my $clip_len = 0;
	if($cigar =~ /$clip_regex/) {
		$clip_len = $1;
	}
	return $clip_len;
}
