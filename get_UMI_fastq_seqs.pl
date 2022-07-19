#!/bin/env perl
# This script is used to relabel paired-end FASTQ files to add the UMI information
# UMIs are the 3' end bases of the I2 sequences
use strict;
use warnings;
use IO::Zlib;
use constant {
	DEFAULT_UMI_LEN => 8
	};

my $usage = "Usage: $0 -i1 FWD-INFILE -i2 REV-INFILE -idx INDEX-INFILE -o1 FWD-OUTFILE -o2 REV-OUTFILE [-l --len UMI-LENGTH]";
my $in1;
my $in2;
my $idx;
my $out1;
my $out2;
my $UMI_len = DEFAULT_UMI_LEN;

# parse options
for(my $i = 0; $i < @ARGV; $i++) {
	if($ARGV[$i] eq '-i1') {
		$in1 = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '-i2') {
		$in2 = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '-idx') {
		$idx = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '-o1') {
		$out1 = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '-o2') {
		$out2 = $ARGV[++$i];
	}
	elsif($ARGV[$i] eq '-l' || $ARGV[$i] eq '--len') {
		$UMI_len = $ARGV[++$i];
	}
	else {
		print STDERR "Warning: unknown option: '", $ARGV[$i], "', ignored\n";
		exit;
	}
}

# check options
unless(defined $in1 && defined $in2 && $UMI_len > 0) {
	print STDERR "$usage\n";
	exit;
}

# open inputs
my ($seqI1, $seqI2, $seqIdx);
my ($seqO1, $seqO2);

if($in1 =~ /\.gz$/) {
	$seqI1 = new IO::Zlib($in1, "rb") || die "Unable to open $in1: $!";
}
else {
	open($seqI1, "<$in1") || die "Unable to open $in1: $!";
}

if($in2 =~ /\.gz$/) {
	$seqI2 = new IO::Zlib($in2, "rb") || die "Unable to open $in2 $!";
}
else {
	open($seqI2, "<$in2") || die "Unable to open $in2: $!";
}

if($idx =~ /\.gz$/) {
	$seqIdx = new IO::Zlib($idx, "rb") || die "Unable to open $idx: $!";
}
else {
	open($seqIdx, "<$idx") || die "Unable to open $idx: $!";
}

# open outputs
if($out1 =~ /\.gz$/) {
	$seqO1 = new IO::Zlib($out1, "wb") || die "Unable to write to $out1: $!";
}
else {
	open($seqO1, ">$out1") || die "Unable to write to $out1: $!";
}

if($out2 =~ /\.gz$/) {
	$seqO2 = new IO::Zlib($out2, "wb") || die "Unable to write to $out2: $!";
}
else {
	open($seqO2, ">$out2") || die "Unable to write to $out2: $!";
}

# Scan index file and get UMIs using 5' bases
my $n_processed = 0;
while(my $line = <$seqIdx>) {
	chomp $line;
	if($n_processed % 4 == 0) { # def line
		#chomp $line1;
		#chomp $line2;
		my $def1 = <$seqI1>;
		my $def2 = <$seqI2>;

		my $seq1 = <$seqI1>; chomp $seq1;
		my $seq2 = <$seqI2>; chomp $seq2;
    my $seqi = <$seqIdx>; chomp $seqi;

		my $sep1 = <$seqI1>;
		my $sep2 = <$seqI2>;

		my $qual1 = <$seqI1>;
		my $qual2 = <$seqI2>;
    my $quali = <$seqIdx>;

		my $UMI = substr($seqi, -$UMI_len); # UMI near the 3' of the I2 read
# update def lines
		$def1 =~ s/^@\S+/$&:UMI:$UMI/;
		$def2 =~ s/^@\S+/$&:UMI:$UMI/;
# write outputs
		print $seqO1 $def1, "$seq1\n", $sep1, $qual1;
		print $seqO2 $def2, "$seq2\n", $sep2, $qual2;
	}
	$n_processed += 4;
}

$seqI1->close();
$seqI2->close();
$seqIdx->close();
$seqO1->close();
$seqO2->close();
