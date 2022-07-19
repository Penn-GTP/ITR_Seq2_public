#!/bin/env perl
# This script is used to relabel paired-end FASTQ files to add the UMI information
# UMIs are the 3' end bases of the I2 sequences
use strict;
use warnings;
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
if($in1 =~ /\.gz$/) {
	open(IN1, "zcat $in1 |") || die "Unable to open $in1: $!";
}
else {
	open(IN1, "<$in1") || die "Unable to open $in1: $!";
}

if($in2 =~ /\.gz$/) {
	open(IN2, "zcat $in2 |") || die "Unable to open $in2 $!";
}
else {
	open(IN2, "<$in2") || die "Unable to open $in2: $!";
}

if($idx =~ /\.gz$/) {
	open(IDX, "zcat $idx |") || die "Unable to open $idx: $!";
}
else {
	open(IDX, "<$idx") || die "Unable to open $idx: $!";
}

# open outputs
if($out1 =~ /\.gz$/) {
	open(OUT1, "| gzip > $out1") || die "Unable to write to $out1: $!";
}
else {
	open(OUT1, ">$out1") || die "Unable to write to $out1: $!";
}

if($out2 =~ /\.gz$/) {
	open(OUT2, "| gzip > $out2") || die "Unable to write to $out2: $!";
}
else {
	open(OUT2, ">$out2") || die "Unable to write to $out2: $!";
}

# Scan index file and get UMIs using 5' bases
my $n_processed = 0;
while(my $line = <IDX>) {
	chomp $line;
	if($n_processed % 4 == 0) { # def line
		#chomp $line1;
		#chomp $line2;
		my $def1 = <IN1>;
		my $def2 = <IN2>;

		my $seq1 = <IN1>; chomp $seq1;
		my $seq2 = <IN2>; chomp $seq2;
    my $seqi = <IDX>; chomp $seqi;

		my $sep1 = <IN1>;
		my $sep2 = <IN2>;

		my $qual1 = <IN1>;
		my $qual2 = <IN2>;
    <IDX>;

		my $UMI = substr($seqi, -$UMI_len); # UMI near the 3' of the I2 read
# update def lines
		$def1 =~ s/^@\S+/$&:UMI:$UMI/;
		$def2 =~ s/^@\S+/$&:UMI:$UMI/;
# write outputs
		print OUT1 $def1, "$seq1\n", $sep1, $qual1;
		print OUT2 $def2, "$seq2\n", $sep2, $qual2;
	}
	$n_processed += 4;
}

close(IN1);
close(IN2);
close(IDX);
close(OUT1);
close(OUT2);
