#!/bin/env perl
# This script is used to relabel paired-end FASTQ files to add the UMI information
# UMIs are the 3' end bases of the I2 sequences
use strict;
use warnings;
use constant {
	DEFAULT_UMI_LEN => 8
};

use Getopt::Long;

my $UMI_len = DEFAULT_UMI_LEN;
my $usage = "Usage: $0 -i1 FWD-INFILE -i2 REV-INFILE -o1 FWD-OUTFILE -o2 REV-OUTFILE [-idx INDEX-INFILE] [-l/--len UMI_LEN ($UMI_len)]";
my $in1;
my $in2;
my $out1;
my $out2;
my $idx;

# parse options
GetOptions(
  "i1=s" => \$in1,
	"i2=s" => \$in2,
	"o1=s" => \$out1,
	"o2=s" => \$out2,
	"idx=s" => \$idx,
	"l|len=i" => \$UMI_len
);

# check options
unless(defined $in1 && defined $in2) {
	print STDERR "$usage\n";
	exit;
}

# open inputs
if($in1 =~ /\.gz$/) {
	open(IN1, "zcat $in1 |") || die "Unable to open $in1: $!";
}
elsif($in1 =~ /\.bz2$/) {
	open(IN1, "bzcat $in1 |") || die "Unable to open $in1: $!";
}
else {
	open(IN1, "<$in1") || die "Unable to open $in1: $!";
}

if($in2 =~ /\.gz$/) {
	open(IN2, "zcat $in2 |") || die "Unable to open $in2 $!";
}
elsif($in2 =~ /\.bz2$/) {
	open(IN2, "bzcat $in2 |") || die "Unable to open $in2: $!";
}
else {
	open(IN2, "<$in2") || die "Unable to open $in2: $!";
}

if(defined $idx) {
	if($idx =~ /\.gz$/) {
		open(IDX, "zcat $idx |") || die "Unable to open $idx: $!";
	}
	elsif($idx =~ /\.bz2$/) {
		open(IDX, "bzcat $idx |") || die "Unable to open $idx: $!";
	}
	else {
		open(IDX, "<$idx") || die "Unable to open $idx: $!";
	}
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

# Scan R1/R2 and optionall I1 to get UMIs
while(!eof(IN1) && !eof(IN2)) {
	my $def1 = <IN1>; chomp $def1;
	my $def2 = <IN2>; chomp $def2;
	my $seq1 = <IN1>;
	my $seq2 = <IN2>;
	my $sep1 = <IN1>;
	my $sep2 = <IN2>;
	my $qual1 = <IN1>;
	my $qual2 = <IN2>;

	my ($name1, $desc1) = $def1 =~ /^@(\S+)(.*)/;
	my ($name2, $desc2) = $def2 =~ /^@(\S+)(.*)/;

	my ($UMI1, $UMI2);
  if(defined $idx) {
		my $defi = <IDX>;
	  my $seqi = <IDX>; chomp $seqi;
		my $sepi = <IDX>;
		my $quali = <IDX>;
		$UMI1 = $UMI2 = substr($seqi, -$UMI_len); # UMI near the 3' of the I2 read
	}
	else {
		my @name_fields1 = split(/:/, $name1);
		my @name_fields2 = split(/:/, $name2);
		($UMI1) = pop @name_fields1;
		($UMI2) = pop @name_fields2;
		$UMI1 =~ s/[^ATCGNatcgn]//g;
		$UMI2 =~ s/[^ATCGNatcgn]//g;
	}

# update def lines
	$def1 = "@" . $name1 .":UMI:$UMI1" . $desc1;
	$def2 = "@" . $name2 .":UMI:$UMI2" . $desc2;

# write outputs
	print OUT1 "$def1\n", $seq1, $sep1, $qual1;
	print OUT2 "$def2\n", $seq2, $sep2, $qual2;
}

close(IN1);
close(IN2);
close(IDX);
close(OUT1);
close(OUT2);
