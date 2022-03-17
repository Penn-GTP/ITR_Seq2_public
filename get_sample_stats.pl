#!/bin/env perl
# Prepare sh script for filtering reference mapping files
our $VERSION = v1.1;

use strict;
use warnings;
use lib '/project/gtplab/pipeline/ITR_Seq';
use MiSeqITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
#my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my @headers = qw(sample_name total_read trimmed_read ref_mapped vec_mapped ref_novec_mapped peak_count peak_read_count peak_UMI_count target_count target_read_count target_UMI_count);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new MiSeqITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
	print STDERR "Error: BASE_DIR $BASE_DIR not exists\n";
	exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $DEMUX_DIR)) {
	  print STDERR "Error: DEMUX $DEMUX_DIR not exists\n";
		  exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT join("\t", @headers), "\n";

foreach my $sample ($design->get_sample_names()) {
	print STDERR "gathering stats for $sample\n";
# get total read
  my $total_read;
	{
		my $in = $design->sample_opt($sample, 'fastq_R1');
		$total_read = $in =~ /\.gz$/ ? `zcat $DEMUX_DIR/$in | wc -l` : `cat $DEMUX_DIR/$in | wc -l`;
		chomp $total_read;
		$total_read /= 4;
	}

# get trimmed read
  my $trimmed_read;
  {
		my $in = $design->get_sample_fwd_ITRtrim_file($sample);
		$trimmed_read = $in =~ /\.gz$/ ? `zcat $WORK_DIR/$in | wc -l` : `cat $DEMUX_DIR/$in | wc -l`;
		chomp $trimmed_read;
		$trimmed_read /= 4;
	}

# get ref mapped
  my $ref_mapped;
	{
		my $in = $design->get_sample_ref_filtered_file($sample);
		$ref_mapped = `samtools view $WORK_DIR/$in | wc -l`;
		chomp $ref_mapped;
	}

# get vec mapped
	my $vec_mapped;
	{
		my $in = $design->get_sample_vec_filtered_file($sample);
		$vec_mapped = `samtools view $WORK_DIR/$in | wc -l`;
		chomp $vec_mapped;
	}

# get novec mapped
  my $novec_mapped;
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		$novec_mapped = `samtools view $WORK_DIR/$in | wc -l`;
		chomp $novec_mapped;
	}


# get peak info
  my ($peak_count, $peak_read_count, $peak_UMI_count) = (0, 0, 0);
	{
		my $in = $design->get_sample_ref_filtered_peak($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			chomp $line;
			$peak_count++;
			my ($rnames) = (split(/\t/, $line))[3];
			my %read_count;
			my %UMI_count;
			foreach my $rname (split(/,/, $rnames)) {
				if($rname =~ /\/1$/) {
					my ($UMI) = $rname =~ /:UMI:(\w+)/;
					$UMI_count{$UMI}++;
				}
				$rname =~ s/\/\d+$//; # remove trailing /1 or /2
				$read_count{$rname}++;
			}
			$peak_read_count += scalar keys %read_count;
			$peak_UMI_count += scalar keys %UMI_count; # count uniq UMI per perak
		}
		close(BED);
	}

# get target info
	my ($target_count, $target_read_count, $target_UMI_count) = (0, 0, 0);
	my $target_file = $design->sample_opt($sample, 'target_file');
	if(!-e $target_file) {
		($target_count, $target_read_count, $target_UMI_count) = qw(NA NA NA);
	}
	else {
		my $in = $design->get_sample_ref_filtered_peak($sample);
		if(-s "$BASE_DIR/$in") { # non-empty peaks found
			open(BED, "$bedtools intersect -a $BASE_DIR/$in -b $target_file -wo |") || die "Unable to open $samtools intersect: $!";
			while(my $line = <BED>) {
				chomp $line;
				$target_count++;
				my ($rnames) = (split(/\t/, $line))[3];
				my %read_count;
				my %UMI_count;
				foreach my $rname (split(/,/, $rnames)) {
					if($rname =~ /\/1$/) {
						my ($UMI) = $rname =~ /:UMI:(\w+)/;
						$UMI_count{$UMI}++;
					}
					$rname =~ s/\/\d+//; # remove tailing /1 or /2
						$read_count{$rname}++;
				}
				$target_read_count += scalar keys %read_count;
				$target_UMI_count += scalar keys %UMI_count;
			}
			close(BED);
		}
	}

# output
  print OUT "$sample\t$total_read\t$trimmed_read\t$ref_mapped\t$vec_mapped\t$novec_mapped\t$peak_count\t$peak_read_count\t$peak_UMI_count\t$target_count\t$target_read_count\t$target_UMI_count\n";
}

close(OUT);
