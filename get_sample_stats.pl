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
my @headers = qw(sample_name total_read trimmed_read ref_mapped ref_mapped_dedup vec_mapped ref_mapped_dedup_novec
peak_count peak_dedup_count target_count target_dedup_count insert_site_count clone_count clone_loc_count clone_loc_freq);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new MiSeqITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');
my $DEFAULT_MIN_CLONE_EXP = 2;

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
		my $in = $design->get_sample_ref_filtered_sorted_file($sample);
		$ref_mapped = `samtools view $WORK_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_mapped;
	}

# get ref dedup
  my $ref_dedup;
	{
		my $in = $design->get_sample_ref_dedup_file($sample);
		$ref_dedup = `samtools view -F 0x400 $WORK_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_dedup;
	}

# get vec mapped
	my $vec_mapped;
	{
		my $in = $design->get_sample_vec_filtered_file($sample);
		$vec_mapped = `samtools view $WORK_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $vec_mapped;
	}

# get ref_novec
  my $ref_novec;
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		$ref_novec = `samtools view -F 0x400 $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_novec;
	}

# get peak info
  my ($peak_count, $peak_dedup_count) = (0, 0);
	{
		my $in = $design->get_sample_ref_filtered_peak($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			chomp $line;
			$peak_count++;
			my ($rnames) = (split(/\t/, $line))[3];
			my %dedup_count;
			foreach my $rname (split(/,/, $rnames)) {
				$rname =~ s/\/\d+$//; # remove trailing /1 or /2
				$dedup_count{$rname}++;
			}
			my $count = scalar keys %dedup_count;
			$peak_dedup_count += $count;
		}
		close(BED);
	}

# get target and clone info
	my ($target_count, $target_dedup_count) = (0, 0);
	my ($site_count, $clone_count, $clone_loc_count) = (0, 0, 0);
	my %clone_loc_freq;
	my $target_file = $design->sample_opt($sample, 'target_file');
	if(-e $target_file) { # a gene editing sample
		($site_count, $clone_count, $clone_loc_count) = qw(NA NA NA);
		my $in = $design->get_sample_ref_filtered_peak($sample);
		if(-s "$BASE_DIR/$in") { # non-empty peaks found
			open(BED, "$bedtools intersect -a $BASE_DIR/$in -b $target_file -wo |") || die "Unable to open $samtools intersect: $!";
			while(my $line = <BED>) {
				chomp $line;
				$target_count++;
				my ($rnames) = (split(/\t/, $line))[3];
				my %dedup_count;
				foreach my $rname (split(/,/, $rnames)) {
					$rname =~ s/\/\d+//;
					$dedup_count{$rname}++;
				}
				my $count = scalar keys %dedup_count;
				$target_dedup_count += $count;
			}
			close(BED);
		}
	}
	else { # a gene therapy sample
		($target_count, $target_dedup_count) = qw(NA NA);
		my $in = $design->get_sample_ref_merged_clone($sample);
		$site_count = `wc -l < $WORK_DIR/$in`; chomp $site_count;

		$in = $design->get_sample_ref_filtered_clone($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			next if($line =~ /^#/);
			chomp $line;
			my ($attrs) = (split(/\t/, $line))[3];
			my ($loc_count) = $attrs =~ /LocCount=(\d+)/;
			$clone_count++;
			$clone_loc_count += $loc_count;
			$clone_loc_freq{$loc_count}++;
		}
		close(BED);
	}

# output
  print OUT "$sample\t$total_read\t$trimmed_read\t$ref_mapped\t$ref_dedup\t$vec_mapped\t$ref_novec\t$peak_count\t$peak_dedup_count\t",
  "$target_count\t$target_dedup_count\t",
  "$site_count\t$clone_count\t$clone_loc_count\t", get_freq_str(%clone_loc_freq), "\n";
}

close(OUT);

# subroutine definitions
sub get_freq_str {
  return 'NA' if(!@_);
  my %freq = @_;
  return join(',', map { "$_:$freq{$_}" } sort {$a <=> $b} keys %freq);
}
