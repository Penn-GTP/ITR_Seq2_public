#!/bin/env perl
# Prepare sh script for filtering reference mapping files
our $VERSION = 'v2.1.2';

use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
#my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $cmd = "$0 " . join(" ", @ARGV);
my @comments = qq(Sample name\tTotal reads\tITR-containing reads\tHost-mapped reads\tHost-mapped deduplexed reads\tVector-mapped reads\tHost-mapped deduplexed non-vector reads\tUnique insert sites\tUnique insert sites mispriming filtered\tMerged insert peaks\tRead abundance of insert peaks\tOn-target insert peaks\tRead abundance of on-target insert peaks\tClonal insert sites\tUMI-locus abundance of clonal insert sites\tFrequency of UMI locus abundance of clonal insert sites);
my @headers = qw(sample_name total_read trimmed_read ref_mapped ref_mapped_dedup vec_mapped ref_mapped_dedup_novec
insert_site_uniq insert_site_filtered
peak_count peak_clone target_count target_clone clonal_count clonal_loc_count clonal_loc_freq);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new ITRSeqExpDesign($infile);
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
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
print OUT "# ", join("\t", @comments), "\n";
print OUT join("\t", @headers), "\n";

foreach my $sample ($design->get_sample_names()) {
	print STDERR "getting stats for $sample\n";
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
		my $in = $design->get_sample_vec_sorted_file($sample);
		$vec_mapped = `samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $vec_mapped;
	}

# get ref_novec
  my $ref_novec;
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		$ref_novec = `samtools view -F 0x400 $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_novec;
	}

# get insert site info
	my ($site_uniq, $site_filtered) = (0, 0);
	{
		my $in = $design->get_sample_ref_insert_site_uniq($sample);
		$site_uniq = `cat $WORK_DIR/$in | wc -l`; chomp $site_uniq;

		$in = $design->get_sample_ref_insert_site_filtered($sample);
		$site_filtered = `cat $WORK_DIR/$in | wc -l`; chomp $site_filtered;
	}

# get peak info
  my ($peak_count, $peak_clone) = (0, 0);
	{
		my $in = $design->get_sample_ref_peak($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			next if($line =~ /^#/);
			chomp $line;
			$peak_count++;
			my $peak_name = (split(/\t/, $line))[3];
			my ($dedup_count) = $peak_name =~ /ReadCount=(\d+)/;
			$peak_clone += $dedup_count;
		}
		close(BED);
	}

# get target info
	my ($target_count, $target_clone) = (0, 0);
	my $target_file = $design->sample_opt($sample, 'target_file');
	if(-e $target_file) { # a gene editing sample
		my $in = $design->get_sample_ref_peak($sample);
		if(-s "$BASE_DIR/$in") { # non-empty peaks found
			open(BED, "$bedtools intersect -a $BASE_DIR/$in -b $target_file -wo |") || die "Unable to open $samtools intersect: $!";
			while(my $line = <BED>) {
				next if($line =~ /^#/);
				chomp $line;
				$target_count++;
				my ($peak_name) = (split(/\t/, $line))[3];
				my ($dedup_count) = $peak_name =~ /ReadCount=(\d+)/;
				$target_clone += $dedup_count;
			}
			close(BED);
		}
	}

# get clone info
	my ($clone_count, $clone_loc_count) = (0, 0, 0);
	my %clone_loc_freq;
	{
		my $in = $design->get_sample_ref_clone($sample);
		open(BED, "<$BASE_DIR/$in") || die "Unable to open $in: $!";
		while(my $line = <BED>) {
			next if($line =~ /^#/);
			chomp $line;
			my ($clone_name) = (split(/\t/, $line))[3];
			my ($loc_count) = $clone_name =~ /LocCount=(\d+)/;
			$clone_count++;
			$clone_loc_count += $loc_count;
			$clone_loc_freq{$loc_count}++;
		}
		close(BED);
	}

# output
  print OUT "$sample\t$total_read\t$trimmed_read\t$ref_mapped\t$ref_dedup\t$vec_mapped\t$ref_novec\t",
	"$site_uniq\t$site_filtered\t",
	"$peak_count\t$peak_clone\t$target_count\t$target_clone\t",
  "$clone_count\t$clone_loc_count\t", get_freq_str(%clone_loc_freq), "\n";
}

close(OUT);

# subroutine definitions
sub get_freq_str {
  return 'NA' if(!@_);
  my %freq = @_;
  return join(',', map { "$_:$freq{$_}" } sort {$a <=> $b} keys %freq);
}
