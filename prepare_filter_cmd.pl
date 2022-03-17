#!/bin/env perl
# Prepare sh script for filtering reference mapping files
our $VERSION = v1.1;
our $ENV_FILE = 'set_filter_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/ITR_Seq';
use MiSeqITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $picard = 'picard.jar';
my $peak_script = 'get_peak_from_merged.pl';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new MiSeqITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
#my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');
my $KEEP_UNPAIR = $design->get_global_opt('KEEP_UNPAIR');
my $KEEP_STRAND = $design->get_global_opt('KEEP_STRAND');

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
	print STDERR "Error: BASE_DIR $BASE_DIR not exists\n";
	exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
	my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');
# prepare vec filter cmd
	{
		my $in = $design->get_sample_vec_map_file($sample);
		my $out = $design->get_sample_vec_filtered_file($sample);
		my $cmd = "$samtools view -F 0x4 -q $min_mapQ $WORK_DIR/$in -b -o $WORK_DIR/$out"; # 0x4 => unmap

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: filtered vector file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare vec sorted cmd
	{
		my $in = $design->get_sample_vec_filtered_file($sample);
		my $out = $design->get_sample_vec_sorted_file($sample);
		my $cmd = "$samtools sort $WORK_DIR/$in -o $BASE_DIR/$out"; # 0x4 => unmap
# index the bam file
   $cmd .= "\n$samtools index $BASE_DIR/$out";

		if(!-e "$BASE_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: filtered sorted vec file exists, won't override\n";
			$cmd =~ s/\n/\n# /sg; # add # after intermediate new-lines
			print OUT "# $cmd\n";
		}
	}

# prepare vec ID cmd
	{
		my $in = $design->get_sample_vec_map_file($sample);
		my $out = $design->get_sample_vec_ID_file($sample);
		my $cmd = "$samtools view -F 0x4 $WORK_DIR/$in | cut -f1 > $WORK_DIR/$out";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: filtered sorted vec ID file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare ref filter cmd
	{
		my $in = $design->get_sample_ref_map_file($sample);
		my $out = $design->get_sample_ref_filtered_file($sample);
		my $cmd = "$samtools view -F 0x4 -q $min_mapQ $WORK_DIR/$in -b -o $WORK_DIR/$out"; # 0x4 => unmap

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: ref filtered file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare ref novec cmd
	{
		my $id = $design->get_sample_vec_ID_file($sample);
		my $in = $design->get_sample_ref_filtered_file($sample);
		my $out = $design->get_sample_ref_novec_file($sample);
		my $cmd = "java -jar $SCRIPT_DIR/$picard FilterSamReads -I $WORK_DIR/$in -O $WORK_DIR/$out --FILTER excludeReadList -RLF $WORK_DIR/$id";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: ref novec file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare ref sort cmd
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		my $out = $design->get_sample_ref_sorted_file($sample);
		my $cmd = "$samtools sort $WORK_DIR/$in -o $BASE_DIR/$out"; # 0x4 => unmap
# index the bam file
   $cmd .= "\n$samtools index $BASE_DIR/$out";

		if(!-e "$BASE_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: filtered sorted map file exists, won't override\n";
			$cmd =~ s/\n/\n# /sg; # add # after intermediate new-lines
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
