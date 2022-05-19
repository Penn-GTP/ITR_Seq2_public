#!/bin/env perl
# Prepare sh script for getting ITR peaks from filtered alignments
our $VERSION = v1.1;
our $ENV_FILE = 'set_peak_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/ITR_Seq';
use MiSeqITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
#my $picard = 'picard.jar';
my $insert_script = 'get_ITR_insertion_site.pl';
my $peak_script = 'get_ITR_peak.pl';

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
# prepare get ITR insert site cmd
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		my $out = $design->get_sample_ref_peak($sample); 
		my $INSERT_SIZE = $design->get_global_opt('INSERT_SIZE');
		my $clip_len = length($design->get_global_opt('ITR_PRIMER'));

		my $cmd = "$SCRIPT_DIR/$insert_script $BASE_DIR/$in $WORK_DIR/$out --insert-size $INSERT_SIZE --min-softclip $clip_len";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: ref peak file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare sort peak cmd
  {
		my $in = $design->get_sample_ref_peak($sample);
		my $out = $design->get_sample_ref_sorted_peak($sample);
		my $cmd = "$bedtools sort -i $WORK_DIR/$in > $WORK_DIR/$out";
		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: sorted map file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare merge peak cmd
  {
		my $in = $design->get_sample_ref_sorted_peak($sample);
		my $out = $design->get_sample_ref_merged_peak($sample);
		my $max_dist = $design->sample_opt($sample, 'max_dist');

		my $cmd = "if [ -s $WORK_DIR/$in ]; then $bedtools merge -d $max_dist -c 4,5,6 -o collapse,sum,collapse -i $WORK_DIR/$in > $WORK_DIR/$out ; else cp $WORK_DIR/$in $WORK_DIR/$out ; fi;";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: merged file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare filter peak cmd
	{
		my $in = $design->get_sample_ref_merged_peak($sample);
		my $out = $design->get_sample_ref_filtered_peak($sample);
		my $cmd = "$SCRIPT_DIR/$peak_script $WORK_DIR/$in $BASE_DIR/$out $KEEP_STRAND";

		if(!-e "$BASE_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: peak file exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
