#!/bin/env perl
# Prepare sh script for mapping ITR-trimmed reads to given reference genome database, using chosen aligner
our $VERSION = v1.1;
our $ENV_FILE = 'set_annotate_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/ITR_Seq';
use MiSeqITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $annotate_prog = 'bedtools';
my $track_script = 'get_BED_track.pl';
my $classify_script = 'get_BEDdetail_gtype_summ.pl';
my $extract_script = 'extract_BED_seq.pl';
#my $summ_script = 'show_gtype_summ.R';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new MiSeqITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
#my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
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

if(!(-e $WORK_DIR)) {
	print STDERR "Error: WORK_DIR $WORK_DIR not exists, creating now\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare anno cmd
	{
		my $ref_in = $design->get_sample_ref_filtered_peak($sample);
		my $vec_in = $design->get_sample_vec_sorted_file($sample);
		my $ref_out = $design->get_sample_ref_peak_anno($sample);
		my $vec_out = $design->get_sample_vec_map_anno($sample);
		my $ref_gff = $design->sample_opt($sample, 'ref_gff');
		my $vec_gff = $design->get_sample_vec_anno($sample);

		my $ref_cmd = "$annotate_prog intersect -a $BASE_DIR/$ref_in -b $ref_gff -wao > $WORK_DIR/$ref_out";
		my $vec_cmd = "$annotate_prog intersect -a $BASE_DIR/$vec_in -b $WORK_DIR/$vec_gff -wao -bed > $WORK_DIR/$vec_out";

		if(!(-e "$WORK_DIR/$ref_out")) {
			print OUT "$ref_cmd\n";
		}
		else {
			print STDERR "Warning: annotated ref peak file already exists, won't override\n";
			print OUT "# $ref_cmd\n";
		}

		if(!(-e "$WORK_DIR/$vec_out")) {
			print OUT "$vec_cmd\n";
		}
		else {
			print STDERR "Warning: annotated vec map file already exists, won't override\n";
			print OUT "# $vec_cmd\n";
		}
	}

# prepare track cmd
  {
		my $ref_in = $design->get_sample_ref_filtered_peak($sample);
		my $ref_out = $design->get_sample_ref_peak_track($sample);

		my $ref_cmd = "$SCRIPT_DIR/$track_script $BASE_DIR/$ref_in $BASE_DIR/$ref_out $sample-ITR-peak";

		if(!(-e "$BASE_DIR/$ref_out")) {
			print OUT "$ref_cmd\n";
		}
		else {
			print STDERR "Warning: annotated ref map file already exists, won't override\n";
			print OUT "# $ref_cmd\n";
		}
	}

# prepare classify anno cmd
	{
		my $ref_in = $design->get_sample_ref_peak_anno($sample);
		my $vec_in = $design->get_sample_vec_map_anno($sample);
		my $ref_out = $design->get_sample_ref_peak_gtype($sample);
		my $vec_out = $design->get_sample_vec_map_gtype($sample);
    
		my $ref_cmd = "$SCRIPT_DIR/$classify_script $WORK_DIR/$ref_in $BASE_DIR/$ref_out";
		my $vec_cmd = "$SCRIPT_DIR/$classify_script $WORK_DIR/$vec_in $BASE_DIR/$vec_out";

		if(!(-e "$BASE_DIR/$ref_out")) {
			print OUT "$ref_cmd\n";
		}
		else {
			print STDERR "Warning: classificaiton summary for ref peak file already exists, won't override\n";
			print OUT "# $ref_cmd\n";
		}
		if(!(-e "$BASE_DIR/$vec_out")) {
			print OUT "$vec_cmd\n";
		}
		else {
			print STDERR "Warning: classification summary for vec map file already exists, won't override\n";
			print OUT "# $vec_cmd\n";
		}
	}

# prepare extract seq cmd
	{
		my $db = $design->sample_opt($sample, 'genome_seq');
		my $in = $design->get_sample_ref_filtered_peak($sample);
		my $out = $design->get_sample_ref_peak_seq($sample);
    
		my $cmd = "$SCRIPT_DIR/$extract_script $db $BASE_DIR/$in $BASE_DIR/$out";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: ref map seq file already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

=pod
# prepare show gtype summ cmd
	{
		my $ref_in = $design->get_sample_map_gtype_file($sample);
		my $vec_in = $design->get_sample_remap_gtype_file($sample);
		my $ref_out = $design->get_sample_map_gtype_figure($sample);
		my $vec_out = $design->get_sample_remap_gtype_figure($sample);
    
		my $ref_cmd = "$summ_script $BASE_DIR/$ref_in $BASE_DIR/$ref_out";
		my $vec_cmd = "$summ_script $BASE_DIR/$vec_in $BASE_DIR/$vec_out";

		if(!(-e "$BASE_DIR/$ref_out")) {
			print OUT "$ref_cmd\n";
		}
		else {
			print STDERR "Warning: gtype summary figure for ref map file already exists, won't override\n";
			print OUT "# $ref_cmd\n";
		}
		if(!(-e "$BASE_DIR/$vec_out")) {
			print OUT "$vec_cmd\n";
		}
		else {
			print STDERR "Warning: summary figure for vector remap file already exists, won't override\n";
			print OUT "# $vec_cmd\n";
		}
	}
=cut

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
