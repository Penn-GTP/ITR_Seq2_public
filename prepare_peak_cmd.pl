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
#my $picard = 'picard.jar';
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
# prepare get peak cmd
	{
		my $in = $design->get_sample_ref_novec_file($sample);
		my $out = $design->get_sample_ref_peak($sample); 
		my $pair_flag = '';
#   add paired alignments
		my $cmd = "$samtools view -f 0x2 -b $WORK_DIR/$in | bedtools bamtobed -i - > $WORK_DIR/$out";
# add unpaired alignment
		if($KEEP_UNPAIR & 0x1) { # keep forward 0x40 => forward
			$cmd .= "\n$samtools view -F 0x2 -f 0x40 -b $WORK_DIR/$in | bedtools bamtobed -i - >> $WORK_DIR/$out";
		}
		if($KEEP_UNPAIR & 0x2) { # keep reverse 0x80 => reverse
			$cmd .= "\n$samtools view -F 0x2 -f 0x80 -b $WORK_DIR/$in | bedtools bamtobed -i - >> $WORK_DIR/$out";
		}
		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: map peak file exists, won't override\n";
			$cmd =~ s/\n/\n# /sg; # add # after intermediate new-lines
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
