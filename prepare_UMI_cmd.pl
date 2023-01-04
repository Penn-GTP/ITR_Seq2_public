#!/bin/env perl
# Prepare sh script for identify and label UMIs for a MiSeq-ITR-Seq run using customized conf file
our $VERSION = v1.1;
our $ENV_FILE = 'set_UMI_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/ITR_Seq2';
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $script_name = 'get_UMI_fastq_seqs.pl';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new ITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
my $UMI_LEN = $design->get_global_opt('UMI_LEN');

# check required directories
if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $DEMUX_DIR)) {
	print STDERR "Error: DEMUX $DEMUX_DIR not exists\n";
	exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Warning: WORK_DIR $WORK_DIR not exists, creating now\n";
	mkdir($WORK_DIR, 0750) || die "Unable to mkdir $WORK_DIR: $!\n";
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
	my $in1 = $design->sample_opt($sample, 'fastq_R1');
	my $in2 = $design->sample_opt($sample, 'fastq_R2');
	my $idx = $design->sample_opt($sample, 'fastq_I2');

  my $out1 = $design->get_sample_fwd_UMI_file($sample);
  my $out2 = $design->get_sample_rev_UMI_file($sample);

  my $cmd = "$SCRIPT_DIR/$script_name -i1 $DEMUX_DIR/$in1 -i2 $DEMUX_DIR/$in2 -idx $DEMUX_DIR/$idx -o1 $WORK_DIR/$out1 -o2 $WORK_DIR/$out2 --len $UMI_LEN";
  if(!(-e "$WORK_DIR/$out1" && -e "$WORK_DIR/$out2")) {
    print OUT "$cmd\n";
  }
  else {
    print STDERR "Warning: UMI labeled FASTQ files already exists, won't override\n";
    print OUT "# $cmd\n";
  }
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
