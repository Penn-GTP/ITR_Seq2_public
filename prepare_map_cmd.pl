#!/bin/env perl
# Prepare bash script for remapping trimmed reads to reference genome database
our $VERSION = 'v2.1.1';
our $ENV_FILE = 'set_map_env.sh';

use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $vec_anno_script = 'get_vector_seq_anno.pl';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $cmd = "$0 " . join(" ", @ARGV);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new ITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $VECTOR_DIR = $design->get_global_opt('VECTOR_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
#my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
	print STDERR "Error: BASE_DIR $BASE_DIR not exists\n";
	exit;
}

if(!(-e $VECTOR_DIR && -d $VECTOR_DIR)) {
	print STDERR "Error: BASE_DIR $VECTOR_DIR not exists\n";
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
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare ref map cmd
	my $aligner = $design->sample_opt($sample, 'aligner');
	my $ref_db = $design->sample_opt($sample, 'ref_db');
	my $align_opts = $design->sample_opt($sample, 'align_opts');
	if(!($aligner eq 'bowtie2')) {
		print STDERR "Error: trim_prog now only support 'bowtie2'\n";
		exit;
	}

  my $in1 = $design->get_sample_fwd_ITRtrim_file($sample);
  my $in2 = $design->get_sample_rev_ITRtrim_file($sample);

  my $out = $design->get_sample_ref_map_file($sample);

  my $cmd;
	if($aligner eq 'bowtie2') {
		$cmd = "bowtie2 -x $ref_db -1 $WORK_DIR/$in1 -2 $WORK_DIR/$in2 -p $NUM_PROC --local $align_opts | $samtools view -b -o $WORK_DIR/$out -";
	}
	else {
		print STDERR "Error: Unsupported aligner '$aligner'\n";
		exit;
	}

  if(!(-e "$WORK_DIR/$out")) {
		print OUT "$cmd\n";
  }
  else {
    print STDERR "Warning: ref alignment file already exists, won't override\n";
		print OUT "# $cmd\n";
  }

# prepare vec map cmd
  {
    my $aligner = $design->sample_opt($sample, 'aligner');
    my $align_opts = $design->sample_opt($sample, 'align_opts');
    my $in = $design->get_sample_vec_seq($sample);
    my $dbname = $design->get_sample_vec_dbname($sample);
    my $in1 = $design->get_sample_fwd_ITRtrim_file($sample);
    my $in2 = $design->get_sample_rev_ITRtrim_file($sample);
    my $out = $design->get_sample_vec_map_file($sample);
    my $cmd;
    if($aligner eq 'bowtie2') {
      $cmd .= "bowtie2 -x $VECTOR_DIR/$dbname -U $WORK_DIR/$in1 -U $WORK_DIR/$in2 -p $NUM_PROC --local $align_opts | $samtools view -b -o $WORK_DIR/$out -";
    }
    else {
      print STDERR "Error: Unsupported aligner '$aligner'\n";
      exit;
    }

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
      print STDERR "Warning: vector map file exists, won't override\n";
      $cmd =~ s/\n/\n# /sg; # comment out to intermediate lines
      print OUT "# $cmd\n";
    }
  }

  print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
