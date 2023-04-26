#!/bin/env perl
# Prepare bash script for trimming/identifying ITR adapter seq from R1 and R2 reads
our $VERSION = 'v1.1.1';
our $ENV_FILE = 'set_trim_env.sh';

use strict;
use warnings;
use Bio::SeqIO;
use lib '/project/gtplab/pipeline/ITR_Seq2';
use ITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $cmd = "$0 " . join(" ", @ARGV);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new ITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
my $VECTOR_DIR = $design->get_global_opt('VECTOR_DIR');
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
	print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
	exit;
}

if(!(-e $VECTOR_DIR)) {
	print STDERR "Error: VECTOR_DIR $VECTOR_DIR not exists\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
	my $trim_prog = $design->sample_opt($sample, 'trim_prog');
	my $max_error = $design->sample_opt($sample, 'max_error_rate');
	my $min_len = $design->sample_opt($sample, 'min_len');
	my $trim_opts = $design->sample_opt($sample, 'trim_opts');

	if(!($trim_prog eq 'cutadapt')) {
		print STDERR "Error: trim_prog now only support 'cutadapt'\n";
		exit;
	}

# prepare adapter trim cmd
  {
		my $primer_fwd = $design->get_global_primer_fwd();
		my $primer_rev = $design->get_global_primer_rev();

    # find min_overlap as the min length of primers
		my $min_overlap = 0;
		my $seq_in = new Bio::SeqIO(-file => "<$BASE_DIR/$primer_fwd", -format => 'fasta');
		while(my $seq_obj = $seq_in->next_seq()) {
			my $seq_len = $seq_obj->length();
			if($min_overlap == 0 || $seq_len < $min_overlap) {
				$min_overlap = $seq_len;
			}
		}

		my $in1 = $design->get_sample_fwd_UMI_file($sample);
		my $in2 = $design->get_sample_rev_UMI_file($sample);

		my $trout1 = $design->get_sample_fwd_ITRtrim_file($sample);
		my $trout2 = $design->get_sample_rev_ITRtrim_file($sample);

		my $unout1 = $design->get_sample_fwd_ITRuntrim_file($sample);
		my $unout2 = $design->get_sample_rev_ITRuntrim_file($sample);

		my $shout1 = $design->get_sample_fwd_ITRshort_file($sample);
		my $shout2 = $design->get_sample_rev_ITRshort_file($sample);

		my $cmd = "$trim_prog -a file:$primer_rev -G file:$primer_fwd -o $WORK_DIR/$trout1 -p $WORK_DIR/$trout2 " .
			"--untrimmed-output $WORK_DIR/$unout1 --untrimmed-paired-output $WORK_DIR/$unout2 " .
			"--too-short-output $WORK_DIR/$shout1 --too-short-paired-output $WORK_DIR/$shout2 " .
			"-O $min_overlap -m $min_len -e $max_error --cores $NUM_PROC --pair-filter both --action lowercase " .
			"$trim_opts $WORK_DIR/$in1 $WORK_DIR/$in2";

		if(!(-e "$WORK_DIR/$trout1" && -e "$WORK_DIR/$trout2")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: ITR-trimmed FASTQ files already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
