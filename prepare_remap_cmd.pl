#!/bin/env perl
# Prepare sh script for remapping failed (filtered out) reads to the vector genome
our $VERSION = v1.1;
our $ENV_FILE = 'set_remap_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/ITR_Seq';
use MiSeqITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $vector_anno_script = 'get_vector_seq_anno.pl';
my $samtools = 'samtools';
my $bedtools = 'bedtools';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new MiSeqITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
#my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
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
	print STDERR "Error: WORK_DIR $VECTOR_DIR not exists\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare get vector seq and annotation cmd
	{
		my $in = $design->sample_opt($sample, 'vector_file');
		my $seq_out = $design->get_sample_vector_seqfile($sample);
		my $anno_out = $design->get_sample_vector_annofile($sample);
		my $cmd = "$SCRIPT_DIR/$vector_anno_script $VECTOR_DIR/$in $WORK_DIR/$seq_out $WORK_DIR/$anno_out";
		if(!(-e "$WORK_DIR/$seq_out" && -e "$WORK_DIR/$anno_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: vector seq and annotation file exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare build vector db and map cmd
  {
		my $aligner = $design->sample_opt($sample, 'aligner');
		my $align_opts = $design->sample_opt($sample, 'align_opts');
		my $in = $design->get_sample_vector_seqfile($sample);
    my $dbname = $design->get_sample_vector_dbname($sample);
		my $in1 = $design->get_sample_map_failed_fwd_file($sample);
		my $in2 = $design->get_sample_map_failed_rev_file($sample);
		my $out = $design->get_sample_remap_file($sample);
    my $cmd;
		if($aligner eq 'bowtie2') {
			$cmd .= "bowtie2-build $WORK_DIR/$in $WORK_DIR/$dbname\n"; # add build cmd
      $cmd .= "bowtie2 -x $WORK_DIR/$dbname -1 $WORK_DIR/$in1 -2 $WORK_DIR/$in2 $align_opts -p $NUM_PROC | samtools view -b -o $WORK_DIR/$out -";
		}
		else {
			print STDERR "Error: Unsupported aligner '$aligner'\n";
			exit;
		}
		
		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: vector remap file exists, won't override\n";
			$cmd =~ s/\n/\n# /sg; # comment out to intermediate lines
			print OUT "# $cmd\n";
		}
	}

# prepare filter & sort remap result
	{
		my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');
		my $in = $design->get_sample_remap_file($sample);
		my $out = $design->get_sample_remap_filtered_file($sample);
#   add filter & sort cmd
    my $cmd = "$samtools view -F 0x4 -q $min_mapQ $WORK_DIR/$in -b | $samtools sort - -o $BASE_DIR/$out\n";
# add index cmd
    $cmd .= "$samtools index $BASE_DIR/$out";
    if(!-e "$BASE_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
      print STDERR "Warning: filtered remap file exists, won't override\n";
      $cmd =~ s/\n/\n# /sg; # add # after intermediate new-lines
      print OUT "# $cmd\n";
    }
	}

# prepare filter remap result
	{
		my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');
		my $in = $design->get_sample_remap_file($sample);
		my $out = $design->get_sample_remap_insert_file($sample);
#   add alignment to output
    my $cmd = "$samtools view -F 0x4 -q $min_mapQ -b $WORK_DIR/$in | bedtools bamtobed -i - > $WORK_DIR/$out";
    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
      print STDERR "Warning: filtered remap file exists, won't override\n";
      $cmd =~ s/\n/\n# /sg; # add # after intermediate new-lines
      print OUT "# $cmd\n";
    }
	}

# prepare sort cmd
  {
    my $in = $design->get_sample_remap_insert_file($sample);
    my $out = $design->get_sample_remap_sorted_file($sample);
    my $cmd = "$bedtools sort -i $WORK_DIR/$in > $WORK_DIR/$out";
    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
      print STDERR "Warning: sorted remap file exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
