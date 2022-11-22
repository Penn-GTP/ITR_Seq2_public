#!/bin/env perl
# Prepare sh script for getting ITR peaks from filtered alignments
our $VERSION = v1.1;
our $ENV_FILE = 'set_peak_env.sh';

use strict;
use warnings;
use Bio::SeqIO;
use lib '/project/gtplab/pipeline/ITR_Seq2';
use MiSeqITRSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
#my $picard = 'picard.jar';
my $insert_script = 'get_ITR_insert_site.pl';
my $extract_script = 'extract_insert_site_flank_seq.pl';
my $aligner = 'water';
my $filter_script = 'filter_insert_site.pl';
my $peak_script = 'get_ITR_peak.pl';
my $clone_script = 'get_ITR_clone.pl';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new MiSeqITRSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
#my $DEMUX_DIR = $design->get_global_opt('DEMUX_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
#my $UMI_LEN = $design->get_global_opt('UMI_LEN');
my $INSERT_SIZE = $design->get_global_opt('INSERT_SIZE');
my $KEEP_UNPAIR = $design->get_global_opt('KEEP_UNPAIR');
my $KEEP_STRAND = $design->get_global_opt('KEEP_STRAND');
my $MAX_PEAK_DIST = $design->get_global_opt('MAX_PEAK_DIST');

my $DEFAULT_MIN_CLONE_LOC = 2;

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
		my $out = $design->get_sample_ref_insert_site($sample); 
		my $primer_file = $design->get_global_primer_fwd();
		my $clip_len = 0;
		my $seq_in = new Bio::SeqIO(-file => "<$BASE_DIR/$primer_file", -format => 'fasta');
		while(my $seq_obj = $seq_in->next_seq()) {
			if($clip_len == 0 || $seq_obj->length() < $clip_len) {
				$clip_len = $seq_obj->length();
			}
		}

		my $cmd = "$SCRIPT_DIR/$insert_script $BASE_DIR/$in $WORK_DIR/$out --insert-size $INSERT_SIZE --min-softclip $clip_len";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site sorted cmd
  {
		my $in = $design->get_sample_ref_insert_site($sample);
		my $out = $design->get_sample_ref_insert_site_sorted($sample);
		my $cmd = "$bedtools sort -i $WORK_DIR/$in > $WORK_DIR/$out";
		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site uniq cmd
  {
		my $in = $design->get_sample_ref_insert_site_sorted($sample);
		my $out = $design->get_sample_ref_insert_site_uniq($sample);

		my $cmd = "if [ -s $WORK_DIR/$in ]; then $bedtools merge -d -$INSERT_SIZE -c 4,5,6 -o collapse,sum,collapse -i $WORK_DIR/$in > $WORK_DIR/$out ; else cp $WORK_DIR/$in $WORK_DIR/$out ; fi;";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare extract insert site flank cmd
	{
		my $db = $design->sample_opt($sample, 'genome_seq');
		my $in = $design->get_sample_ref_insert_site_uniq($sample);
		my $seq_out = $design->get_sample_ref_insert_site_flank_seq($sample);
		my $name_out = $design->get_sample_ref_insert_site_name2loc($sample);
		my $ext_len = $design->get_global_opt('PRIMER_FLANK');
		my $cmd = "$SCRIPT_DIR/$extract_script $db $WORK_DIR/$in $WORK_DIR/$seq_out $WORK_DIR/$name_out --ext-len $ext_len";

		if(!(-e "$WORK_DIR/$seq_out" && -e "$WORK_DIR/$name_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$seq_out and $WORK_DIR/$name_out exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site flank primer align cmd
	{
		my $primer_fwd = $design->get_global_primer_fwd();
		my $primer_rev = $design->get_global_primer_rev();
		my $opts = $design->get_global_opt('PRIMER_ALN_OPTS');

		my $in = $design->get_sample_ref_insert_site_flank_seq($sample);
		my $out_fwd = $design->get_sample_ref_insert_site_flank_fwd_align($sample);
		my $out_rev = $design->get_sample_ref_insert_site_flank_rev_align($sample);
		my $cmd = "> $WORK_DIR/$out_fwd \n> $WORK_DIR/$out_rev\n";

		my $in_fwd = new Bio::SeqIO(-file => "<$primer_fwd", -format => 'fasta', -alphabet => 'dna');
		while(my $seq_fwd = $in_fwd->next_seq()) {
			my $id = $seq_fwd->id();
			my $seq = $seq_fwd->seq();
			$cmd .= "if [ -s $WORK_DIR/$in ]; then echo '$seq' | $aligner -asequence stdin -bsequence $WORK_DIR/$in -auto -sid1 $id -sformat1 plain -stdout $opts >> $WORK_DIR/$out_fwd;";
			$cmd .= "\nelse >> $WORK_DIR/$out_fwd; fi;\n";
		}
	
		my $in_rev = new Bio::SeqIO(-file => "<$primer_rev", -format => 'fasta', -alphabet => 'dna');
		while(my $seq_rev = $in_rev->next_seq()) {
			my $id = $seq_rev->id();
			my $seq = $seq_rev->seq();
			$cmd .= "if [ -s $WORK_DIR/$in ]; then echo '$seq' | $aligner -asequence stdin -bsequence $WORK_DIR/$in -auto -sid1 $id -sformat1 plain -stdout $opts >> $WORK_DIR/$out_rev;";
			$cmd .= "\nelse >> $WORK_DIR/$out_rev; fi;\n";
		}
	
		if(!(-e "$WORK_DIR/$out_fwd" && -e "$WORK_DIR/$out_rev")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out_fwd and $WORK_DIR/$out_rev exist, won't override\n";
			$cmd =~ s/\n/\n# /sg;
			print OUT "# $cmd\n";
		}
	}

# prepare insert site filtered cmd
	{
		my $site_in = $design->get_sample_ref_insert_site_uniq($sample);
		my $name_in = $design->get_sample_ref_insert_site_name2loc($sample);
		my $fwd_aln = $design->get_sample_ref_insert_site_flank_fwd_align($sample);
		my $rev_aln = $design->get_sample_ref_insert_site_flank_rev_align($sample);
		my $site_out = $design->get_sample_ref_insert_site_filtered($sample);
		my $info_out = $design->get_sample_ref_insert_site_align_info($sample);

		my $flank_size = $design->get_global_opt('PRIMER_FLANK');
		my $seed_len = $design->get_global_opt('PRIMER_SEED_LEN');
		my $primer_file = $design->get_global_primer_fwd();
		my $max_seed_error = $design->get_global_opt('PRIMER_MAX_SEED_ERROR');
		my $min_match = $design->get_global_opt('PRIMER_MIN_MATCH');
		my $cmd = "$SCRIPT_DIR/$filter_script $WORK_DIR/$site_in $WORK_DIR/$name_in $WORK_DIR/$fwd_aln $WORK_DIR/$rev_aln "
		. "$WORK_DIR/$site_out $WORK_DIR/$info_out "
		. "--flank-size $flank_size --seed-len $seed_len --primer-file $BASE_DIR/$primer_file --max-seed-error $max_seed_error --min-match $min_match";

		if(!(-e "$WORK_DIR/$site_out" && -e "$WORK_DIR/$info_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$site_out and $WORK_DIR/$info_out exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare insert site merged cmd
  {
		my $in = $design->get_sample_ref_insert_site_filtered($sample);
		my $out = $design->get_sample_ref_insert_site_merged($sample);

		my $cmd = "if [ -s $WORK_DIR/$in ]; then $bedtools merge -d $MAX_PEAK_DIST -c 4,5,6 -o collapse,sum,collapse -i $WORK_DIR/$in > $WORK_DIR/$out ; else cp $WORK_DIR/$in $WORK_DIR/$out ; fi;";

		if(!-e "$WORK_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare call peak cmd
	{
		my $in = $design->get_sample_ref_insert_site_merged($sample);
		my $out = $design->get_sample_ref_peak($sample);
		my $cmd = "$SCRIPT_DIR/$peak_script $WORK_DIR/$in $BASE_DIR/$out --keep-strand $KEEP_STRAND";

		if(!-e "$BASE_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare call clone cmd
  {
		my $site_in = $design->get_sample_ref_insert_site_filtered($sample);
		my $aln_in = $design->get_sample_ref_novec_file($sample);
		my $out = $design->get_sample_ref_clone($sample);
		my $min_clone_loc = $design->sample_opt($sample, 'min_clone_loc') ? $design->sample_opt($sample, 'min_clone_loc') : $DEFAULT_MIN_CLONE_LOC;
		my $cmd = "$SCRIPT_DIR/$clone_script $WORK_DIR/$site_in $BASE_DIR/$aln_in $BASE_DIR/$out --min-loc $min_clone_loc";

		if(!-e "$BASE_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;

sub revcom {
  my $seq = shift;
  $seq = reverse $seq;
  $seq =~ tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
  return $seq;
}

