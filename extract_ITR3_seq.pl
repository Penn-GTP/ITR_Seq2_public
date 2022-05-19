#!/bin/env perl
# This script used to extract the fwd and revcom sequences of the 3' ITR from the VECTOR GenBank file and
# use it (and its reverse-complement) as trimming adapters for the R1 and R2 reads
use strict;
use warnings;

use Bio::SeqIO;

my $usage = "Usage: $0 GenBank-INFILE ITR-FWD-OUTFILE ITR-REV-OUTFILE";
my $infile = shift or die $usage;
my $outFwdfile = shift or die $usage;
my $outRevfile = shift or die $usage;

# open GenBank input
my $in = new Bio::SeqIO(-file => "<$infile", -format => 'genbank');

# open FASTA output
my $outFwd = new Bio::SeqIO(-file => ">$outFwdfile", -format => "fasta", -alphabet => 'dna');
my $outRev = new Bio::SeqIO(-file => ">$outRevfile", -format => "fasta", -alphabet => 'dna');

# read each GenBank seq
while(my $seqobj = $in->next_seq()) {
	my $id = $seqobj->display_id();
# search each FEATURE entity
  my $ITR3_feat; # as the 3' most ITR feature
	foreach my $feat ($seqobj->get_SeqFeatures()) {
		if($feat->has_tag('label')) {
#$feat->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
			foreach my $val ($feat->get_tag_values('label')) {
				if($val =~ /ITR/) { # is an ITR seq feature
					if(!defined($ITR3_feat) || $feat->start() > $ITR3_feat->start()) { # a more 3' ITR found
						$ITR3_feat = $feat;
					}
				}
			}
		}
	}
	my $ITRfwd_seqobj = $ITR3_feat->seq(); # original seq of the ITR'3 seq
	my $ITRrev_seqobj = $ITR3_feat->seq()->revcom(); # revcom seq of the ITR'3 seq
# update ITR seq id and desc
	my $ITR3_start = $ITR3_feat->start();
	my $ITR3_end = $ITR3_feat->end();
	my $ITR3_len = $ITR3_feat->length();
	$ITRfwd_seqobj->display_id("$id:ITR3");
	$ITRfwd_seqobj->desc("$ITR3_start-$ITR3_end:$ITR3_len");

	$ITRrev_seqobj->display_id("$id:ITR3");
	$ITRrev_seqobj->desc("$ITR3_start-$ITR3_end:$ITR3_len");

# output
	$outFwd->write_seq($ITRfwd_seqobj);
	$outRev->write_seq($ITRrev_seqobj); # use the revcom seq as the adapter sequence
}

$in->close();
$outFwd->close();
$outRev->close();
