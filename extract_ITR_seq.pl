#!/bin/env perl
# This script used to extract the 3' ITR sequence from the VECTOR GenBank file and
# use it (and its reverse-complement) as trimming adapters for the R1 and R2 reads
use strict;
use warnings;

use Bio::SeqIO;

my $usage = "Usage: $0 GenBank-INFILE ITR5-OUTFILE ITR3-OUTFILE";
my $infile = shift or die $usage;
my $out3file = shift or die $usage;
my $out5file = shift or die $usage;

# open GenBank input
my $in = new Bio::SeqIO(-file => "<$infile", -format => 'genbank');

# open FASTA output
my $out5 = new Bio::SeqIO(-file => ">$out5file", -format => "fasta", -alphabet => 'dna');
my $out3 = new Bio::SeqIO(-file => ">$out3file", -format => "fasta", -alphabet => 'dna');

# read each GenBank seq
while(my $seqobj = $in->next_seq()) {
	my $id = $seqobj->display_id();
# search each FEATURE entity
	my $ITR5_feat; # as the 5' most ITR feature
	my $ITR3_feat; # as the 3' most ITR feature
	foreach my $feat ($seqobj->get_SeqFeatures()) {
		if($feat->has_tag('label')) {
#$feat->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
			foreach my $val ($feat->get_tag_values('label')) {
				if($val =~ /ITR/) { # is an ITR seq feature
					if(!defined($ITR5_feat) || $feat->start() < $ITR5_feat->start()) { # a more 5' ITR found
						$ITR5_feat = $feat;
					}
					if(!defined($ITR3_feat) || $feat->start() > $ITR3_feat->start()) { # a more 3' ITR found
						$ITR3_feat = $feat;
					}
				}
			}
		}
	}
	my $ITR5_seqobj = $ITR5_feat->seq();
	my $ITR3_seqobj = $ITR3_feat->seq();
# update ITR seq id and desc
	my $ITR5_start = $ITR5_feat->start();
	my $ITR5_end = $ITR5_feat->end();
	my $ITR5_len = $ITR5_feat->length();
	$ITR5_seqobj->display_id("$id:ITR5");
	$ITR5_seqobj->desc("$ITR5_start-$ITR5_end:$ITR5_len");

	my $ITR3_start = $ITR3_feat->start();
	my $ITR3_end = $ITR3_feat->end();
	my $ITR3_len = $ITR3_feat->length();
	$ITR3_seqobj->display_id("$id:ITR3");
	$ITR3_seqobj->desc("$ITR3_start-$ITR3_end:$ITR3_len");

# output
	$out5->write_seq($ITR5_seqobj);
	$out3->write_seq($ITR3_seqobj->revcom()); # use the revcom seq as the adapter sequence
}

$in->close();
$out5->close();
$out3->close();
