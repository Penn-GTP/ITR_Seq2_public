#!/usr/bin/env perl
# This script is ued to extract seq from DB using pos in BED file
use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

my $ext_len = 0;
my $options = qq([OPTIONS]
OPTIONS:
  --ext-len INT: extend INT length both up-stream and down-stream of the regions [$ext_len]
);

my $usage = "Usage: $0 DB-PATH INFILE OUTFILE $options";
my $db_path = shift or die $usage;
my $infile = shift or die $usage;
my $seq_out = shift or die $usage;
my $name_out = shift or die $usage;

GetOptions(
"ext-len=i" => \$ext_len)
or die "Error in command line arguments, usage: $usage";

my $db = new Bio::DB::Fasta($db_path);
open(IN, "<$infile") || die "Unable to open $infile: $!";
my $out = new Bio::SeqIO(-file => ">$seq_out", -format => 'fasta', -alphabet => 'dna');
open(NAME, ">$name_out") || die "Unable to write to $name_out: $!";

# read BED file
my $id = 0;
while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end) = split(/\t/, $line);
	my $loc = "$chr:$start-$end";
	if($ext_len > 0) {
		$start -= $ext_len;
		$end += $ext_len - 1;
		$start = 0 if($start < 0);
	}
	my $seq = $db->seq($chr, $start + 1 => $end);
	my $name = "loc" . (++$id);
	my $seq_obj = new Bio::Seq(-seq => $seq, -display_id => $name, -desc => $loc);
	$out->write_seq($seq_obj);
	print NAME "$name\t$loc\n";
}

close(IN);
$out->close();
close(NAME);
