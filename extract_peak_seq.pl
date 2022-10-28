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
my $outfile = shift or die $usage;

GetOptions(
"ext-len=i" => \$ext_len)
or die "Error in command line arguments, usage: $usage";

my $db = new Bio::DB::Fasta($db_path);
open(IN, "<$infile") || die "Unable to open $infile: $!";
my $out = new Bio::SeqIO(-file => ">$outfile", -format => 'fasta', -alphabet => 'dna');

# read BED file
my %name2loc;

while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end, $peak_name) = split(/\t/, $line);
	my ($name) = $peak_name =~ /Name=([^;]+)/;
	if(!exists $name2loc{$name}) { # not seen yet
		$name2loc{$name} = "$chr:$start-$end";
    if($ext_len > 0) {
      $start -= $ext_len;
      $end += $ext_len - 1;
			$start = 0 if($start < 0);
    }
		my $seq = $db->seq($chr, $start + 1 => $end);
		my $seq_obj = new Bio::Seq(-seq => $seq, -display_id => $name, -desc => $name2loc{$name});
		$out->write_seq($seq_obj);
	}
}

close(IN);
$out->close();
