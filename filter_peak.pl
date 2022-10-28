#!/usr/bin/env perl
# Filter peaks based on primer vs. peak/clone flanking-seq alignments (fwd and rev)
use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SimpleAlign;

my $flank_size = 50;
my $seed_len = 10;
my $primer_len = 29;
my $max_seed_error = 0.1;
my $min_match = 12;

my $options = qq(Options:
  --flank-size INT: flanking size of the target sequence extracted for alignment [$flank_size]
  --seed-len INT: seed length of primer [$seed_len]
  --primer-len INT: overal length of primer [$primer_len]
	--max-seed-error FLOAT: max error of seed region (3' of primer, 5' of revcom) [$max_seed_error]
	--min-match FLOAT: min matched bases of the primer [$min_match]
);
  
my $usage = "Usage: $0 INFILE FWD-ALN REV-ALN BED-OUTFILE INFO-OUTFILE $options";

my $infile = shift or die $usage;
my $fwd_aln = shift or die $usage;
my $rev_aln = shift or die $usage;
my $bed_out = shift or die $usage;
my $info_out = shift or die $usage;

my $opts = join(" ", @ARGV);

GetOptions(
"flank-size=i" => \$flank_size,
"seed-len=i" => \$seed_len,
"primer-len=i" => \$primer_len,
"max-seed-error=f" => \$max_seed_error,
"min-match=i" => \$min_match)
or die "Error in command line arguments, usage: $usage";

if(!($seed_len > 0 && $primer_len >= $seed_len && 0 <= $max_seed_error && $max_seed_error <= 1 && $min_match > 0)) {
	print STDERR $options;
	exit;
}

open(IN, "<$infile") || die "Unable to open $infile: $!";
my $fwd_in = new Bio::AlignIO(-file => $fwd_aln, -format => 'emboss');
my $rev_in = new Bio::AlignIO(-file => $rev_aln, -format => 'emboss');
open(BED, ">$bed_out") || die "Unable to write to $bed_out: $!";
open(INFO, ">$info_out") || die "Unable to write to $info_out: $!";

# read in fwd and rev alignments
my @peakNames;
my %peak2strand;
my %peak2aln;

while(my $aln = $fwd_in->next_aln()) {
	my $peakName = $aln->get_seq_by_pos(2)->id();
	push(@peakNames, $peakName);
	$peak2strand{$peakName} = '+';
	$peak2aln{$peakName} = $aln;
}

while(my $aln = $rev_in->next_aln()) {
	my $peakName = $aln->get_seq_by_pos(2)->id;
	if(!exists $peak2aln{$peakName} || $aln->score() > $peak2aln{$peakName}->score()) {
		$peak2strand{$peakName} = '-';
		$peak2aln{$peakName} = $aln;
	}
}

# scan best alignments
my %peak2flag;
my %peak2attrs;

print INFO "peakName\tqname\tstrand\tflank_size\tseed_len\tprimer_len\talign_len\tqstart\tqend\ttstart\ttend\tseed_error\tprimer_match\tis_mispriming\n";
foreach my $peakName (@peakNames) {
  my $strand = $peak2strand{$peakName};
	my $aln = $peak2aln{$peakName};
  my $aln_len = $aln->length();
	my $GAP_CHAR = $aln->gap_char();

	my $query = $aln->get_seq_by_pos(1);
	my $qname = $query->id();
	my $qseq = uc $query->seq();
	my $qstart = $query->start();
	my $qend = $query->end();

	my $target = $aln->get_seq_by_pos(2);
	my $tname = $target->id();
	my $tseq = uc $target->seq();
	my $tstart = $target->start();
	my $tend = $target->end();

  my $seed_error = 0;
  my $match = 0;

	for(my $i = 1; $i <= $primer_len; $i++) {
		my $is_seed = $strand eq '+' && $i > $primer_len - $seed_len || $strand eq '-' && $i <= $seed_len;
# outside aligned region
		if($i < $qstart || $i > $qend) {
			$seed_error++ if($is_seed);
		}
# in alignment
		else {
			my $j = $aln->column_from_residue_number($qname, $i);
			my $qch = substr($qseq, $j - 1, 1);
			my $tch = substr($tseq, $j - 1, 1);
			next unless($qch ne $GAP_CHAR);
			if(substr($qseq, $j - 1, 1) eq substr($tseq, $j - 1, 1)) {
				$match++;
			}
			else {
				$seed_error++ if($is_seed);
			}
		}
	}

	my $flag = $seed_error / $seed_len <= $max_seed_error && $match >= $min_match ? 1 : 0;
	print INFO "$peakName\t$qname\t$strand\t$flank_size\t$seed_len\t$primer_len\t$aln_len\t$qstart\t$qend\t$tstart\t$tend\t$seed_error\t$match\t$flag\n";
	$peak2attrs{$peakName} = qq(PrimerStrand=$strand;AlignLen=$aln_len;QueryRange=$qstart-$qend;TargetRange=$tstart-$tend;SeedError=$seed_error;PrimerMatch=$match;);
# set mispriming flag
	$peak2flag{$peakName} = $flag;
}

# add headers
print BED "# Options invoked: $opts\n";
while(my $line = <IN>) {
	chomp $line;
	my @fields = split(/\t/, $line);
	my ($name) = $fields[3] =~ /Name=([^;]+)/;
	if(!$peak2flag{$name}) {
		$fields[3] .= $peak2attrs{$name};
		print BED join("\t", @fields), "\n";
	}
}

close(IN);
close(BED);
close(INFO);
