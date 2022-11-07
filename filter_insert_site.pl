#!/usr/bin/env perl
# Filter insert site based on primer vs. site/clone flanking-seq alignments (fwd and rev)
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

my $site_in = shift or die $usage;
my $name_in = shift or die $usage;
my $fwd_aln = shift or die $usage;
my $rev_aln = shift or die $usage;
my $site_out = shift or die $usage;
my $info_out = shift or die $usage;

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

open(IN, "<$site_in") || die "Unable to open $site_in: $!";
open(NAME, "<$name_in") || die "Unable to open $name_in: $!";
my $fwd_in = new Bio::AlignIO(-file => $fwd_aln, -format => 'emboss');
my $rev_in = new Bio::AlignIO(-file => $rev_aln, -format => 'emboss');
open(SITE, ">$site_out") || die "Unable to write to $site_out: $!";
open(INFO, ">$info_out") || die "Unable to write to $info_out: $!";

# read in name2loc
my @siteNames;
my %name2loc;
while(my $line = <NAME>) {
	chomp $line;
	my ($name, $loc) = split(/\t/, $line);
	push(@siteNames, $name);
	$name2loc{$name} = $loc;
}

# read in fwd and rev alignments
my %site2strand;
my %site2aln;

while(my $aln = $fwd_in->next_aln()) {
	my $siteName = $aln->get_seq_by_pos(2)->id();
	$site2strand{$siteName} = '+';
	$site2aln{$siteName} = $aln;
}

while(my $aln = $rev_in->next_aln()) {
	my $siteName = $aln->get_seq_by_pos(2)->id;
	if(!exists $site2aln{$siteName} || $aln->score() > $site2aln{$siteName}->score()) {
		$site2strand{$siteName} = '-';
		$site2aln{$siteName} = $aln;
	}
}

# scan best alignments
my %loc2flag;
print INFO "siteName\tloc\tqstrand\tflank_size\tseed_len\tprimer_len\talign_len\tqstart\tqend\ttstart\ttend\tseed_error\tprimer_match\tis_mispriming\n";
foreach my $siteName (@siteNames) {
	my $loc = $name2loc{$siteName};
  my $strand = $site2strand{$siteName};
	my $aln = $site2aln{$siteName};
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

  # scan pre-alignment region
	for(my $i = 1; $i < $qstart; $i++) {
		my $is_seed = $strand eq '+' && $i > $primer_len - $seed_len || $strand eq '-' && $i <= $seed_len;
		$seed_error++ if($is_seed);
	}

  # scan in-alignment region
	for(my ($k, $i, $j) = (0, $qstart, $tstart); $k < $aln_len; $k++) {
		my $is_seed = $strand eq '+' && $i > $primer_len - $seed_len || $strand eq '-' && $i <= $seed_len;
		
		my $qch = substr($qseq, $k, 1);
		my $tch = substr($tseq, $k, 1);

		if($qch eq $tch) {
			$match++;
		}
		else {
			$seed_error++ if($is_seed);
		}
		$i++ if($qch ne $GAP_CHAR);
		$j++ if($tch ne $GAP_CHAR);
	}

  # scan post-alignment region
	for(my $i = $qend + 1; $i <= $primer_len; $i++) {
		my $is_seed = $strand eq '+' && $i > $primer_len - $seed_len || $strand eq '-' && $i <= $seed_len;
		$seed_error++ if($is_seed);
	}

	my $flag = $seed_error / $seed_len <= $max_seed_error && $match >= $min_match ? 1 : 0;
	print INFO "$siteName\t$loc\t$strand\t$flank_size\t$seed_len\t$primer_len\t$aln_len\t$qstart\t$qend\t$tstart\t$tend\t$seed_error\t$match\t$flag\n";
# set mispriming flag
	$loc2flag{$loc} = $flag;
}

# output insert site filtered
while(my $line = <IN>) {
	chomp $line;
	my ($chr, $start, $end) = split(/\t/, $line);
	my $loc = "$chr:$start-$end";
	if(!$loc2flag{$loc}) {
		print SITE "$line\n";
	}
}

close(IN);
close(NAME);
close(SITE);
close(INFO);
