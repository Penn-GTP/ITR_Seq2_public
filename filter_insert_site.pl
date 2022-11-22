#!/usr/bin/env perl
# Filter insert site based on primer vs. site/clone flanking-seq alignments (fwd and rev)
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;

my $flank_size = 50;
my $seed_len = 10;
my $primer_file = 'ITR_primer.fa';
my $max_seed_error = 0.1;
my $min_match = 12;

my $options = qq(Options:
  --flank-size INT: flanking size of the target sequence extracted for alignment [$flank_size]
  --seed-len INT: seed length of primer [$seed_len]
  --primer-file FILE: primer (fwd) file [$primer_file]
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
"primer-file=s" => \$primer_file,
"max-seed-error=f" => \$max_seed_error,
"min-match=i" => \$min_match)
or die "Error in command line arguments, usage: $usage";

if(!($seed_len > 0 && -e $primer_file && 0 <= $max_seed_error && $max_seed_error <= 1 && $min_match > 0)) {
	print STDERR $options;
	exit;
}

my $primer_in = new Bio::SeqIO(-file => "<$primer_file", -format => 'fasta');
open(IN, "<$site_in") || die "Unable to open $site_in: $!";
open(NAME, "<$name_in") || die "Unable to open $name_in: $!";
my $fwd_in = new Bio::AlignIO(-file => $fwd_aln, -format => 'emboss');
my $rev_in = new Bio::AlignIO(-file => $rev_aln, -format => 'emboss');
open(SITE, ">$site_out") || die "Unable to write to $site_out: $!";
open(INFO, ">$info_out") || die "Unable to write to $info_out: $!";

# read in primer length
my @primer_names;
my %primer2len;
while(my $seq_obj = $primer_in->next_seq()) {
	my $primer_id = $seq_obj->id();
	push(@primer_names, $primer_id);
	$primer2len{$primer_id} = $seq_obj->length();
}

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
my %aln2strand;
my %aln2info;

while(my $aln = $fwd_in->next_aln()) {
	my $qname = $aln->get_seq_by_pos(1)->id();
	my $tname = $aln->get_seq_by_pos(2)->id();
	$aln2strand{$qname}{$tname} = '+';
	$aln2info{$qname}{$tname} = $aln;
}

while(my $aln = $rev_in->next_aln()) {
	my $qname = $aln->get_seq_by_pos(1)->id();
	my $tname = $aln->get_seq_by_pos(2)->id();
	if(!exists $aln2info{$qname}{$tname} || $aln->score() > $aln2info{$qname}{$tname}->score()) {
		$aln2strand{$qname}{$tname} = '-';
		$aln2info{$qname}{$tname} = $aln;
	}
}

# scan best alignments
my %loc2flag;
print INFO "siteName\tloc\tqname\tqstrand\tflank_size\tseed_len\tqlen\talign_len\tqstart\tqend\ttstart\ttend\tscore\tseed_error\tprimer_match\tis_mispriming\n";
foreach my $siteName (@siteNames) {
	my $loc = $name2loc{$siteName};
	foreach my $qname (@primer_names) {
		my $strand = $aln2strand{$qname}{$siteName};
		my $aln = $aln2info{$qname}{$siteName};
		next unless(defined $aln);

		my $aln_len = $aln->length();
		my $GAP_CHAR = $aln->gap_char();

		my $query = $aln->get_seq_by_pos(1);
		my $qlen = $primer2len{$qname};
		my $qseq = uc $query->seq();
		my $qstart = $query->start();
		my $qend = $query->end();

		my $target = $aln->get_seq_by_pos(2);
		my $tname = $target->id();
		my $tseq = uc $target->seq();
		my $tstart = $target->start();
		my $tend = $target->end();

    my $score = $aln->score();
		my $seed_error = 0;
		my $match = 0;

# scan pre-alignment region
		for(my $i = 1; $i < $qstart; $i++) {
			my $is_seed = $strand eq '+' && $i > $qlen - $seed_len || $strand eq '-' && $i <= $seed_len;
			$seed_error++ if($is_seed);
		}

# scan in-alignment region
		for(my ($k, $i, $j) = (0, $qstart, $tstart); $k < $aln_len; $k++) {
			my $is_seed = $strand eq '+' && $i > $qlen - $seed_len || $strand eq '-' && $i <= $seed_len;

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
		for(my $i = $qend + 1; $i <= $qlen; $i++) {
			my $is_seed = $strand eq '+' && $i > $qlen - $seed_len || $strand eq '-' && $i <= $seed_len;
			$seed_error++ if($is_seed);
		}

		my $flag = $seed_error / $seed_len <= $max_seed_error && $match >= $min_match ? 1 : 0;
		print INFO "$siteName\t$loc\t$qname\t$strand\t$flank_size\t$seed_len\t$qlen\t$aln_len\t$qstart\t$qend\t$tstart\t$tend\t$score\t$seed_error\t$match\t$flag\n";
# set mispriming flag
		$loc2flag{$loc} |= $flag;
	}
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
