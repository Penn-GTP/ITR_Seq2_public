package MiSeqITRSeqExpDesign;
use strict;
use warnings;
use File::Basename;

our $VERSION = v1.1;
# This module is used to do option parsing given a filehandle of an experimental design file of ITR_Seq pipeline
# Author: Qi Zheng
# Since: 02/01/2022

# GLOBAL options and default values
our %GLOBAL_OPTS = (
  RUN_STAT => 'ITR_Seq_run_stats.tsv',
  MAX_PROC => 8,
	BASE_DIR => '.',
	WORK_DIR => 'WORK',
	SCRIPT_DIR => 'scripts',
	DEMUX_DIR => 'fastqs/demux',
# FASTQ_DIR => 'FASTQ',
#	SRA_DIR => 'SRA_submission',
	VECTOR_DIR => 'AAV_vector',
	UMI_LEN => 8,
	UMI_MM => 0,
	INSERT_SIZE => 20,
	KEEP_UNPAIR => 1,
	KEEP_STRAND => 3,
	MAX_PEAK_DIST => 44,
	MAX_CLONE_DIST => -2
);
  
# Constructor taking a filehandle or a filename
sub new {
  my $class = shift;
  my $file = shift;
  my $self = { };
  my %global_opt;
  my @sample_names;  # all samples that include these SampleIDs
  my @opt_names;  # per-sample option names
  my %sample_opt;
  # Set default global opts
	while(my ($key, $val) = each %GLOBAL_OPTS) {
		$global_opt{$key} = $val;
	}

  # read the experimental design file
  if(-f $file || -l $file) { # if $file is a file or link
		open(IN, "<$file") || die "Unable to open $file: $!";
  }
  elsif(ref $file eq 'IO') { # is a filehandle reference
	  *IN = $file; # file is a filehandle
  }
  else {
	  print STDERR "new() must take a filename or a filehandle!\n";
	  exit;
	}
  while(my $line = <IN>) {
	  chomp $line;
	  if($line =~ /^#/) {
		  if($line =~ /^## (\w+)=(\S+):/) { # global opt
			  $global_opt{$1} = $2;
		  }
		  elsif($line =~ /^## (\w+):/) { # opt name description line
			  push(@opt_names, $1);
		  }
		  else {
			  next; # ignore
		  }
	  }
	  else { # opt value line
		  my @values = split(/\t/, $line, -1);
		  if(@opt_names != @values) {
			  print STDERR "Incorrect field number at line $.: $line\n",
				  "Found ", scalar @values, " fields, but required ", scalar @opt_names, " fields\n";
			  exit;
		  }
		  my $sample = $values[0];
		  push(@sample_names, $sample);
# pair opt names with values
		  for(my $i = 0; $i < @opt_names; $i++) {
			  my $val = $values[$i];
			  $val = '' if(!defined $val);
				$sample_opt{$sample}{$opt_names[$i]} = $val;
			}
		}
	} # end each line of experiment design file

	if(-f $file) {
		close(IN);
	}

# Record variables
	%{$self->{'global_opt'}} = %global_opt;
	@{$self->{'sample_names'}} = @sample_names;
	@{$self->{'opt_names'}} = @opt_names;
	%{$self->{'sample_opt'}} = %sample_opt; 
	return bless $self, $class;
}

# Method to get global opt given opt name
sub get_global_opt {
	my ($self, $name) = @_;
	return $self->{'global_opt'}{$name};
}

# Method to get all sample_names
sub get_sample_names {
	my $self = shift;
	return @{$self->{'sample_names'}};
}

# Method to get all opt_names
sub get_opt_names {
	my $self = shift;
	return @{$self->{'opt_names'}};
}

# Method to get or set single per-sample opt value
sub sample_opt {
	my $self = shift;
	my $sample = shift;
	my $opt = shift;
	if(@_) {
		$self->{'sample_opt'}{$sample}{$opt} = shift;
	}
	return $self->{'sample_opt'}{$sample}{$opt};
}

# Method to get one or more per-sample opt values
sub get_sample_opts {
	my ($self, $sample, @opt_names) = @_;
	my @values;
	foreach my $opt (@opt_names) {
		push(@values, $self->{'sample_opt'}{$sample}{$opt});
	}
	return @values;
}

# get per-sample UMI labled forward fastq output file
sub get_sample_fwd_UMI_file {
	my ($self, $sample) = @_;
	return "$sample\_R1_UMI.fastq.gz";
}

# get per-sample UMI labled reverse fastq output file
sub get_sample_rev_UMI_file {
	my ($self, $sample) = @_;
	return "$sample\_R2_UMI.fastq.gz";
}

# get per-sample forward ITR-seq trimmed fastq output file
sub get_sample_fwd_ITRtrim_file {
	my ($self, $sample) = @_;
	return "$sample\_R1_trimmed.fastq.gz";
}

# get per-sample reverse ITR-seq trimmed fastq output file
sub get_sample_rev_ITRtrim_file {
	my ($self, $sample) = @_;
	return "$sample\_R2_trimmed.fastq.gz";
}

# get per-sample forward ITR-seq untrimmed fastq output file
sub get_sample_fwd_ITRuntrim_file {
	my ($self, $sample) = @_;
	return "$sample\_R1_untrimmed.fastq.gz";
}

# get per-sample reverse ITR-seq trimmed fastq output file
sub get_sample_rev_ITRuntrim_file {
	my ($self, $sample) = @_;
	return "$sample\_R2_untrimmed.fastq.gz";
}

# get per-sample forward ITR-seq short fastq output file
sub get_sample_fwd_ITRshort_file {
	my ($self, $sample) = @_;
	return "$sample\_R1_short.fastq.gz";
}

# get per-sample reverse ITR-seq short fastq output file
sub get_sample_rev_ITRshort_file {
	my ($self, $sample) = @_;
	return "$sample\_R2_short.fastq.gz";
}

# get per-sample ref map file
sub get_sample_ref_map_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map.bam";
}

# get per-sample ref filtered sorted file
sub get_sample_ref_filtered_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_filtered_sorted.bam";
}

# get per-sample ref deduplicated file
sub get_sample_ref_dedup_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_filtered_sorted_dedup.bam";
}

# get per-sample ref deduplicated log
sub get_sample_ref_dedup_log {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_filtered_sorted_dedup.log";
}

# get per-sample ref novec file
sub get_sample_ref_novec_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_filtered_sorted_dedup_novec.bam";
}

# get per-sample vec map file
sub get_sample_vec_map_file {
	my ($self, $sample) = @_;
	return "$sample\_vec_map.bam";
}

# get per-sample vec filtered file
sub get_sample_vec_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_vec_map_filtered.bam";
}

# get per-sample vec sorted file
sub get_sample_vec_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_vec_map_filtered_sorted.bam";
}

# get per-sample vec id file
sub get_sample_vec_ID_file {
	my ($self, $sample) = @_;
	return "$sample\_vec_map_filtered_sorted_ID.txt";
}

# get per-sample peak file
sub get_sample_ref_peak {
	my ($self, $sample) = @_;
	return "$sample\_ref_peak.bed";
}

# get per-sample sorted peak file
sub get_sample_ref_sorted_peak {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_peak.bed";
}

# get per-sample ref merged peak file
sub get_sample_ref_merged_peak {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_peak.bed";
}

# get per-sample filtered peak file
sub get_sample_ref_filtered_peak {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_filtered_peak.bed";
}

# get per-sample ref merged clone file
sub get_sample_ref_merged_clone {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_clone.bed";
}

# get per-sample filtered clone file
sub get_sample_ref_filtered_clone {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_filtered_clone.bed";
}

# get per-sample vec seq file
sub get_sample_vec_seq {
	my ($self, $sample) = @_;
	my $name = basename($self->sample_opt($sample, 'vector_file'), qw(.gb .gbk));
	return "$name\_vec_seq.fasta";
}

# get per-sample vec annotation
sub get_sample_vec_anno {
	my ($self, $sample) = @_;
	my $name = basename($self->sample_opt($sample, 'vector_file'), qw(.gb .gbk));
	return "$name\_vec_anno.gff3";
}

# get per-sample vec seq masked file
sub get_sample_vec_seq_masked {
	my ($self, $sample) = @_;
	my $name = basename($self->sample_opt($sample, 'vector_file'), qw(.gb .gbk));
	return "$name\_vec_seq_masked.fasta";
}

# get per-sample vec dbname
sub get_sample_vec_dbname {
	my ($self, $sample) = @_;
	my $name = basename($self->sample_opt($sample, 'vector_file'), qw(.gb .gbk));
	return "$name\_vec_seq_masked";
}

# get per-sample ref peak track file
sub get_sample_ref_peak_track {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_filtered_peak_track.bed";
}

# get per-sample ref clone track file
sub get_sample_ref_clone_track {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_filtered_clone_track.bed";
}

# get per-sample clone info
sub get_sample_ref_clone_info {
  my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_filtered_clone_info.txt";
}

# get per-sample ref peak seq
sub get_sample_ref_peak_seq {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_seq.fasta";
}

# get per-sample ref peak anno
sub get_sample_ref_peak_anno {
	my ($self, $sample) = @_;
	return "$sample\_ref_sorted_merged_filtered_peak_anno.tsv";
}

# get per-exp stats
sub get_exp_stats_file {
  my ($self, $exp_file) = @_;
  my $stats_file = basename($exp_file, qw(.conf .txt .tsv));
  return "$stats_file\_sample_stats.tsv";
}

# get per-exp stats
sub get_exp_clone_loc_distrib_fig {
  my ($self, $exp_file) = @_;
  my $stats_file = basename($exp_file, qw(.conf .txt .tsv));
  return "$stats_file\_sample_clone_loc_distrib.pdf";
}

1;
