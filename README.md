# ITR_Seq2 pipeline manual

This pipeline guideline was initiated by: Qi Zheng (<zhengqi@pennmedicine.upenn.edu>).

This manual is for the updated ITR_Seq2 pipeline.

We thank Kelly Martins (<kelly.martins@pennmedicine.upenn.edu>) and Jenny Greig (<greg@upenn.edu>) for their help of understanding the original ITR_Seq pipeline.

ITR sequencing assay is a method to identify genome-wide DNA editing sites *in vivo* following the adeno-associated viral vector-mediated genome editing.

---

## Citations
Please cite the following research if you are using this ITR_Seq pipeline.
> Breton, C., et.al., *ITR-Seq, a next-generation sequencing assay, identifies genome-wide DNA editing sites in vivo following adeno-associated viral vector-mediated genome editing.*
[BMC Genomics. 2020 Mar 17;21(1):239.](https://pubmed.ncbi.nlm.nih.gov/32183699/)

---

## Dependency and pre-requirement
You need working Perl and Java distributions and available in your *PATH* to run this pipeline.
Both Perl and Java are generally available under most model Operation Systems,
including Windows, Mac OS X, and Unix/Linux, and is likely pre-install on most Unix/Linux based systems.

Besides, this pipeline also depends the following programs/tools pre-installed and have them available in your PATH.

If you are running this pipeline on PMACS HPC, all the following dependencies are pre-installed and pre-configured,
so you can ignore this step and jump to the next step.

1. **Samtools** - basic tool for manipulating standard alignment files in SAM/BAM format.
If you are running this pipeline on PMACS HPC, it is already installed and configured and 
You can download and install **Samtools** by following the instructions at <https://www.htslib.org/>.

2. **Bowtie2** - Default NGS read aligner used by this pipeline.
You can download and install **Bowtie2** by following the instructions at <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>.

3. **Cutadapt** - Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
You can download and install **Cutadapt** by following the instructions at <https://cutadapt.readthedocs.io/>.

4. **Picard tools** - Picard is a set of Java based command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
This pipeline ships a JAR copy of **Picard tools** with it, but you can replace it with your own copy, you can achieve this by
simply download and replace the latest `picard.jar` file from https://broadinstitute.github.io/picard/ in this pipeline.

5. **Je** - tool suite used for removing PCR/UMI duplication (dedup) in NGS alignments.
You can download and install **Je** by following the instructions at <https://github.com/gbcs-embl/Je>.
 
## Getting Started

This pipeline's scripts only "prepare" the running commands (cmds) and save them into the shell/bash files.
After running these "prepare" scripts, you will need to run them either by directly typing their names,
or "submitting" them via an LSF system such as `bsub` that is available on PMACS/HPC.
The "prepare" scripts will also check whether the output files already exist,
and if yes, it will not override the existing files and instead generate "comment-out" commands that you can feel free to enable/disable them manually.

**Note:** If you are using this pipeine in an HPC environment, be sure to do any of these steps in an interactive node,
   initiated by typing `bsub -Is "bash"`.

### Create the experimental design config file

You need to create a special tab-delimited config file to specify all the experimental designs for a given project (**PROJECT_ID**)
You should do this by copying the template file `Template_ITR_Seq_experimental_design.conf` in this pipeline into your working directory,
and fill/edit it with a text-editor or Excel/Libre-Office (choose csv/tsv mode),
and save it to your working directory (assuming it is saved as `PROJECT_ID_EXP_DESIGN.conf`)
 
Each project-specific config file contains global options and per-sample local options,
both were explained as comment lines at the beginning of the template,
while the global options are in all upper-cases, and the values are given in **GLOBAL\_OPT=GLOBAL\_VAL** format.

**Note:** If you are setting any of the working dir options below to values other than the current dir (./),
you will have to create or link the directories manually before running this pipeline.

#### Global options

- **NUM\_PROC**: # of processors (cores/threads) to use [default 8]
- **BASE\_DIR**: base directory for generating major result files (sequences, alignments, tables, figures, etc.) [default ./]
- **WORK\_DIR**: working directory for all intermediate files [default WORK/]
- **SCRIPT\_DIR**: path to this pipeline feel free to link it to your working directory [default scripts/]
- **DEMUX\_DIR**: path to demultiplexed per-sample FASTQ read files [default fastqs/demux/]
- **VECTOR\_DIR**: path to the AAV vector/plasmid map files in GenBank format, also used to build vector databases [default AAV_vector/]
- **UMI\_LEN**: Unique Molecular Identifier (UMI) length that were built into the P5 Index primer [default 8]
- **UMI\_MM**: mismatch allowed when mark UMI duplicates [default 0]
- **ITR\_PRIMER**: ITR-specific primer sequence that are used as the P7 primer for ITR-Seq
- **INSERT\_SIZE**: ITR insertion recognition site size by the meganuclease, used to call insertion site [default 2]
- **KEEP\_UNPAIR**: beside paired mapped reference reads, which additional strand to keep, 0 for none, 1 for forward, 2 for reverse and 3 for both [default 1]
- **KEEP\_STRAND**: strand(s) required for peaks, 0: no requirement, 1 forward, 2: reverse, 3: both, we recommend use 3 for gene-editing and 0 for gene-therapy [default 3]
- **MAX\_PEAK\_DIST**: maximum distance allowed for merging multiple insertion sites into peaks, we recommend use 2 X known target size [default 44]
- **MAX\_CLONE\_DIST**: maximum distance allowed for merging unique clone sites, ignored for gene-editing samples,
use negative values to enforce overlapping; We recommend -INSERT_SIZE [default -2]

#### Per-sample options
- **sample\_name**: unique sample name
- **fastq\_R1**: filename of forward FASTQ read (R1), .gz or .bz2 files are accepted
- **fastq\_R2**: filename of reverse FASTQ read (R2), .gz or .bz2 files are accepted
- **fastq\_I1**: filename of forward FASTQ index (I1), .gz or .bz2 files are accepted
- **fastq\_I2**: filename of reverse FASTQ index (I2), .gz or .bz2 files are accepted

- **vector\_file**: Vector sequence(s) file in GenBank format
- **trim\_prog**: program used to identify ITR-primer containing reads (5' of R2, 3' of R1), now only supports 'cutadapt'
- **max\_error\_rate**: maximum error rate allowed for identifying ITR-primer containing reads (recommend 0.1)
- **min\_len**: minimum length after trimming (recommend 18)
- **trim\_opts**: additional options to invoke the trimming program, or leave blank if none
- **aligner**: NGS aligner to use, now supports 'bowtie2' and 'bwa'
- **align\_opts**: additional options to invoke the NGS aligner, or leave blank if none
- **ref\_db**: path/name to pre-built reference (host) genome database, e.g.: `/project/gtplab/pub_data/genomes/Macaca_mulatta/Bowtie2_index/Mmul_10/MMul_10`
- **min\_mapQ**: minimum mapping quality (mapQ) values required (recommend 30)
- **ref\_gff**: path to reference annotation file(s) in GFF3/GTF format, e.g.: /project/gtplab/pub_data/genomes/Macaca_mulatta/annotation/Macaca_mulatta.Mmul_10.105.clean.gff3, multiple files can be separated by space
- **genome\_seq**: path to reference seq file/dir, e.g. /project/gtplab/pub_data/genomes/Macaca_mulatta/fasta/Mmul_10/
- **target\_file**: path to the known gene editing target region(s) in BED format, leave blank if this is a gene therapy sample
- **min\_clone\_loc**: minimum distinct UMI R1 locus required to define a clone expansion, ignored if `target_file` is provided (gene-editing sample) (default 2)

---

## Step REF -- Build Reference Database files for the project
- Prepare REF cmds by running:
`SCRIPT_DIR/prepare_ref_cmd.pl PROJECT_ID_EXP_DESIGN.conf ref.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `ref.sh`
   
- Run REF cmds in `ref.sh` by running:
   - On a Linux cluster: `bsub -J REF -o ref.log ./ref.sh`
   - On a regular Linux server: `./ref.sh > ref.log 2>&1`

- Output: REF database files in `VECTOR_DIR`

## Step UMI -- Append UMI information from index 2 read (I2) to the read names of de-multiplexed FASTQ files
- Prepare UMI cmds by running:
`SCRIPT_DIR/prepare_UMI_cmd.pl PROJECT_ID_EXP_DESIGN.conf UMI.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `UMI.sh`
   
- Run UMI cmds in `UMI.sh` by running:
   - On a Linux cluster (recommended): `bsub -J UMI -o UMI.log ./UMI.sh`

- Output for each **SAMPLE**: `WORK_DIR/SAMPLE_R1_UMI.fastq.gz` and `WORK_DIR/SAMPLE_R2_UMI.fastq.gz`

## Step TRIM -- Trim (identify) ITR-primer containing reads
- Prepare TRIM cmds by running:
`SCRIPT_DIR/prepare_trim_cmd.pl PROJECT_ID_EXP_DESIGN.conf trim.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `trim.sh`
   
- Run TRIM cmds in `trim.sh` by running:
   - On a Linux cluster (recommended): `bsub -J TRIM -o trim.log -n 24 -M 16000 ./trim.sh`
 
- Output for each **SAMPLE**:
       - ITR-primer containing FASTQ files: `WORK_DIR/SAMPLE_R1_trimmed.fastq.gz` and `WORK_DIR/SAMPLE_R2_trimmed.fastq.gz`
       - ITR-primer non-containing FASTQ files: `WORK_DIR/SAMPLE_R1_untrimmed.fastq.gz` and `WORK_DIR/SAMPLE_R2_untrimmed.fastq.gz`
       - Too short reads after ITR-primer trimming: `WORK_DIR/SAMPLE_R1_short.fastq.gz` and `WORK_DIR/SAMPLE_R2_short.fastq.gz`
       
## Step MAP -- Map ITR-containing reads to host and vector sequences
- Prepare MAP cmds by running:
`SCRIPT_DIR/prepare_map_cmd.pl PROJECT_ID_EXP_DESIGN.conf map.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `map.sh`
   
- Run MAP cmds in `map.sh` by running:
   - On a Linux cluster (recommended): `bsub -J MAP -o map.log -n 24 -M 24000 ./map.sh`
   
- Output for each **SAMPLE**:
       - Host reference (ref) mapped alignment file: `WORK_DIR/SAMPLE_ref_map.bam`
       - Vector (vec) mapped alignment file: `WORK_DIR/SAMPLE_vec_map.bam`
       
## Step FILTER -- Filter ref and vec mapped alignments
- Prepare FILTER cmds by running:
`SCRIPT_DIR/prepare_filter_cmd.pl PROJECT_ID_EXP_DESIGN.conf filter.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `filter.sh`
   
- Run FILTER cmds in `filter.sh` by running:
   - On a Linux cluster: `bsub -J FILTER -o filter.log -n 4 -M 24000 ./filter.sh`
   - On a regular Linux server: `./filter.sh > filter.log 2>&1`
   
- Output for each **SAMPLE**:
       - ref mapped, deduplexed and non-vec mapped alignment file: `BASE_DIR/SAMPLE_ref_map_filtered_sorted_dedup_novec.bam`

## Step PEAK -- call AAV insert sites/peaks
- Prepare PEAK cmds by running:
`SCRIPT_DIR/prepare_peak_cmd.pl PROJECT_ID_EXP_DESIGN.conf peak.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `peak.sh`
   
- Run PEAK cmds in `peak.sh` by running:
   - On a Linux cluster: `bsub -J PEAK -o peak.log -n 4 -M 24000 ./peak.sh`
   - On a regular Linux server: `./peak.sh > peak.log 2>&1`
   
- Output for each **SAMPLE**:
       - filtered peaks in BED format: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_peak.bed`
       - (Optional) for gene-therapy only: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_clone.bed`

## Step ANNOTATE -- annotate called peaks/clones
- Prepare ANNO cmds by running:
`SCRIPT_DIR/prepare_annotate_cmd.pl PROJECT_ID_EXP_DESIGN.conf annotate.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `annotate.sh`
   
- Run ANNOTATE cmds in `annotate.sh` by running:
   - On a Linux cluster: `bsub -J ANNOTATE -o annotate.log -n 4 -M 24000 ./annotate.sh`
   - On a regular Linux server: `./annotate.sh > annotate.log 2>&1`
   
- Output for each **SAMPLE**:
   - peak track file in BED format for being displayed in IGV: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_peak_track.bed`
   - peak annotation file in BED format with overlapping gene annotations: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_peak_anno.bed`
   - peak info file in TSV format: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_peak_info.txt`
   - peak seq file in FASTA format: `BASE_DIR/SAMPLE_ref_peak_seq.fasta`
   - (optional for gene-therapy only) clone track file in BED format for being displayed in IGV: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_clone_track.bed`
   - (optional for gene-therapy only) clone annotation file in BED format with overlapping gene annotations: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_clone_anno.bed`
   - (optional for gene-therapy only) clone info file: `BASE_DIR/SAMPLE_ref_sorted_merged_filtered_clone_info.txt`
   
- Output for this entire experiment/run:
   - sample statistics summary in TSV format: `PROJECT_ID_EXP_DESIGN_sample_stats.tsv`

 
## Notes:
- Common options can be used for the LSF `bsub` system, try `man bsub` for more options and details:
   - `-J`: Job name (for display only)
   - `-o`: redirect job logs (stdout and stderr) to the given file 
   - `-n`: # of cores/CPUs requests for this job
   - `-M`: memory required for this job

---
