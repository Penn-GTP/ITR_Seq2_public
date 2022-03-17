#!/bin/bash

## Load required modules on HPC
module load java/openjdk-1.8.0
module load samtools/1.11 # SAMtools
#module load bowtie2/2.3.4.1 # Bowtie2 aligner
#module load bwa-0.7.10 # BWA aligner
#module load picard/1.96 # PICARD tool for manipulating SAM/BAM files
module load R/4.0.2 # R for statistics and visulization

# set envs
export PATH=/project/gtplab/apps/bin:$PATH
#export PATH=$PATH:/project/gtplab/apps/cutadapt-3.4/bin
