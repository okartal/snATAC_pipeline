#!/bin/bash

seq_dir="/home/zhc268/scratch/seqdata/P56/demultiplexed/"
prefix="p56.rep1"
WORKDIR="/home/zhc268/scratch/others/"
OUTDIR="${WORKDIR}/${prefix}_test"

fastq1="${seq_dir}${prefix}.R1.fastq.gz"
fastq2="${seq_dir}${prefix}.R2.fastq.gz"
barcode="/home/zhc268/data/software/scATAC/barcodes/"
#genome="/projects/ps-epigen/GENOME/hg19/bowtie2_index/male.hg19.fa"
genome="/projects/ps-epigen/GENOME/mm10/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.fasta"

source activate bds_scATAC
mkdir -p $OUTDIR

cd $OUTDIR

scATAC.bds -v -threads 8 -r1 $fastq1 -r2 $fastq2 -barcode_dir $barcode  -make_barcode_mismatch 2 \
       -mark_duplicate $(which picard) -bowtie2_idx $genome -prefix $prefix -min_read 500 

