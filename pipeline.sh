#!/bin/sh

if [ "$#" != "5" ]; then
  echo "[!] Usage : ./pipeline.sh <reference.fa> <dataset-directory> <normal.bam> <tumor.bam> <number of threads>"
  exit
fi

REFERENCE_FASTA=$1
REFERENCE_BASE=$(basename $1)
DATASET_DIR=$2
NORMAL_BAM=$3
TUMOR_BAM=$4
NUM_THREADS=$5

SCRIPT_DIR=$(cd $(dirname $0); pwd)

echo "!! Note that samtools/bwa/mafft needs to be installed !!"
echo "This program creates intermediate files in the $DATASET_DIR during the pipeline."
date

### indexing for bwa
if [ ! -e "$REFERENCE_FASTA.index.amb" ] || [ ! -e "$REFERENCE_FASTA.index.ann" ] || [ ! -e "$REFERENCE_FASTA.index.bwt" ] || [ ! -e "$REFERENCE_FASTA.index.pac" ] || [ ! -e "$REFERENCE_FASTA.index.sa" ]; then
  if [ ! -e "$REFERENCE_FASTA" ]; then
    echo "[!] $REFERENCE_FASTA not found. "
    exit
  fi
  echo "[pipeline] $REFERENCE_FASTA.index is not found. Indexing process started ... "
  bwa index -p $REFERENCE_FASTA.index $REFERENCE_FASTA
  echo "[pipeline] bwa index Done."
  date
fi

if [ ! -e "$REFERENCE_FASTA.fixed.fa" ]; then
  if [ ! -e "$REFERENCE_FASTA" ]; then
    echo "[!] $REFERENCE_FASTA not found. "
    exit
  fi
  echo "[pipeline] $REFERENCE_FASTA.fixed.fa is not found. Fixing process started ... "
  $SCRIPT_DIR/bin/fix_reference_for_indexing $REFERENCE_FASTA $REFERENCE_FASTA.fixed.fa
  echo "[pipeline] fix_reference_for_indexing Done."
  date
fi

### filter read from BAM
if [ ! -e "$DATASET_DIR/$NORMAL_BAM.filtered.fastq" ]; then
  if [ ! -e "$DATASET_DIR/$NORMAL_BAM" ]; then
    echo "[!] $DATASET_DIR/$NORMAL_BAM not found. "
    exit
  fi
  echo "[pipeline] $DATASET_DIR/$NORMAL_BAM.filtered.fastq is not found. Filtering process started ... "
  $SCRIPT_DIR/scripts/filter_read_from_sam.py $DATASET_DIR/$NORMAL_BAM $DATASET_DIR/$NORMAL_BAM.filtered.bam
  samtools fastq $DATASET_DIR/$NORMAL_BAM.filtered.bam > $DATASET_DIR/$NORMAL_BAM.filtered.fastq
  echo "[pipeline] filter_read_from_sam (normal) Done."
  date
fi

if [ ! -e "$DATASET_DIR/$TUMOR_BAM.filtered.fastq" ]; then
  if [ ! -e "$DATASET_DIR/$TUMOR_BAM" ]; then
    echo "[!] $DATASET_DIR/$TUMOR_BAM not found. "
    exit
  fi
  echo "[pipeline] $DATASET_DIR/$TUMOR_BAM.filtered.fastq is not found. Filtering process started ... "
  $SCRIPT_DIR/scripts/filter_read_from_sam.py $DATASET_DIR/$TUMOR_BAM $DATASET_DIR/$TUMOR_BAM.filtered.bam
  samtools fastq $DATASET_DIR/$TUMOR_BAM.filtered.bam > $DATASET_DIR/$TUMOR_BAM.filtered.fastq
  echo "[pipeline] filter_read_from_sam (tumor) Done."
  date
fi

### indexing
if [ ! -e "$DATASET_DIR/$REFERENCE_BASE.$NORMAL_BAM.index.cst" ] || [ ! -e "$DATASET_DIR/$REFERENCE_BASE.$NORMAL_BAM.index.refbv" ]; then
  echo "[pipeline] $DATASET_DIR/$REFERENCE_BASE.$NORMAL_BAM.index is not found. Indexing process started ... "
  $SCRIPT_DIR/bin/construct_index $REFERENCE_FASTA.fixed.fa $DATASET_DIR/$NORMAL_BAM.filtered.fastq $DATASET_DIR/$REFERENCE_BASE.$NORMAL_BAM.index $NUM_THREADS
  rm $DATASET_DIR/$NORMAL_BAM.filtered.fastq_reference_read_text.tmp
  echo "[pipeline] construct_index Done."
  date
fi

### alignment on graph
if [ ! -e "$DATASET_DIR/$TUMOR_BAM.extracted.fastq" ]; then
  echo "[pipeline] Main filtering process started ... "
  $SCRIPT_DIR/bin/rf_mapping $DATASET_DIR/$TUMOR_BAM.filtered.fastq $DATASET_DIR/$REFERENCE_BASE.$NORMAL_BAM.index $DATASET_DIR/$TUMOR_BAM.extracted.fastq $NUM_THREADS
  echo "[pipeline] pruning_read Done."
  date
fi

### report on reference
if [ ! -e "$DATASET_DIR/$TUMOR_BAM.result.bedpe" ]; then
  echo "[pipeline] Creating report ... "
  bwa mem $REFERENCE_FASTA.index $DATASET_DIR/$TUMOR_BAM.extracted.fastq > $DATASET_DIR/$TUMOR_BAM.extracted.sam
  $SCRIPT_DIR/bin/contig_assembly $DATASET_DIR/$TUMOR_BAM.extracted.fastq $DATASET_DIR/$TUMOR_BAM.extracted.sam $DATASET_DIR/$TUMOR_BAM.extracted.contig.fastq
  bwa mem $REFERENCE_FASTA.index $DATASET_DIR/$TUMOR_BAM.extracted.contig.fastq > $DATASET_DIR/$TUMOR_BAM.extracted.contig.sam
  $SCRIPT_DIR/bin/summary_info $DATASET_DIR/$TUMOR_BAM.extracted.contig.fastq $DATASET_DIR/$TUMOR_BAM.extracted.contig.sam $DATASET_DIR/$TUMOR_BAM.result.bedpe
  echo "[pipeline] contig_assembly and summary_info Done."
  echo "[pipeline] The result is in $DATASET_DIR/$TUMOR_BAM.result.bedpe "
  date
fi

echo "[pipeline] Finished. "
date
