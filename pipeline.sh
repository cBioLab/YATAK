#!/bin/sh

if [ "$#" != "6" ] && [ "$#" != "7" ]; then
  echo "[!] Usage : ./pipeline.sh <reference-directory> <reference.fa> <dataset-directory> <normal.bam> <tumor.bam> <number of threads> (<clear-map/clear-index>)"
  exit
fi

if [ "$#" = "7" ]; then
  if [ "$7" = "clear-map" ]; then
    rm $3/$5.extracted.fastq
    rm $3/$5.summary.txt
  elif [ "$7" = "clear-index" ]; then
    rm $3/$2.$4.index.cst
    rm $3/$2.$4.index.refbv
    rm $3/$5.extracted.fastq
    rm $3/$5.summary.txt
  fi
fi

echo "!! Note that samtools/bwa/mafft needs to be installed. !!"
date

### indexing for bwa
if [ ! -e "$1/$2.index.amb" ] || [ ! -e "$1/$2.index.ann" ] || [ ! -e "$1/$2.index.bwt" ] || [ ! -e "$1/$2.index.pac" ] || [ ! -e "$1/$2.index.sa" ]; then
  if [ ! -e "$1/$2" ]; then
    echo "[!] $1/$2 not found. "
    exit
  fi
  echo "[pipeline] $1/$2.index is not found. Indexing process started ... "
  bwa index -p $1/$2.index $1/$2
  echo "[pipeline] bwa index Done."
  date
fi

if [ ! -e "$1/$2.fixed.fa" ]; then
  if [ ! -e "$1/$2" ]; then
    echo "[!] $1/$2 not found. "
    exit
  fi
  echo "[pipeline] $1/$2.fixed.fa is not found. Fixing process started ... "
  ./bin/fix_reference_for_indexing $1/$2 $1/$2.fixed.fa
  echo "[pipeline] fix_reference_for_indexing Done."
  date
fi

### filter read from SAM
if [ ! -e "$3/$4.filtered.fastq" ]; then
  if [ ! -e "$3/$4" ]; then
    echo "[!] $3/$4 not found. "
    exit
  fi
  echo "[pipeline] $3/$4.filtered.fastq is not found. Filtering process started ... "
  ./scripts/filter_read_from_sam.py $3/$4 $3/$4.filtered.bam
  samtools fastq $3/$4.filtered.bam > $3/$4.filtered.fastq
  echo "[pipeline] filter_read_from_sam (normal) Done."
  date
fi

if [ ! -e "$3/$5.filtered.fastq" ]; then
  if [ ! -e "$3/$5" ]; then
    echo "[!] $3/$5 not found. "
    exit
  fi
  echo "[pipeline] $3/$5.filtered.fastq is not found. Filtering process started ... "
  ./scripts/filter_read_from_sam.py $3/$5 $3/$5.filtered.bam
  samtools fastq $3/$5.filtered.bam > $3/$5.filtered.fastq
  echo "[pipeline] filter_read_from_sam (tumor) Done."
  date
fi

### indexing
if [ ! -e "$3/$2.$4.index.cst" ] || [ ! -e "$3/$2.$4.index.refbv" ]; then
  echo "[pipeline] $3/$2.$4.index is not found. Indexing process started ... "
  ./bin/construct_index $1/$2.fixed.fa $3/$4.filtered.fastq $3/$2.$4.index $6
  rm $3/$4.filtered.fastq_reference_read_text.tmp
  echo "[pipeline] construct_index Done."
  date
fi

### alignment on graph
if [ ! -e "$3/$5.extracted.fastq" ]; then
  echo "[pipeline] Main filtering process started ... "
  ./bin/rf_mapping $3/$5.filtered.fastq $3/$2.$4.index $3/$5.extracted.fastq $6
  echo "[pipeline] pruning_read Done."
  date
fi

### report on reference
if [ ! -e "$3/$5.result.bedpe" ]; then
  echo "[pipeline] Creating report ... "
  bwa mem $1/$2.index $3/$5.extracted.fastq > $3/$5.extracted.sam
  ./bin/contig_assembly $3/$5.extracted.fastq $3/$5.extracted.sam $3/$5.extracted.contig.fastq
  bwa mem $1/$2.index $3/$5.extracted.contig.fastq > $3/$5.extracted.contig.sam
  ./bin/summary_info $3/$5.extracted.contig.fastq $3/$5.extracted.contig.sam $3/$5.result.bedpe
  echo "[pipeline] contig_assembly and summary_info Done."
  date
fi

echo "[pipeline] Finished. "
date
