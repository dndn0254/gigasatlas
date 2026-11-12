#!/bin/bash
# Step 3: download one public run from SRA (skip this if you have your own FASTQ).
# usage: bash 02_get_reads.sh <SRR accession> <out_dir> [threads]
# needs SRA Toolkit (prefetch, fasterq-dump). 10x runs give two reads: R1 = cell barcode + UMI, R2 = cDNA.
set -euo pipefail
ACC=${1:?accession}; OUT=${2:?out_dir}; THR=${3:-8}; mkdir -p "$OUT"; cd "$OUT"
prefetch --max-size 200G "$ACC"
fasterq-dump --split-files --include-technical --threads "$THR" "$ACC"
for f in ${ACC}_*.fastq; do gzip "$f"; done
ls -l ${ACC}_*.fastq.gz
# check which file is which: the barcode read is 26 nt (v2) or 28 nt (v3); rename to ${ACC}_R1.fastq.gz and ${ACC}_R2.fastq.gz
for f in ${ACC}_*.fastq.gz; do echo "$f read length: $(zcat "$f" | head -2 | tail -1 | tr -d '\n' | wc -c)"; done
