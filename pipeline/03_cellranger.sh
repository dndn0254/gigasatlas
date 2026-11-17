#!/bin/bash
# Step 4 (default): count one 10x library with Cell Ranger, exactly as for every 10x library of the Gigas Cell Atlas.
# usage: bash 03_cellranger.sh <fastq_dir> <sample> <ref_dir from step 1> <out_dir> <cells|nuclei> [threads] [mem_gb]
#   fastq_dir must hold <sample>_S1_L001_R1_001.fastq.gz (barcode + UMI) and <sample>_S1_L001_R2_001.fastq.gz (cDNA);
#   more lanes or runs of the same library: L002, L003 ... (Cell Ranger adds them up).
#   SRA downloads (02_get_reads.sh) have other names: link them, do not copy, for example
#     ln -s $PWD/fastq/SRR22554543_1.fastq.gz fq/E1_S1_L001_R1_001.fastq.gz
#     ln -s $PWD/fastq/SRR22554543_2.fastq.gz fq/E1_S1_L001_R2_001.fastq.gz
#   cells  : whole cells (exons only, --include-introns=false), as for the embryo and larva libraries of the atlas
#   nuclei : single-nucleus or tissue cells with many unspliced reads (exons + introns), as for hemocyte, gonad, ganglia
# output: <out_dir>/<sample>/outs/filtered_feature_bc_matrix  (read it with 04_seurat.R)
# atlas version: Cell Ranger 7.0.1 (newer versions run with the same options but were not tested by us).
# From 7.0 the default is --include-introns=true, so always set it explicitly.
set -euo pipefail
FQ=${1:?fastq_dir}; S=${2:?sample}; REF=${3:?ref_dir}; OUT=${4:?out_dir}; KIND=${5:?cells_or_nuclei}; THR=${6:-16}; MEM=${7:-54}
case "$KIND" in cells) INT=false ;; nuclei) INT=true ;; *) echo "kind must be cells or nuclei"; exit 1 ;; esac
mkdir -p "$OUT"; cd "$OUT"
cellranger count --id="$S" --transcriptome="$REF/xbMagGiga1_1_rs2024_06" --fastqs="$FQ" --sample="$S" \
  --chemistry=auto --include-introns="$INT" --no-bam --localcores="$THR" --localmem="$MEM" --disable-ui
grep -o '"Estimated Number of Cells","[^"]*"' "$S/outs/metrics_summary.csv" || head -2 "$S/outs/metrics_summary.csv"
# chemistry: 'auto' lets Cell Ranger detect 10x v2 / v3 / v3.1 (the atlas checked each call against the read length).
# If auto fails on a short or trimmed barcode read, set it: --chemistry=SC3Pv2 (26 nt R1) or SC3Pv3 (28 nt R1).
