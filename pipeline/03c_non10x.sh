#!/bin/bash
# Step 4 for a non-10x droplet library whose R1 has the 10x v3 layout (16 nt barcode + 12 nt UMI) but barcodes that are
# not on the 10x list, for example BMKMANU DG1000 (check the read layout of other platforms first: R1 length 28 or more,
# barcode in nt 1-16, UMI in nt 17-28). Cell Ranger alone fails here ("extremely low rate of correct barcodes",
# 0.9% for the atlas mantle library).
# Two routes, both used for the atlas mantle library; run both and compare them with 05_compare_counters.R:
#   A (Cell Ranger): find the library's barcodes, give each a 10x barcode, rewrite R1, run Cell Ranger, rename back
#   C (STARsolo)   : the same barcode list as STARsolo whitelist, one-mismatch correction, Cell Ranger-like settings
# usage: bash 03c_non10x.sh <R1.fastq.gz> <R2.fastq.gz> <ref_dir (01_get_reference.sh ... star)> <out_dir> [threads]
set -euo pipefail
R1=$(readlink -f ${1:?R1}); R2=$(readlink -f ${2:?R2}); REF=$(readlink -f ${3:?ref_dir}); OUT=${4:?out_dir}; THR=${5:-16}
HERE=$(dirname "$(readlink -f "$0")"); mkdir -p "$OUT"; cd "$OUT"
SOLO="--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 --soloStrand Forward
      --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --soloFeatures Gene GeneFull
      --soloCellFilter EmptyDrops_CR --soloMultiMappers Unique --outSAMtype None --readFilesCommand zcat"
# 1. barcodes seen in the library: STARsolo without a whitelist (raw output = every barcode with a counted read)
mkdir -p seen; STAR --runThreadN "$THR" --genomeDir "$REF/star_index" --readFilesIn "$R2" "$R1" $SOLO --soloCBwhitelist None --outFileNamePrefix seen/
# 2. barcode map (barcodes with >= 10 UMIs; the rest cannot become cells)
python3 "$HERE/bmk_to_10x.py" map seen/Solo.out/GeneFull/raw "$REF/3M-february-2018.txt" map.tsv 10
# --- route A: Cell Ranger
mkdir -p fq10x; python3 "$HERE/bmk_to_10x.py" conv map.tsv "$R1" fq10x/M_S1_L001_R1_001.fastq.gz; ln -sf "$R2" fq10x/M_S1_L001_R2_001.fastq.gz
cellranger count --id=A_cellranger --transcriptome="$REF/xbMagGiga1_1_rs2024_06" --fastqs=fq10x --sample=M --chemistry=SC3Pv3 \
  --include-introns=true --no-bam --localcores="$THR" --disable-ui
F=A_cellranger/outs/filtered_feature_bc_matrix; mkdir -p A_filtered_original_barcodes
ln -sf "$PWD/$F/matrix.mtx.gz" "$PWD/$F/features.tsv.gz" A_filtered_original_barcodes/
python3 "$HERE/bmk_to_10x.py" back map.tsv $F/barcodes.tsv.gz A_filtered_original_barcodes/barcodes.tsv.gz
# --- route C: STARsolo with the library's own whitelist
mkdir -p C_starsolo; tail -n +2 map.tsv | cut -f1 > C_starsolo/whitelist.txt
STAR --runThreadN "$THR" --genomeDir "$REF/star_index" --readFilesIn "$R2" "$R1" $SOLO --soloCBwhitelist C_starsolo/whitelist.txt \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outFileNamePrefix C_starsolo/
# --- compare (same barcodes now): cells, depth, pseudobulk, clusters
Rscript "$HERE/05_compare_counters.R" A_filtered_original_barcodes C_starsolo/Solo.out/GeneFull/filtered compare CellRanger STARsolo
echo "use A_filtered_original_barcodes (default) or C_starsolo/Solo.out/GeneFull/filtered with 04_seurat.R"
