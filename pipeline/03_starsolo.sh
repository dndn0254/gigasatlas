#!/bin/bash
# Step 4: align and count one 10x library with STARsolo, using the Cell Ranger-matched settings of the Gigas Cell Atlas.
# usage: bash 03_starsolo.sh <R1.fastq.gz> <R2.fastq.gz> <ref_dir from step 1> <out_dir> <v2|v3> [threads]
#   R1 = barcode + UMI, R2 = cDNA. v2: 16 nt barcode + 10 nt UMI; v3 / v3.1: 16 nt barcode + 12 nt UMI.
# output: <out_dir>/Solo.out/Gene/filtered (whole cells) and Solo.out/GeneFull/filtered (nuclei: exons + introns)
set -euo pipefail
R1=${1:?R1}; R2=${2:?R2}; REF=${3:?ref_dir}; OUT=${4:?out_dir}; CHEM=${5:?v2_or_v3}; THR=${6:-8}
if [ "$CHEM" = v2 ]; then UMI=10; WL=$REF/737K-august-2016.txt; else UMI=12; WL=$REF/3M-february-2018.txt; fi
mkdir -p "$OUT"
STAR --runThreadN "$THR" --genomeDir "$REF/star_index" --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
     --outFileNamePrefix "$OUT/" --outSAMtype None \
     --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen "$UMI" --soloBarcodeReadLength 0 \
     --soloCBwhitelist "$WL" --soloStrand Forward \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
     --clipAdapterType CellRanger4 --soloFeatures Gene GeneFull --soloCellFilter EmptyDrops_CR --soloMultiMappers Unique
grep "Estimated Number of Cells" "$OUT/Solo.out/Gene/Summary.csv"
# Non-10x chemistry (for example BMK): set --soloCBlen / --soloUMIstart / --soloUMIlen to the read layout and use
# --soloCBwhitelist None (the atlas did this for the BMK mantle library: barcode 1-16, UMI 17-28).
