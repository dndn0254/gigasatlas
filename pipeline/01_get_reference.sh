#!/bin/bash
# Step 1-2: download the current Pacific oyster genome (NCBI RefSeq xbMagGiga1.1, annotation RS_2024_06) and build
# the Cell Ranger reference (default). Add 'star' to also build the STAR index for the STARsolo option (03b_starsolo.sh).
# usage: bash 01_get_reference.sh <out_dir> [threads] [star]
set -euo pipefail
OUT=${1:?out_dir}; THR=${2:-8}; STAR_TOO=${3:-}; mkdir -p "$OUT"; cd "$OUT"
U=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1
P=GCF_963853765.1_xbMagGiga1.1
wget -c "$U/${P}_genomic.fna.gz" "$U/${P}_genomic.gtf.gz" "$U/md5checksums.txt"
md5sum -c --ignore-missing md5checksums.txt
gunzip -kf ${P}_genomic.fna.gz ${P}_genomic.gtf.gz
# Cell Ranger uses exon lines only; keep exons that belong to a transcript (33,068 genes, the same genes as STARsolo).
# No gene-type filter (cellranger mkgtf) is applied, so non-coding genes stay in, as in the atlas.
awk -F'\t' '$0 ~ /^#/ || ($3 == "exon" && $9 !~ /transcript_id ""/)' ${P}_genomic.gtf > genes.exon.gtf
cellranger mkref --genome=xbMagGiga1_1_rs2024_06 --fasta=${P}_genomic.fna --genes=genes.exon.gtf --nthreads="$THR"
echo "Cell Ranger reference ready: $OUT/xbMagGiga1_1_rs2024_06"
if [ "$STAR_TOO" = star ]; then
  # small-genome rule from the STAR manual: genomeSAindexNbases = min(14, log2(genome length) / 2 - 1)
  L=$(grep -v '>' ${P}_genomic.fna | tr -d '\n' | wc -c)
  NB=$(awk -v l="$L" 'BEGIN{n=int(log(l)/log(2)/2-1); print (n<14?n:14)}')
  mkdir -p star_index
  STAR --runMode genomeGenerate --runThreadN "$THR" --genomeDir star_index --genomeFastaFiles ${P}_genomic.fna \
       --sjdbGTFfile ${P}_genomic.gtf --sjdbOverhang 100 --genomeSAindexNbases "$NB"
  # 10x barcode whitelists for STARsolo (v2: 737K-august-2016, v3/v3.1: 3M-february-2018)
  wget -c https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz
  wget -c https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz
  gunzip -kf 737K-august-2016.txt.gz 3M-february-2018.txt.gz
  echo "STAR index ready: $OUT/star_index"
fi
