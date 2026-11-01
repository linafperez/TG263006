#!/bin/bash
#SBATCH --job-name=STAR_INDEX
#SBATCH --output=STAR_INDEX_%j.out
#SBATCH --error=STAR_INDEX_%j.err
#SBATCH --partition=short
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=04:00:00

set -euo pipefail

module load singularity/3.7.1

BASE="$HOME/Metatranscriptomics/data/references/human_GRCh38_p14_GENCODE_v44"

singularity exec \
  docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_7 \
  STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir "$BASE/STAR_IDX" \
    --genomeFastaFiles "$BASE/GRCh38.p14.genome.fa" \
    --sjdbGTFfile "$BASE/gencode.v44.annotation.gtf" \
    --sjdbOverhang 149
