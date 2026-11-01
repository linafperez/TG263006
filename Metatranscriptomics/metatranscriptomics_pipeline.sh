#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-}"

usage() {
    echo "Usage:"
    echo "  ./metatranscriptomics_pipeline.sh stub"
    echo "  ./metatranscriptomics_pipeline.sh real"
    echo "  ./metatranscriptomics_pipeline.sh resume"
}

if [[ -z "${MODE}" ]]; then
    usage
    exit 1
fi

cd "${PROJECT_DIR}"

case "${MODE}" in

    stub)
        echo "Submitting isolated stub validation..."
        mkdir -p "${PROJECT_DIR}/stub"
        sbatch "${PROJECT_DIR}/run_stub_pipeline.sbatch"
        ;;

    real)
        if [[ -e "${PROJECT_DIR}/results/95_samples" ]]; then
            echo "ERROR: results/95_samples already exists."
            echo "Use 'resume' if this belongs to an interrupted production run."
            exit 1
        fi

        echo "Submitting fresh 95-sample production run..."

        sbatch \
            --job-name=mt_real_95 \
            --partition=medium \
            --time=5-00:00:00 \
            --cpus-per-task=8 \
            --mem=16G \
            --output=metatranscriptomics_real_%j.out \
            --error=metatranscriptomics_real_%j.err \
            --wrap="bash -lc 'cd \"${PROJECT_DIR}\" && module load singularity/3.7.1 && source /hpcfs/apps/conda4.12.0/etc/profile.d/conda.sh && conda activate nextflow && export PATH=/hpcfs/apps/conda4.12.0/bin:\$PATH && export NXF_OPTS=\"-Xms256m -Xmx4g\" && nextflow run main.nf -profile slurm,singularity -work-dir work --outdir results/95_samples -ansi-log false'"
        ;;

    resume)
        echo "Submitting resumed 95-sample production run..."

        sbatch \
            --job-name=mt_real_95 \
            --partition=medium \
            --time=5-00:00:00 \
            --cpus-per-task=8 \
            --mem=16G \
            --output=metatranscriptomics_real_%j.out \
            --error=metatranscriptomics_real_%j.err \
            --wrap="bash -lc 'cd \"${PROJECT_DIR}\" && module load singularity/3.7.1 && source /hpcfs/apps/conda4.12.0/etc/profile.d/conda.sh && conda activate nextflow && export PATH=/hpcfs/apps/conda4.12.0/bin:\$PATH && export NXF_OPTS=\"-Xms256m -Xmx4g\" && nextflow run main.nf -profile slurm,singularity -work-dir work --outdir results/95_samples -resume -ansi-log false'"
        ;;

    *)
        echo "ERROR: unknown mode '${MODE}'"
        usage
        exit 1
        ;;
esac
