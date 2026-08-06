#!/bin/bash
#SBATCH --job-name=lowpass_nf
#SBATCH --nodelist=node4
#SBATCH --cpus-per-task=77
#SBATCH --mem=305G
#SBATCH --time=20-00:00:00
#SBATCH --output=logs/nextflow_%j.out
#SBATCH --error=logs/nextflow_%j.err

set -euo pipefail

module purge
module load miniconda3
module load apptainer

source ~/.bashrc
conda activate nextflow

nextflow run main.nf \
    -profile singularity \
    -resume

conda deactivate