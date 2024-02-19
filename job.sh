#!/bin/bash
#SBATCH --job-name=als_job
#SBATCH --output=als_job2.out
#SBATCH --error=als_job2.out
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --partition=batch
python -u STAGE_2/2.3_gsea_per_celltype.py