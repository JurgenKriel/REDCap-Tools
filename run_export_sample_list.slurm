#!/bin/bash
#SBATCH --job-name=redcap_sample_list
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/GLIMMER/scripts/logs/export_sample_list_%j.log
#SBATCH --error=/vast/projects/GLIMMER/scripts/logs/export_sample_list_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kriel.j@wehi.edu.au

echo "Job started: $(date)"
echo "Running on: $(hostname)"

source /vast/projects/GLIMMER/spatial_upload/venv/bin/activate

DATE=$(date +%Y-%m-%d)
OUTPUT="/vast/projects/GLIMMER/scripts/logs/sample_list_${DATE}.csv"

# Optional: add --prefixes GL GX JJ PRISM to broaden the search
python3 /vast/projects/GLIMMER/scripts/export_redcap_sample_list.py \
    --output "${OUTPUT}"

echo ""
echo "Job finished: $(date)"
echo "Output: ${OUTPUT}"
