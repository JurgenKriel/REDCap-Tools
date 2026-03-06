#!/bin/bash
#SBATCH --job-name=redcap_sample_list
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=/path/to/logs/export_sample_list_%j.log
#SBATCH --error=/path/to/logs/export_sample_list_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

echo "Job started: $(date)"
echo "Running on: $(hostname)"

source /path/to/venv/bin/activate

DATE=$(date +%Y-%m-%d)
OUTPUT="/path/to/logs/sample_list_${DATE}.csv"

# Optional: add --prefixes GL GX JJ PRISM to broaden the search
python3 /path/to/scripts/export_redcap_sample_list.py \
    --output "${OUTPUT}"

echo ""
echo "Job finished: $(date)"
echo "Output: ${OUTPUT}"
