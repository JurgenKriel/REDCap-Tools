#!/bin/bash
#SBATCH --job-name=redcap_report
#SBATCH --output=/path/to/logs/report_instruments_%j.log
#SBATCH --error=/path/to/logs/report_instruments_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=regular
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

echo "Job started: $(date)"
echo "Running on: $(hostname)"

source /path/to/venv/bin/activate

# ── Output paths ──────────────────────────────────────────────────────────────
DATE=$(date +%Y-%m-%d)
OUT_DIR="/path/to/logs"
OUTPUT="${OUT_DIR}/instrument_report_${DATE}.csv"
# With --detail, also writes:
#   instrument_report_<DATE>_detail.json
#   instrument_report_<DATE>_instances.csv

python3 /path/to/scripts/report_instruments.py \
    #--detail 
    \ --output "${OUTPUT}"

echo ""
echo "Job finished: $(date)"
echo "Outputs:"
echo "  Summary  : ${OUTPUT}"
echo "  Detail   : ${OUT_DIR}/instrument_report_${DATE}_detail.json"
echo "  Instances: ${OUT_DIR}/instrument_report_${DATE}_instances.csv"
