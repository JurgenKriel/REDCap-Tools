#!/bin/bash
#SBATCH --job-name=redcap_upload_v6
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/vast/projects/GLIMMER/scripts/logs/redcap_upload_v6_%j.log
#SBATCH --error=/vast/projects/GLIMMER/scripts/logs/redcap_upload_v6_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kriel.j@wehi.edu.au
#SBATCH --chdir=/vast/projects/GLIMMER/spatial_upload/

# ============================================================
# PER-RUN CONFIGURATION — edit these for each submission batch
# ============================================================
TIME_POINT_MAP="/vast/projects/GLIMMER/spatial_upload/timepoint_map.csv"
LAB_CONTACT="Jurgen Kriel"
DATA_DIR="/stornext/Projects/GLIMMER/data_raw/ST/ST_251212"
CURRENT_INSTANCE="1"
FACILITY_CONTACT="2"                # Radio value from REDCap codebook
FACILITY_PROTOCOL="1"               # Radio value from REDCap codebook
PANEL="1"                            # Optional: set if needed, e.g. "1"
MULTIPLEXING="0"                    # 0=No, 1=Yes
#SAMPLE_TYPE_SPT="TFF"               # TFF=Tumour in FFPE block, TFR=Tumour fresh/frozen
#ANALYTE_SPT="TIS"                   # TIS=Tissue, RNA=RNA
DO_UPLOAD="1"                       # Set to 1 to upload to REDCap, 0 for dry-run
# ============================================================

# Read API token from file
TOKEN_FILE="/vast/projects/GLIMMER/scripts/redcap_api_token.txt"
if [ ! -f "$TOKEN_FILE" ]; then
    echo "Error: API token file not found: $TOKEN_FILE"
    exit 1
fi
REDCAP_API_TOKEN=$(cat "$TOKEN_FILE")
if [ -z "$REDCAP_API_TOKEN" ]; then
    echo "Error: API token file is empty: $TOKEN_FILE"
    exit 1
fi

# Validate time point map exists
if [ ! -f "$TIME_POINT_MAP" ]; then
    echo "Error: Time point map not found: $TIME_POINT_MAP"
    echo "Expected format (CSV):"
    echo "  gl_id,time_point,histology_block"
    echo "  GL0229,T0,HB-001"
    echo "  GL0229-2,T1,"
    exit 1
fi

# Load virtual environment
source /vast/projects/GLIMMER/spatial_upload/venv/bin/activate

# Print run configuration
echo "============================================================"
echo " REDCap ST Upload — Run Configuration (v6)"
echo "============================================================"
echo "  Time point map:    $TIME_POINT_MAP"
echo "  Lab contact:       $LAB_CONTACT"
echo "  Data directory:    $DATA_DIR"
echo "  Current instance:  $CURRENT_INSTANCE"
echo "  Facility contact:  $FACILITY_CONTACT"
echo "  Facility protocol: $FACILITY_PROTOCOL"
echo "  Panel:             ${PANEL:-not set}"
echo "  Multiplexing:      $MULTIPLEXING"
echo "  (Submission date derived from output folder name)"
echo "  (Data received date derived from output folder mtime)"
echo "============================================================"

# Construct command args
ARGS="--token $REDCAP_API_TOKEN"
ARGS="$ARGS --time-point-map $TIME_POINT_MAP"
ARGS="$ARGS --lab-contact \"$LAB_CONTACT\""
ARGS="$ARGS --data-dir $DATA_DIR"
ARGS="$ARGS --current-instance $CURRENT_INSTANCE"
ARGS="$ARGS --facility-contact $FACILITY_CONTACT"
ARGS="$ARGS --facility-protocol $FACILITY_PROTOCOL"
ARGS="$ARGS --multiplexing $MULTIPLEXING"
ARGS="$ARGS --sample-type-spt $SAMPLE_TYPE_SPT"
ARGS="$ARGS --analyte-spt $ANALYTE_SPT"

# Optional args — only added if variable is set and non-empty
[ -n "${PANEL:-}" ] && ARGS="$ARGS --panel \"$PANEL\""

if [ "${DO_UPLOAD}" == "1" ]; then
    ARGS="$ARGS --upload"
    echo "Upload mode ENABLED — records will be pushed to REDCap."
    echo "(Existing instances: time_point_spt, current_instance_spt, facility_st + sample_type_analysed_st___3 will be patched — all other fields left untouched.)"
else
    echo "Dry-run mode — CSV will be generated but NOT uploaded."
    echo "To upload, set DO_UPLOAD=\"1\" in the PER-RUN CONFIGURATION block."
fi

echo ""
echo "Starting REDCap upload preparation (v6)..."
eval python3 /vast/projects/GLIMMER/scripts/generate_redcap_st_upload_v6.py $ARGS
echo "Script finished."