# GLIMMER REDCap Tools

Python and SLURM scripts for managing spatial transcriptomics (ST) data submissions and queries against the GLIMMER REDCap project at WEHI.

---

## Repository Contents

| Script | Purpose |
|---|---|
| `lookup_gl.py` | Quick command-line TRU-ID lookup for one or more GL/GX/JJ IDs |
| `export_redcap_sample_list.py` | Export all samples + TRU-IDs across all DAGs to a CSV |
| `run_export_sample_list.slurm` | SLURM wrapper for `export_redcap_sample_list.py` |
| `report_instruments.py` | Report REDCap instrument completion / decoded field values for specific samples |
| `run_report_instruments.slurm` | SLURM wrapper for `report_instruments.py` |
| `generate_redcap_st_upload_v6.py` | Scan a Xenium ST batch directory and generate / upload REDCap records |
| `submit_redcap_upload_v6.sh` | SLURM wrapper for `generate_redcap_st_upload_v6.py` |
| `timepoint_map_template.csv` | Template for per-sample metadata (copy and populate before running the upload script) |

---

## Prerequisites

- Python 3.8+
- SLURM HPC cluster (Petrichor / WEHI HPC) for the `.slurm` / `.sh` wrappers
- REDCap API access token for the GLIMMER project

---

## Environment Setup

A shared virtual environment is expected at `/vast/projects/GLIMMER/spatial_upload/venv/`.
To create it from scratch:

```bash
python3 -m venv /vast/projects/GLIMMER/spatial_upload/venv
source /vast/projects/GLIMMER/spatial_upload/venv/bin/activate
pip install -r requirements.txt
```

### Python Dependencies

| Package | Purpose |
|---|---|
| `PyCap` | REDCap API client (`from redcap import Project`) |
| `pandas` | DataFrame / CSV handling in the ST upload script |

---

## Configuration

### API Token

Store your REDCap API token as plain text:

```
/vast/projects/GLIMMER/scripts/redcap_api_token.txt
```

All scripts read from this path by default; override with `--token-file`.
**Do not commit this file** — it is listed in `.gitignore`.

### Log / Output Directory

Default output directory for generated CSVs and SLURM logs:

```
/vast/projects/GLIMMER/scripts/logs/
```

Create it if it does not exist:

```bash
mkdir -p /vast/projects/GLIMMER/scripts/logs
```

---

## Script Reference

### `lookup_gl.py` — Quick TRU-ID Lookup

Looks up one or more GL/GX/JJ tissue bank IDs and prints their REDCap `record_id` (TRU-ID).

```bash
python3 lookup_gl.py GL0229
python3 lookup_gl.py GL0229 GL0201 GX0004
```

Sample output:

```
GL Number       BCRL ID (TRU-ID)     tissue_bank_id in REDCap
------------------------------------------------------------
GL0229          1243                 GL0229
GL0201          NOT FOUND            -
```

---

### `export_redcap_sample_list.py` — Export All Samples

Performs a bulk export across **all accessible DAGs** to produce a complete `tissue_bank_id → record_id` mapping as a CSV. Records are deduplicated by `record_id` and sorted alphabetically/numerically.

```bash
# All prefixes, output to logs/sample_list_YYYY-MM-DD.csv
python3 export_redcap_sample_list.py

# Custom output path
python3 export_redcap_sample_list.py --output /path/my_list.csv

# Restrict to specific prefixes
python3 export_redcap_sample_list.py --prefixes GL GX JJ PRISM
```

**Output columns:** `tissue_bank_id`, `record_id`

**Via SLURM:**
```bash
sbatch run_export_sample_list.slurm
```

The script prints a per-DAG breakdown at runtime, e.g.:

```
Fetching records from REDCap across all DAGs...
    default context           15 new  (running total: 15)
    bpop_all_samples         515 new  (running total: 530)
    bcrl_only                108 new  (running total: 638)
    ...
  638 total sample(s) found (all prefixes).
```

---

### `report_instruments.py` — Instrument Status Report

Reports REDCap instrument completion status and optionally prints decoded field values for a list of samples.

```bash
# Summary for built-in DEFAULT_SAMPLE_IDS list
python3 report_instruments.py

# Specific samples
python3 report_instruments.py --gl GL0283 GX0003 JJ14

# Detailed decoded field values (stdout)
python3 report_instruments.py --gl GL0283 --detail

# Save summary CSV + detail JSON + instances CSV
python3 report_instruments.py --gl GL0283 --detail --output /path/report.csv
```

**Via SLURM:**
```bash
sbatch run_report_instruments.slurm
```

Edit `DEFAULT_SAMPLE_IDS` at the top of the Python file, or pass IDs via `--gl`, to control which samples are queried.

---

### `generate_redcap_st_upload_v6.py` + `submit_redcap_upload_v6.sh` — ST Data Upload

Scans a Xenium ST batch directory, matches each output folder to a REDCap record via `tissue_bank_id`, and either:
- Generates upload-ready CSVs for manual review, or
- Pushes records directly to REDCap via the API.

Two CSV types are produced:
- `*_new.csv` — full field set for brand-new ST instances
- `*_backfill.csv` — minimal patch (time point + facility fields only) for instances that already exist

#### Per-Run Configuration

Edit the **PER-RUN CONFIGURATION** block at the top of `submit_redcap_upload_v6.sh` before each batch:

```bash
DATA_DIR="/stornext/Projects/GLIMMER/data_raw/ST/ST_XXXXXX"  # path to the batch
TIME_POINT_MAP="/vast/projects/GLIMMER/spatial_upload/timepoint_map.csv"
LAB_CONTACT="Jurgen Kriel"
CURRENT_INSTANCE="1"
FACILITY_CONTACT="2"    # radio value from REDCap codebook
FACILITY_PROTOCOL="1"   # radio value from REDCap codebook
PANEL="1"               # optional
MULTIPLEXING="0"        # 0=No, 1=Yes
SAMPLE_TYPE_SPT="TFF"   # TFF=Tumour FFPE, TFR=Tumour fresh/frozen
ANALYTE_SPT="TIS"       # TIS=Tissue, RNA=RNA
DO_UPLOAD="0"           # 0=dry-run (CSV only), 1=push to REDCap
```

#### Recommended Workflow

```bash
# 1. Dry-run — generates CSVs, does NOT upload
#    Set DO_UPLOAD="0" in the script, then:
sbatch submit_redcap_upload_v6.sh

# 2. Review the generated *_new.csv and *_backfill.csv

# 3. Upload — set DO_UPLOAD="1", then resubmit:
sbatch submit_redcap_upload_v6.sh
```

#### Supported Folder Naming Schemes

The parser automatically handles all Xenium output folder variants:

| Format | Example | Parsed As |
|---|---|---|
| GL/GX flat (original) | `output-XETG__SLIDE__GL0229-2__DATE__TIME` | `GL0229-2` |
| GL/GX hyphenated prefix | `output-XETG__SLIDE__GL-0198-1__DATE__TIME` | `GL0198-1` |
| GL/GX underscore | `output-XETG__SLIDE__GL_0175_1__DATE__TIME` | `GL0175-1` |
| Sample appended after time | `output-XETG__SLIDE__Region_3__DATE__TIME_GL_0290-1` | `GL0290-1` |
| Non-GL with embedded time point | `output-XETG__SLIDE__JJ73_1375_T1_Region_2__DATE__TIME` | `JJ73_1375`, T1 |
| Old Region style (parent dir) | `GL0004/output-XETG__SLIDE__Region_2__DATE__TIME` | `GL0004` |

---

## `timepoint_map.csv` Format

Required by `generate_redcap_st_upload_v6.py`. Copy `timepoint_map_template.csv` to `timepoint_map.csv` and populate it before running the upload script. **Do not commit the populated file** (it is in `.gitignore`).

| Column | Required | Description |
|---|---|---|
| `gl_id` | Yes | Sample ID as it appears in the folder name (e.g. `GL0229`, `GL0229-2`, `JJ73_1375_T1`) |
| `time_point` | Yes | REDCap time point code (`T0`, `T1`, `T2`, …) |
| `histology_block` | No | Histology block identifier |
| `redcap_id` | No | Override `tissue_bank_id` for REDCap lookup — use when the folder name differs from the REDCap registration (e.g. `JJ73_1375` is registered as `JJ73`) |

**Tip:** For samples with an embedded time point in the folder name (e.g. `JJ73_1375_T1_Region_2`), create one row per time point keyed as `<sample>_T<N>` so each time point gets its correct histology block:

```csv
JJ73_1375,T0,HB-T0,JJ73
JJ73_1375_T1,T1,HB-T1,JJ73
JJ73_1375_T2,T2,HB-T2,JJ73
```

---

## Notes

- The `glimmer` DAG switch may return a permissions error — this is expected; those records are already returned by the default export context.
- Scripts with DAG switching (`export_redcap_sample_list.py`, `report_instruments.py`) gracefully skip DAGs they cannot access.
- All SLURM scripts source the venv at `/vast/projects/GLIMMER/spatial_upload/venv/bin/activate` and send mail to `kriel.j@wehi.edu.au` on job end/failure.
