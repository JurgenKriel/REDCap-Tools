import os
import re
import pandas as pd
import argparse
import urllib3
from datetime import datetime
from redcap import Project, RedcapError

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Configuration
DEFAULT_DATA_DIR = "/path/to/ST_batch_directory/"
OUTPUT_FILE      = "st_upload_{date}_v6_api{suffix}.csv"
DEFAULT_API_URL  = "https://your-redcap-instance.org/api/"

# ── REDCap @DEFAULT values (codebook fields 848-857) ──────────────────────────
REDCAP_DEFAULTS = {
    "sample_type_analysed_st___3": "1",    # Tumour checkbox
    "facility_st":                 "3",    # SCORE
}

# Sample-type SPT codes  (sample_type_analysed_spt_2)
SAMPLE_TYPE_SPT_FFPE   = "TFF"   # Tumour in FFPE block
SAMPLE_TYPE_SPT_FROZEN = "TFR"   # Tumour fresh or frozen

# Analyte SPT codes  (analyte_spt_2)
ANALYTE_SPT_TISSUE = "TIS"   # Tissue (FFPE block runs)
ANALYTE_SPT_RNA    = "RNA"   # RNA    (fresh/frozen runs)

# ── Directory name patterns ────────────────────────────────────────────────────
# Old style: output dir uses Region_X; GL/GX ID is the parent folder name
# e.g. GX0004/output-XETG00068__0043198__Region_2__20250211__053448
DIR_PATTERN_OLD = re.compile(
    r"output-[^_]+__[^_]+__Region_\d+__(?P<date>\d{8})__"
)

# Batch subdirectory: named with a leading date stamp
# e.g. 20250922__040656__20250922_G543_Batch7B_XEN
BATCH_SUBDIR_RE = re.compile(r"^\d{8}__")

# General output-folder extractor — works for all naming schemes.
# Captures the sample-field (everything between 2nd and 3rd '__') and the date.
OUTPUT_FOLDER_RE = re.compile(
    r"output-[^_]+__[^_]+__(?P<sf>.+)__(?P<date>\d{8})__(?P<time>\d+)(?:_(?P<trail>.+))?"
)

# Tissue suffixes that indicate non-tumour (normal/peri) tissue
NORMAL_SUFFIXES = {"Peri", "Cortex", "Non_enh"}


def load_time_point_map(filepath):
    """
    Load per-sample metadata from a CSV with columns:
      gl_id, time_point, histology_block (optional)
    Returns {gl_id: {"time_point": ..., "histology_block": ...}}
    """
    if not filepath:
        return {}
    if not os.path.exists(filepath):
        print(f"WARNING: Time point map file not found: {filepath}")
        return {}
    df = pd.read_csv(filepath, dtype=str).fillna("")
    mapping = {}
    for _, row in df.iterrows():
        gl_id = row["gl_id"].strip()
        tp    = row["time_point"].strip()
        if gl_id:
            mapping[gl_id] = {
                "time_point":      tp,
                "histology_block": row["histology_block"].strip() if "histology_block" in df.columns else "",
                # Optional redcap_id: alternative tissue_bank_id to use for REDCap
                # lookup when the sample is registered under a different ID
                # (e.g. JJ73_1375 → GX0008).
                "redcap_id":       row["redcap_id"].strip() if "redcap_id" in df.columns else "",
            }
    print(f"Loaded {len(mapping)} sample mappings from {filepath}")
    return mapping


def resolve_sample_info(gl_id, time_point_map):
    """
    Returns the metadata dict for a GL/GX ID.
    Tries full ID first, then falls back to base ID (before first '-').
    """
    if gl_id in time_point_map:
        return time_point_map[gl_id]
    base_id = gl_id.split("-")[0]
    if base_id in time_point_map:
        return time_point_map[base_id]
    return None


def parse_date_ymd(date_str):
    """Convert YYYYMMDD → YYYY-MM-DD for REDCap API."""
    try:
        return datetime.strptime(date_str, "%Y%m%d").strftime("%Y-%m-%d")
    except Exception:
        return ""


def reformat_date_to_ymd(date_str):
    """
    Accept DD-MM-YYYY (CLI input) or YYYY-MM-DD and normalise to YYYY-MM-DD
    for the REDCap API.
    """
    if not date_str:
        return ""
    for fmt in ("%d-%m-%Y", "%Y-%m-%d"):
        try:
            return datetime.strptime(date_str, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    print(f"WARNING: Could not parse date '{date_str}', leaving blank.")
    return ""


def get_redcap_project(api_url, api_token):
    try:
        try:
            project = Project(api_url, api_token, verify_ssl=False)
        except TypeError:
            project = Project(api_url, api_token)
        print(f"Connected to REDCap: {project.metadata[0]['field_label'] if project.metadata else 'Unknown'}")
        return project
    except Exception as e:
        print(f"Error connecting to REDCap: {e}")
        return None


def get_redcap_mapping(project):
    """Returns {tissue_bank_id: record_id}."""
    try:
        records = project.export_records(fields=["record_id", "tissue_bank_id"])
        mapping = {
            r["tissue_bank_id"]: r["record_id"]
            for r in records if r.get("tissue_bank_id")
        }
        print(f"Fetched {len(records)} records; {len(mapping)} with tissue_bank_id.")
        return mapping
    except RedcapError as e:
        print(f"REDCap API Error: {e}")
        return None


def get_existing_st_instances(project):
    """
    Fetch existing st_spatial_transcriptomics instances from REDCap.

    Returns a tuple of:
      existing_by_sample : {(record_id, sample_id_st): [instance_int, ...]}
          All known instance numbers per (record, sample_id) pair.
          Multiple values indicate duplicate entries that need manual review.
      max_inst_per_record: {record_id: max_instance_int}
          The highest instance number seen for each record, used to
          compute the next safe integer when creating new instances.

    BUG FIX (v6.1):
      Previously used a plain dict keyed by (rec_id, sample_id), so when a
      patient had multiple ST instances with the same sample_id the earlier
      instance numbers were silently overwritten. Now stores a list per key
      so duplicates are detected and flagged rather than lost.
    """
    try:
        records = project.export_records(
            forms=["st_spatial_transcriptomics"],
            fields=["record_id", "sample_id_st"],
            export_survey_fields=False,
            export_data_access_groups=False,
        )
        existing_by_sample  = {}
        max_inst_per_record = {}

        for r in records:
            if r.get("redcap_repeat_instrument") != "st_spatial_transcriptomics":
                continue
            rec_id    = str(r.get("record_id", "")).strip()
            sample_id = r.get("sample_id_st", "").strip()
            try:
                inst_num = int(r.get("redcap_repeat_instance", ""))
            except (ValueError, TypeError):
                continue

            # Track highest instance per record for next-integer computation
            max_inst_per_record[rec_id] = max(
                max_inst_per_record.get(rec_id, 0), inst_num
            )

            # Track all instances per (record, sample_id) pair
            if rec_id and sample_id:
                key = (rec_id, sample_id)
                existing_by_sample.setdefault(key, []).append(inst_num)

        # Warn about any ambiguous entries that will need manual resolution
        ambiguous = {k: v for k, v in existing_by_sample.items() if len(v) > 1}
        if ambiguous:
            print(f"\n  WARNING: {len(ambiguous)} (record, sample_id) pair(s) have multiple "
                  f"existing instances — these will be SKIPPED and need manual review:")
            for (rec_id, sid), insts in ambiguous.items():
                print(f"    record_id={rec_id}, sample_id_st={sid} → instances {insts}")

        n_entries = sum(len(v) for v in existing_by_sample.values())
        print(f"Found {n_entries} existing ST instance(s) across "
              f"{len(existing_by_sample)} unique (record, sample_id) pairs.")
        return existing_by_sample, max_inst_per_record

    except RedcapError as e:
        print(f"WARNING: Could not fetch existing ST instances: {e}")
        return {}, {}


def resolve_tru_id(gl_id, mapping):
    """Direct lookup, then fall back to base ID (strip suffix after first '-')."""
    if gl_id in mapping:
        return mapping[gl_id]
    base_id = gl_id.split("-")[0]
    if base_id in mapping:
        return mapping[base_id]
    return None


def get_tissue_suffix(gl_id):
    """
    Extract tissue-type suffix from a GL/GX ID.
    Returns empty string for numeric suffixes (timepoint indicators),
    e.g. GX0002-1 → '', GL0212-Tumour → 'Tumour'.
    """
    parts = gl_id.split("-", 1)
    if len(parts) > 1 and not parts[1].isdigit():
        return parts[1]
    return ""


def parse_output_folder(entry):
    """
    Parse a Xenium output folder name.

    Returns (sample_id, date_str, time_point_override) where:
      sample_id          — normalised sample identifier (None if unrecognised)
      date_str           — YYYYMMDD string from folder name
      time_point_override — time point extracted from folder name (e.g. 'T0'),
                           or None when the time point must come from the map

    Recognised naming schemes
    ─────────────────────────
    Original GL/GX (flat):
      output-XETG__SLIDE__GL0229-2__DATE__TIME
      output-XETG__SLIDE__GL0212-Tumour__DATE__TIME

    GL/GX with underscores (new):
      output-XETG__SLIDE__GL_0175_Region_2__DATE__TIME      → GL0175
      output-XETG__SLIDE__GL_0175_1_Region_4__DATE__TIME    → GL0175-1
      output-XETG__SLIDE__GL_0175_Cortex_Region_2__DATE__TIME → GL0175-Cortex

    Non-GL samples with embedded time point (JJ73, etc.):
      output-XETG__SLIDE__JJ73_1375_T3__DATE__TIME          → JJ73_1375, T3
      output-XETG__SLIDE__JJ73_1375_T0_Region_3__DATE__TIME → JJ73_1375, T0

    PRISM / other named samples:
      output-XETG__SLIDE__PRISM_11_Region_1__DATE__TIME     → PRISM_11

    Old-style Region (GL/GX ID is the parent directory):
      output-XETG__SLIDE__Region_2__DATE__TIME              → returns "REGION_ONLY"
    """
    if not entry.startswith("output-"):
        return None, None, None

    m = OUTPUT_FOLDER_RE.match(entry)
    if not m:
        return None, None, None

    sample_field = m.group("sf")
    date_str     = m.group("date")

    # ── Old-style Region: sample ID lives in the parent folder name ───────────
    # New batches append the sample ID AFTER the time field:
    #   output-XETG__SLIDE__Region_N__DATE__TIME_GL_0290-1
    # In that case `trail` holds the sample ID; fall back to REGION_ONLY otherwise.
    trail = m.group("trail")
    if re.match(r"Region_\d+$", sample_field):
        if trail:
            sample_field = trail   # e.g. "GL_0290-1", "Biopsy_25.2395.1"
        else:
            return "REGION_ONLY", date_str, None

    # ── Strip trailing _Region_N ──────────────────────────────────────────────
    field = re.sub(r"_Region_\d+$", "", sample_field)

    # ── Strip trailing embedded time point _TN ────────────────────────────────
    tp_match = re.search(r"_(?P<tp>T\d+)$", field)
    time_point_override = None
    if tp_match:
        time_point_override = tp_match.group("tp")
        field = field[:tp_match.start()]

    # ── New GL/GX underscore format: GL_0175, GL_0175_1, GL_0175_Cortex ───────
    m2 = re.match(r"^(?P<pfx>G[LX])_(?P<num>\d+)(?:[_-](?P<suf>.+))?$", field)
    if m2:
        sample_id = f"{m2.group('pfx')}{m2.group('num')}"
        if m2.group("suf"):
            sample_id = f"{sample_id}-{m2.group('suf')}"
        return sample_id, date_str, time_point_override

    # ── GL/GX hyphen format: GL-0198, GL-0198-1, GL-0198-Cortex ────────────────
    # Folder names sometimes include a hyphen between the prefix and digits
    # (e.g. GL-0198-1).  REDCap tissue_bank_id uses no hyphen (GL0198-1), so
    # normalise here so that all downstream lookups succeed.
    m3 = re.match(r"^(?P<pfx>G[LX])-(?P<num>\d+)(?:-(?P<suf>.+))?$", field)
    if m3:
        sample_id = f"{m3.group('pfx')}{m3.group('num')}"
        if m3.group("suf"):
            sample_id = f"{sample_id}-{m3.group('suf')}"
        return sample_id, date_str, time_point_override

    # ── Original GL/GX format: GL0229-2, GL0175-Tumour, GX0002 ───────────────
    if re.match(r"G[LX]\d+", field):
        return field, date_str, time_point_override

    # ── Any other named sample: PRISM_11, JJ73_1375, etc. ────────────────────
    if field and not field.startswith("Region"):
        return field, date_str, time_point_override

    return None, None, None


def scan_data_dir(data_dir):
    """
    Yields (sample_id, date_str, full_output_path, time_point_override) for all
    Xenium output dirs found.

    Handles three layouts:
      Layout 1 (flat):        DATA_DIR/output-XETG...__GL0229-2__DATE__...
      Layout 2 (GL folder):   DATA_DIR/GX0004/output-XETG...__Region_2__DATE__...
      Layout 3 (batch subdir):DATA_DIR/YYYYMMDD__BATCH/output-XETG...__GL_0175__DATE__...
    """
    for top_entry in sorted(os.listdir(data_dir)):
        top_path = os.path.join(data_dir, top_entry)
        if not os.path.isdir(top_path):
            continue

        # ── Layout 1: output folder directly in DATA_DIR ──────────────────────
        if top_entry.startswith("output-"):
            sample_id, date_str, tp_override = parse_output_folder(top_entry)
            if sample_id and sample_id != "REGION_ONLY":
                yield sample_id, date_str, top_path, tp_override
            else:
                print(f"  Skipping (no pattern match): {top_entry}")
            continue

        # ── Layout 3: batch subdir (YYYYMMDD__...) → descend one level ────────
        if BATCH_SUBDIR_RE.match(top_entry):
            for sub_entry in sorted(os.listdir(top_path)):
                sub_path = os.path.join(top_path, sub_entry)
                if not os.path.isdir(sub_path) or not sub_entry.startswith("output-"):
                    continue
                sample_id, date_str, tp_override = parse_output_folder(sub_entry)
                if sample_id and sample_id != "REGION_ONLY":
                    yield sample_id, date_str, sub_path, tp_override
                else:
                    print(f"  Skipping (no pattern match): {top_entry}/{sub_entry}")
            continue

        # ── Layout 2: GL/GX parent folder → Region_N subdirs ─────────────────
        if re.match(r"G[LX]\d+", top_entry):
            gl_id = top_entry
            for sub_entry in sorted(os.listdir(top_path)):
                sub_path = os.path.join(top_path, sub_entry)
                if not os.path.isdir(sub_path):
                    continue
                sub_match = DIR_PATTERN_OLD.search(sub_entry)
                if sub_match:
                    yield gl_id, sub_match.group("date"), sub_path, None
                else:
                    print(f"  Skipping (no pattern match): {top_entry}/{sub_entry}")
            continue

        print(f"  Skipping (unrecognised top-level entry): {top_entry}")


def build_backfill_row(tru_id, gl_id, instance, time_point, current_instance):
    """
    Minimal update row for patching an EXISTING ST instance during backfill.

    Sets time_point_spt, current_instance_spt (needed for TRU-ID calc), and
    the REDCAP_DEFAULTS fields (facility_st, sample_type_analysed_st___3) so
    that @DEFAULT values are consistently applied even when the instance
    already existed before this script version.

    Intentionally does NOT touch preservation type, dates, analyte
    selections, or any other field not listed here.

    BUG FIX (v6.1):
      Previously build_row() was used for both creates AND updates, sending
      all fields and therefore overwriting existing preservation type, dates,
      analyte selections, etc.  This minimal row avoids that.
    BUG FIX (v6.2):
      REDCAP_DEFAULTS fields were omitted from backfill rows, so existing
      instances never received facility_st / sample_type_analysed_st___3.
    """
    suffix    = get_tissue_suffix(gl_id)
    is_normal = suffix in NORMAL_SUFFIXES
    return {
        "record_id":                    tru_id,
        "redcap_repeat_instrument":     "st_spatial_transcriptomics",
        "redcap_repeat_instance":       instance,
        "time_point_spt":               time_point,
        "current_instance_spt":         current_instance,
        # @DEFAULT fields — ensure these are set on every instance
        "sample_type_analysed_st___3":  "0" if is_normal else REDCAP_DEFAULTS["sample_type_analysed_st___3"],
        "facility_st":                  REDCAP_DEFAULTS["facility_st"],
    }


def build_row(tru_id, gl_id, full_path, date_str, time_point, histology_block, instance, args):
    date_ymd   = parse_date_ymd(date_str)
    suffix     = get_tissue_suffix(gl_id)
    is_normal  = suffix in NORMAL_SUFFIXES

    # CLI dates reformatted to YYYY-MM-DD; fall back to date parsed from dir name
    submission_date    = reformat_date_to_ymd(args.submission_date)    or date_ymd
    data_received_date = reformat_date_to_ymd(args.data_received_date) or date_ymd

    # Resolve sample-type and analyte SPT codes from CLI args
    # sample_type_spt: TFF = Tumour in FFPE block, TFR = Tumour fresh/frozen
    # analyte_spt:     TIS = Tissue,               RNA = RNA
    sample_type_spt = args.sample_type_spt if not is_normal else "NOR"
    analyte_spt     = args.analyte_spt

    is_ffpe = (sample_type_spt == SAMPLE_TYPE_SPT_FFPE)
    is_tissue = (analyte_spt == ANALYTE_SPT_TISSUE)

    row = {
        # ── Repeat instrument identifiers ────────────────────────────────────
        "record_id":                           tru_id,
        "redcap_repeat_instrument":            "st_spatial_transcriptomics",
        "redcap_repeat_instance":              instance,

        # ── Time point & instance (required for TRU-ID calc) ─────────────────
        "time_point_spt":                      time_point,
        "current_instance_spt":                args.current_instance,
        "tru_id_spt":                          "",  # @CALCTEXT — REDCap computes this

        # ── Sample identification ─────────────────────────────────────────────
        "sample_id_st":                        gl_id,
        "lab_contact_st":                      args.lab_contact,

        # ── Sample type checkboxes ────────────────────────────────────────────
        "sample_type_analysed_st___1":         "0",  # Blood
        "sample_type_analysed_st___2":         "0",  # CSF
        "sample_type_analysed_st___3":         "0" if is_normal else REDCAP_DEFAULTS["sample_type_analysed_st___3"],  # T = Tumour
        "sample_type_analysed_st___4":         "1" if is_normal else "0",  # N = Normal
        "sample_type_analysed_st___5":         "0",
        "sample_type_analysed_st___6":         "0",
        "sample_type_analysed_st___7":         "0",
        "sample_type_analysed_st___8":         "0",
        # SPT dropdown: NOR for normal tissue, otherwise use CLI value (TFF or TFR)
        "sample_type_analysed_spt_2":          sample_type_spt,

        # ── Tumour preservation ───────────────────────────────────────────────
        # Fr (1) for fresh/frozen; FFPE (6) for FFPE block
        "tumour_preservation_st___1":          "0" if (is_normal or is_ffpe) else "1",  # Fr
        "tumour_preservation_st___2":          "0",
        "tumour_preservation_st___3":          "0",
        "tumour_preservation_st___4":          "0",
        "tumour_preservation_st___5":          "0",
        "tumour_preservation_st___6":          "1" if (is_ffpe and not is_normal) else "0",  # FFPE

        # ── Normal preservation (conditional on sample type N) ────────────────
        "normal_preservation_st___1":          "1" if is_normal else "0",  # Fr
        "normal_preservation_st___2":          "0",
        "normal_preservation_st___3":          "0",
        "normal_preservation_st___4":          "0",
        "normal_preservation_st___5":          "0",
        "normal_preservation_st___6":          "0",

        # ── Analyte checkboxes ────────────────────────────────────────────────
        # Tissue (FFPE runs) → check FFPE (10); RNA (fresh/frozen) → check R (2)
        "analyte_st___1":                      "0",
        "analyte_st___2":                      "0" if is_tissue else "1",   # R = RNA
        "analyte_st___3":                      "0",
        "analyte_st___4":                      "0",
        "analyte_st___5":                      "0",
        "analyte_st___6":                      "0",
        "analyte_st___7":                      "0",
        "analyte_st___8":                      "0",
        "analyte_st___9":                      "0",
        "analyte_st___10":                     "1" if is_tissue else "0",   # FFPE = Tissue
        # SPT dropdown: TIS for Tissue, RNA for RNA
        "analyte_spt_2":                       analyte_spt,

        # ── Histology / analyte prep ──────────────────────────────────────────
        "histology_block_st":                  histology_block,
        "analyte_prep_st":                     "",

        # ── Facility & protocol ───────────────────────────────────────────────
        "unique_id_st":                        "",
        "facility_st":                         REDCAP_DEFAULTS["facility_st"],
        "submission_date_st":                  submission_date,
        "facility_sample_name_st":             os.path.basename(full_path),
        "facility_contact_st":                 args.facility_contact or "",
        "facility_protocol_st":                args.facility_protocol or "",
        "panel_st":                            args.panel or "",
        "multiplexing_st":                     args.multiplexing or "",
        "raw_sample_name_st":                  gl_id,

        # ── Results ───────────────────────────────────────────────────────────
        "data_received_date_st":               data_received_date,
        "data_processed_date_st":              "",
        "additional_details_st":               full_path,
        "findings_st":                         "",
        "additional_comments_st":              suffix if suffix else "",

        # ── Form status ───────────────────────────────────────────────────────
        "st_spatial_transcriptomics_complete": "2",
    }
    return row


def main():
    parser = argparse.ArgumentParser(description="Generate REDCap upload CSV for Spatial Transcriptomics.")
    parser.add_argument("--token",               help="REDCap API Token")
    parser.add_argument("--url",                 help="REDCap API URL", default=DEFAULT_API_URL)
    parser.add_argument("--upload",              help="Upload directly to REDCap", action="store_true")
    parser.add_argument("--time-point-map",      help="Path to CSV with columns: gl_id, time_point, histology_block (optional)", default=None)
    parser.add_argument("--current-instance",    help="Current instance (default: 1)", default="1")
    parser.add_argument("--lab-contact",         help="Lab contact name", default=None)
    parser.add_argument("--data-dir",            help="Path to ST batch directory (e.g. ST_260223)", default=DEFAULT_DATA_DIR)
    parser.add_argument("--submission-date",     help="Submission date (DD-MM-YYYY or YYYY-MM-DD)", default=None)
    parser.add_argument("--data-received-date",  help="Data received date (DD-MM-YYYY or YYYY-MM-DD)", default=None)
    parser.add_argument("--facility-contact",    help="Facility contact (radio value)", default=None)
    parser.add_argument("--facility-protocol",   help="Facility protocol (radio value)", default=None)
    parser.add_argument("--panel",               help="Panel used (optional)", default=None)
    parser.add_argument("--multiplexing",        help="Multiplexing: 0=No, 1=Yes (optional)", default=None)
    parser.add_argument("--sample-type-spt",    help="Sample type SPT code: TFF=Tumour FFPE, TFR=Tumour fresh/frozen (default: TFF)", default=SAMPLE_TYPE_SPT_FFPE)
    parser.add_argument("--analyte-spt",        help="Analyte SPT code: TIS=Tissue, RNA=RNA (default: TIS)", default=ANALYTE_SPT_TISSUE)
    args = parser.parse_args()

    api_token = args.token or os.environ.get("REDCAP_API_TOKEN")
    if not api_token:
        print("Error: REDCap API Token is required.")
        return

    # Load per-sample map (time point + histology block)
    time_point_map = load_time_point_map(args.time_point_map)
    if not time_point_map:
        print("WARNING: No sample map loaded. All samples will be skipped.")

    print("Connecting to REDCap...")
    project = get_redcap_project(args.url, api_token)
    if not project:
        return

    gl_to_tru = get_redcap_mapping(project)
    if not gl_to_tru:
        print("Failed to fetch mapping. Exiting.")
        return

    # ── BUG FIX (v6.1): get_existing_st_instances now returns two structures:
    #   existing_by_sample   — list-based dict to detect duplicates safely
    #   max_inst_per_record  — per-record max, used to compute next valid integer
    existing_by_sample, max_inst_per_record = get_existing_st_instances(project)

    # Mutable counter so each new record within this run gets a unique integer
    next_inst_counter = dict(max_inst_per_record)

    data_dir = args.data_dir
    if not os.path.exists(data_dir):
        print(f"Error: Data directory '{data_dir}' not found.")
        return

    backfill_rows     = []   # Minimal updates for existing instances
    new_rows          = []   # Full rows for brand-new instances
    skipped           = []
    consumed_existing = set()  # (tru_id, gl_id) pairs whose existing instance has
                               # already been claimed this run; subsequent scans of
                               # the same gl_id (e.g. Region_1 + Region_2 dirs) will
                               # create a new instance instead of re-backfilling.

    print(f"\nScanning {data_dir}...")

    for gl_id, date_str, full_path, time_point_override in scan_data_dir(data_dir):

        sample_info = resolve_sample_info(gl_id, time_point_map)

        # When the time point is embedded in the folder name, prefer a map entry
        # keyed as <gl_id>_<T#> (e.g. JJ73_1375_T1) so that per-timepoint
        # histology blocks and redcap_id overrides are respected.
        if time_point_override:
            tp_specific_key = f"{gl_id}_{time_point_override}"
            if tp_specific_key in time_point_map:
                sample_info = time_point_map[tp_specific_key]

        # Use redcap_id alias for the REDCap lookup if one is defined in the map
        # (e.g. JJ73_1375 is registered in REDCap as GX0008)
        lookup_id = (sample_info.get("redcap_id") or gl_id) if sample_info else gl_id
        tru_id = resolve_tru_id(lookup_id, gl_to_tru)
        if not tru_id:
            alias_note = f" (looked up as {lookup_id})" if lookup_id != gl_id else ""
            print(f"  Warning: No TRU-ID found for {gl_id}{alias_note} (path: {full_path})")
            skipped.append(gl_id)
            continue

        if time_point_override:
            # Time point is embedded in the folder name (e.g. JJ73_1375_T0_Region_3).
            # The map is still consulted for histology_block but is not required.
            time_point      = time_point_override
            histology_block = sample_info["histology_block"] if sample_info else ""
        else:
            # Time point must come from the map (GL, PRISM, etc.)
            if not sample_info:
                print(f"  Skipping (not in map): {gl_id} — add to timepoint_map.csv")
                skipped.append(gl_id)
                continue
            time_point      = sample_info["time_point"]
            histology_block = sample_info["histology_block"]

        suffix = get_tissue_suffix(gl_id)

        instance_key = (str(tru_id), gl_id)
        existing     = existing_by_sample.get(instance_key, [])

        if existing and instance_key not in consumed_existing:
            # ── Ambiguous: multiple instances for the same (record, sample_id) ─
            if len(existing) > 1:
                print(f"  SKIPPING (ambiguous — instances {existing} for "
                      f"{gl_id}, TRU-ID {tru_id}): resolve manually in REDCap")
                skipped.append(gl_id)
                continue

            # ── Safe update: patch ONLY time_point_spt + current_instance_spt ──
            # Only the two fields needed to generate the TRU-ID are touched;
            # all other existing field data is left completely untouched.
            instance = str(existing[0])
            row      = build_backfill_row(tru_id, gl_id, instance, time_point, args.current_instance)
            backfill_rows.append(row)
            consumed_existing.add(instance_key)  # Claim this existing instance
            print(f"  [BACKFILL instance {instance}] {gl_id} -> TRU-ID: {tru_id}, "
                  f"Time point: {time_point}"
                  + (f", Histology: {histology_block}" if histology_block else "")
                  + (f", Tissue: {suffix}" if suffix else ""))

        else:
            # ── New instance: assign next valid integer ───────────────────────
            # Covers three cases:
            #   1. No existing instance in REDCap for this (record, gl_id)
            #   2. A second (or third…) output dir for the same gl_id, e.g.
            #      Region_1 + Region_2 both yielding the same GL ID — the first
            #      scan consumed the existing instance above; subsequent scans
            #      land here and get a fresh instance number.
            next_inst = next_inst_counter.get(str(tru_id), 0) + 1
            next_inst_counter[str(tru_id)] = next_inst
            instance  = str(next_inst)
            row       = build_row(tru_id, gl_id, full_path, date_str,
                                  time_point, histology_block, instance, args)
            new_rows.append(row)
            print(f"  [CREATE instance {instance}] {gl_id} -> TRU-ID: {tru_id}, "
                  f"Time point: {time_point}"
                  + (f", Histology: {histology_block}" if histology_block else "")
                  + (f", Tissue: {suffix}" if suffix else ""))

    if not backfill_rows and not new_rows:
        print("\nNo matching records found. Nothing to save.")
        return

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n  {len(backfill_rows)} existing instance(s) will be BACKFILLED "
          f"(time_point_spt + current_instance_spt only — no other fields changed).")
    print(f"  {len(new_rows)} new instance(s) will be CREATED with full field set.")

    # ── Column order for full rows (matches st_upload_template.xlsx) ─────────
    template_cols = [
        "record_id", "redcap_repeat_instrument", "redcap_repeat_instance",
        "time_point_spt", "current_instance_spt", "tru_id_spt",
        "sample_id_st", "lab_contact_st",
        "sample_type_analysed_st___1", "sample_type_analysed_st___2",
        "sample_type_analysed_st___3", "sample_type_analysed_st___4",
        "sample_type_analysed_st___5", "sample_type_analysed_st___6",
        "sample_type_analysed_st___7", "sample_type_analysed_st___8",
        "sample_type_analysed_spt_2",
        "tumour_preservation_st___1", "tumour_preservation_st___2",
        "tumour_preservation_st___3", "tumour_preservation_st___4",
        "tumour_preservation_st___5", "tumour_preservation_st___6",
        "normal_preservation_st___1", "normal_preservation_st___2",
        "normal_preservation_st___3", "normal_preservation_st___4",
        "normal_preservation_st___5", "normal_preservation_st___6",
        "analyte_st___1",  "analyte_st___2",  "analyte_st___3",  "analyte_st___4",
        "analyte_st___5",  "analyte_st___6",  "analyte_st___7",  "analyte_st___8",
        "analyte_st___9",  "analyte_st___10",
        "analyte_spt_2", "histology_block_st", "analyte_prep_st",
        "unique_id_st", "facility_st", "submission_date_st",
        "facility_sample_name_st", "facility_contact_st", "facility_protocol_st",
        "panel_st", "multiplexing_st", "raw_sample_name_st",
        "data_received_date_st", "data_processed_date_st",
        "additional_details_st", "findings_st", "additional_comments_st",
        "st_spatial_transcriptomics_complete",
    ]

    # Backfill columns — the 7 fields that are actually being patched
    backfill_cols = [
        "record_id", "redcap_repeat_instrument", "redcap_repeat_instance",
        "time_point_spt", "current_instance_spt",
        "sample_type_analysed_st___3", "facility_st",
    ]

    today = datetime.today().strftime("%Y-%m-%d")

    # ── Save backfill CSV ─────────────────────────────────────────────────────
    if backfill_rows:
        bf_file = OUTPUT_FILE.format(date=today, suffix="_backfill")
        df_bf   = pd.DataFrame(backfill_rows)[backfill_cols]
        df_bf.to_csv(bf_file, index=False)
        print(f"\nSaved {len(df_bf)} backfill row(s) to: {bf_file}")
        print("  (These rows update ONLY time_point_spt and current_instance_spt.)")

    # ── Save new-records CSV ──────────────────────────────────────────────────
    if new_rows:
        new_file = OUTPUT_FILE.format(date=today, suffix="_new")
        df_new   = pd.DataFrame(new_rows)
        for col in template_cols:
            if col not in df_new.columns:
                df_new[col] = ""
        df_new = df_new[template_cols]
        df_new.to_csv(new_file, index=False)
        print(f"Saved {len(df_new)} new record(s) to: {new_file}")

    if skipped:
        print(f"\nSkipped {len(skipped)} GL/GX IDs: {', '.join(skipped)}")

    # ── Upload ────────────────────────────────────────────────────────────────
    if args.upload:
        if backfill_rows:
            print(f"\nUploading {len(backfill_rows)} backfill row(s) to REDCap...")
            try:
                response = project.import_records(backfill_rows)
                print(f"  Backfill upload successful! Response: {response}")
            except Exception as e:
                print(f"  Error uploading backfill rows: {e}")

        if new_rows:
            print(f"\nUploading {len(new_rows)} new record(s) to REDCap...")
            try:
                response = project.import_records(new_rows)
                print(f"  New records upload successful! Response: {response}")
            except Exception as e:
                print(f"  Error uploading new records: {e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
