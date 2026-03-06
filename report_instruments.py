#!/usr/bin/env python3
"""
Report available REDCap instruments for a list of GL/GX/JJ sample IDs.

Summary mode (default):
  Shows completion status or instance count per instrument.

Detail mode (--detail):
  A) Prints decoded field values per instance to stdout.
  B) When --output is also given, writes alongside it:
       <name>_detail.json      — structured JSON per sample/instrument/instance
       <name>_instances.csv    — long-format CSV (one row per field value)

Usage:
    python3 report_instruments.py
    python3 report_instruments.py --detail
    python3 report_instruments.py --gl GL0283 GX0003 JJ14 --detail
    python3 report_instruments.py --gl GL0283 --detail --output /path/report.csv
"""

import sys
import re
import csv
import json
import urllib3
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from redcap import Project, RedcapError

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

API_URL    = "https://redcap.wehi.edu.au/api/"
TOKEN_FILE = "/vast/projects/GLIMMER/scripts/redcap_api_token.txt"

DEFAULT_SAMPLE_IDS = [
    "GL0405","GL0231","GL0233","GL0218-1","GL0233-1","GL0231-1", 
]

COMPLETION_LABELS = {"0": "Incomplete", "1": "Unverified", "2": "Complete", "": "No data"}
SKIP_FIELDS       = {"redcap_repeat_instrument", "redcap_repeat_instance"}

# All DAGs the API token has read access to.
# Used as fallback when a sample is not found in the default DAG context.
ACCESSIBLE_DAGS = ["bpop_all_samples", "bcrl_only", "anheart", "glimmer"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def strip_html(text):
    return re.sub(r'<[^>]+>', '', text or '').strip()


def connect(api_url, token):
    try:
        try:
            project = Project(api_url, token, verify_ssl=False)
        except TypeError:
            project = Project(api_url, token)
        print("Connected to REDCap project.")
        return project
    except Exception as e:
        print(f"ERROR connecting to REDCap: {e}")
        sys.exit(1)


def get_instruments(project):
    instruments = project.export_instruments()
    print(f"  {len(instruments)} instrument(s) found.")
    return instruments


def get_id_mapping(project):
    """
    Build {tissue_bank_id: record_id} from a bulk export.
    Uses bare parameters (matching the ST upload script) to maximise
    the number of records returned within the current DAG access.
    Returns an empty dict on failure — resolve_record_id will fall back
    to per-sample filter_logic in that case.
    """
    print("Fetching ID mapping from REDCap...")
    try:
        rows = project.export_records(fields=["record_id", "tissue_bank_id"])
        mapping = {
            r["tissue_bank_id"]: r["record_id"]
            for r in rows
            if r.get("tissue_bank_id")
        }
        print(f"  {len(mapping)} records loaded into mapping.")
        return mapping
    except RedcapError as e:
        print(f"  WARNING: Bulk mapping failed ({e}). Will use per-sample lookup.")
        return {}


def resolve_record_id(project, sample_id, mapping):
    """
    Returns (record_id, matched_id, dag_used).
    dag_used is None when found without a DAG switch, or the DAG name when
    a switch_dag() call was required. The caller must re-apply switch_dag(dag_used)
    before fetching record data so the correct DAG context is active.

    Lookup order:
    1. Pre-built bulk mapping (fast, no API call).
    2. filter_logic in current DAG context.
    3. switch_dag() through each ACCESSIBLE_DAGS entry.
    """
    for query in [sample_id, sample_id.split("-")[0]]:
        if query in mapping:
            return mapping[query], query, None

    # filter_logic fallback — searches within currently active DAG
    for query in [sample_id, sample_id.split("-")[0]]:
        try:
            records = project.export_records(
                fields=["record_id", "tissue_bank_id"],
                filter_logic=f'[tissue_bank_id] = "{query}"',
            )
            for r in records:
                if r.get("redcap_repeat_instrument") == "":
                    return r["record_id"], r["tissue_bank_id"], None
        except RedcapError:
            pass

    # DAG-switch fallback — try each accessible DAG in turn
    for dag in ACCESSIBLE_DAGS:
        try:
            project.switch_dag(dag)
        except Exception:
            continue
        for query in [sample_id, sample_id.split("-")[0]]:
            try:
                records = project.export_records(
                    fields=["record_id", "tissue_bank_id"],
                    filter_logic=f'[tissue_bank_id] = "{query}"',
                )
                for r in records:
                    if r.get("redcap_repeat_instrument") == "":
                        # Return the dag so the caller can re-apply it before fetching
                        return r["record_id"], r["tissue_bank_id"], dag
            except RedcapError:
                pass

    return None, None, None


def build_field_info(project):
    """
    Returns {field_name: {label, form, type, choices}} from project metadata.
    Falls back to None if metadata export fails (detail will use raw values).
    """
    try:
        meta = project.export_metadata()
        field_info = {}
        for f in meta:
            choices = {}
            raw = f.get("select_choices_or_calculations", "")
            if raw and f["field_type"] in ("radio", "dropdown", "checkbox"):
                for part in raw.split("|"):
                    kv = part.strip().split(",", 1)
                    if len(kv) == 2:
                        choices[kv[0].strip()] = kv[1].strip()
            field_info[f["field_name"]] = {
                "label":   strip_html(f["field_label"]),
                "form":    f["form_name"],
                "type":    f["field_type"],
                "choices": choices,
            }
        print(f"  {len(field_info)} fields loaded from metadata.")
        return field_info
    except Exception as e:
        print(f"  WARNING: Could not load metadata ({e}). Detail view will show raw field names.")
        return None


def fetch_record_data(project, record_id):
    try:
        return project.export_records(
            records=[str(record_id)],
            export_survey_fields=False,
            export_data_access_groups=False,
        )
    except RedcapError as e:
        print(f"  WARNING: Could not fetch data for record {record_id}: {e}")
        return []


# ── Summarise ─────────────────────────────────────────────────────────────────

def summarise_instruments(rows, instruments):
    """
    Returns (summary_dict, repeat_rows_by_instrument, base_row).
    Repeating detection is derived from the row data itself.
    """
    repeat_counts = defaultdict(int)
    repeat_rows   = defaultdict(list)
    base_complete = {}
    base_row      = {}

    for row in rows:
        ri = row.get("redcap_repeat_instrument", "")
        if ri:
            repeat_counts[ri] += 1
            repeat_rows[ri].append(row)
        else:
            base_row = row
            for k, v in row.items():
                if k.endswith("_complete"):
                    base_complete[k[:-len("_complete")]] = v

    summary = {}
    for instr in instruments:
        name = instr["instrument_name"]
        if name in repeat_counts:
            summary[name] = f"{repeat_counts[name]} instance(s)"
        elif name in base_complete:
            summary[name] = COMPLETION_LABELS.get(base_complete[name], f"code={base_complete[name]}")
        else:
            summary[name] = "No data"

    return summary, dict(repeat_rows), base_row


# ── Field decoding ────────────────────────────────────────────────────────────

def decode_row(row, form_name, field_info):
    """
    Return list of (label, value) for all non-empty fields in the row
    that belong to form_name.  Falls back to (field_name, raw_value) when
    field_info is None.  Checkbox sub-fields are grouped into one entry.
    """
    output        = []
    checkbox_grps = defaultdict(list)

    for field_name, raw_value in row.items():
        if field_name in SKIP_FIELDS or field_name.endswith("_complete"):
            continue
        if raw_value == "" or raw_value is None:
            continue

        # ── No metadata: raw fallback ─────────────────────────────────────
        if field_info is None:
            output.append((field_name, raw_value))
            continue

        # ── Checkbox sub-field ────────────────────────────────────────────
        if "___" in field_name:
            base, code = field_name.rsplit("___", 1)
            fi = field_info.get(base)
            if fi and fi["type"] == "checkbox" and fi["form"] == form_name and raw_value == "1":
                checkbox_grps[base].append(fi["choices"].get(code, code))
            continue

        fi = field_info.get(field_name)
        if not fi or fi["form"] != form_name:
            continue
        if fi["type"] in ("descriptive", "file"):
            continue

        decoded = fi["choices"].get(raw_value, raw_value) if fi["type"] in ("radio", "dropdown") else raw_value
        output.append((fi["label"] or field_name, decoded))

    for base, selected in checkbox_grps.items():
        fi = field_info.get(base) if field_info else None
        output.append((fi["label"] if fi else base, ", ".join(selected)))

    return output


# ── Output A: stdout ──────────────────────────────────────────────────────────

def print_summary_table(results, instruments):
    instr_w = max(len(i["instrument_label"]) for i in instruments)
    sep = "-" * (instr_w + 24)
    for r in results:
        print(f"\n{sep}")
        print(f"  Sample ID : {r['sample_id']}")
        print(f"  TRU-ID    : {r['record_id'] or 'NOT FOUND'}")
        if not r["record_id"]:
            continue
        print(f"  {'Instrument':<{instr_w}}  Status / Instances")
        print(f"  {'-'*instr_w}  {'-'*14}")
        for instr in instruments:
            print(f"  {instr['instrument_label']:<{instr_w}}  {r['instrument_summary'].get(instr['instrument_name'], '-')}")
    print(sep)


def print_detail_table(results, instruments, field_info):
    sep = "=" * 70
    for r in results:
        print(f"\n{sep}")
        print(f"  Sample ID : {r['sample_id']}")
        print(f"  TRU-ID    : {r['record_id'] or 'NOT FOUND'}")
        if not r["record_id"]:
            continue

        for instr in instruments:
            name  = instr["instrument_name"]
            label = instr["instrument_label"]
            val   = r["instrument_summary"].get(name, "No data")
            if val == "No data":
                continue

            print(f"\n  ── {label} ──")
            if name in r["repeat_rows"]:
                for row in r["repeat_rows"][name]:
                    inst_num = row.get("redcap_repeat_instance", "?")
                    print(f"    [Instance {inst_num}]")
                    fields = decode_row(row, name, field_info)
                    if fields:
                        lw = max(len(l) for l, _ in fields)
                        for lbl, v in fields:
                            print(f"      {lbl:<{lw}}  :  {v}")
                    else:
                        print("      (no decoded fields)")
            else:
                fields = decode_row(r["base_row"], name, field_info)
                if fields:
                    lw = max(len(l) for l, _ in fields)
                    for lbl, v in fields:
                        print(f"    {lbl:<{lw}}  :  {v}")
                else:
                    print(f"    Status: {val}")
        print(f"\n{sep}")


# ── Output B: files ───────────────────────────────────────────────────────────

def build_detail_structure(results, instruments, field_info):
    """Build a structured list suitable for JSON export."""
    instr_label = {i["instrument_name"]: i["instrument_label"] for i in instruments}
    output = []
    for r in results:
        entry = {
            "sample_id": r["sample_id"],
            "record_id": r["record_id"],
            "instruments": {},
        }
        if not r["record_id"]:
            output.append(entry)
            continue

        for instr in instruments:
            name  = instr["instrument_name"]
            label = instr["instrument_label"]
            val   = r["instrument_summary"].get(name, "No data")
            if val == "No data":
                continue

            instr_entry = {"status": val, "instances": []}
            if name in r["repeat_rows"]:
                for row in r["repeat_rows"][name]:
                    inst_num = row.get("redcap_repeat_instance", "?")
                    fields   = decode_row(row, name, field_info)
                    instr_entry["instances"].append({
                        "instance_number": inst_num,
                        "fields": {lbl: v for lbl, v in fields},
                    })
            else:
                fields = decode_row(r["base_row"], name, field_info)
                instr_entry["instances"].append({
                    "instance_number": "base",
                    "fields": {lbl: v for lbl, v in fields},
                })
            entry["instruments"][label] = instr_entry

        output.append(entry)
    return output


def save_summary_csv(results, instruments, path):
    fieldnames = ["sample_id", "record_id"] + [i["instrument_name"] for i in instruments]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            row = {"sample_id": r["sample_id"], "record_id": r["record_id"] or "NOT FOUND"}
            if r["record_id"]:
                row.update(r["instrument_summary"])
            w.writerow(row)
    print(f"  Summary CSV  : {path}")


def save_detail_json(detail_data, path):
    with open(path, "w") as f:
        json.dump(detail_data, f, indent=2)
    print(f"  Detail JSON  : {path}")


def save_instances_csv(detail_data, path):
    """Long-format CSV: one row per field value across all instances."""
    rows = []
    for sample in detail_data:
        sid = sample["sample_id"]
        rid = sample["record_id"] or "NOT FOUND"
        for instr_label, instr_data in sample.get("instruments", {}).items():
            for inst in instr_data.get("instances", []):
                inst_num = inst["instance_number"]
                for field_label, value in inst.get("fields", {}).items():
                    rows.append({
                        "sample_id":       sid,
                        "record_id":       rid,
                        "instrument":      instr_label,
                        "instance_number": inst_num,
                        "field":           field_label,
                        "value":           value,
                    })

    fieldnames = ["sample_id", "record_id", "instrument", "instance_number", "field", "value"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"  Instances CSV: {path}  ({len(rows)} rows)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Report REDCap instruments for GL/GX/JJ sample IDs."
    )
    parser.add_argument("--gl", nargs="+", metavar="SAMPLE_ID",
                        help="GL/GX/JJ IDs to query (default: built-in list)")
    parser.add_argument("--detail", "-d", action="store_true",
                        help="Show decoded field values per instance (stdout + files if --output given)")
    parser.add_argument("--output", "-o", metavar="FILE",
                        help="Save summary CSV; with --detail also writes _detail.json and _instances.csv")
    parser.add_argument("--url",        default=API_URL)
    parser.add_argument("--token-file", default=TOKEN_FILE)
    args = parser.parse_args()

    sample_ids = args.gl or DEFAULT_SAMPLE_IDS

    try:
        with open(args.token_file) as f:
            token = f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: Token file not found: {args.token_file}")
        sys.exit(1)

    project    = connect(args.url, token)
    instruments = get_instruments(project)
    id_mapping  = get_id_mapping(project)

    field_info = None
    if args.detail:
        print("Fetching field metadata for detail view...")
        field_info = build_field_info(project)

    results = []
    print(f"\nLooking up {len(sample_ids)} sample ID(s)...\n")

    for sample_id in sample_ids:
        record_id, matched_id, dag_used = resolve_record_id(project, sample_id, id_mapping)
        if not record_id:
            print(f"  {sample_id:<15} NOT FOUND in REDCap")
            results.append({"sample_id": sample_id, "record_id": None,
                             "instrument_summary": {}, "repeat_rows": {}, "base_row": {}})
            continue

        # If the record was found via switch_dag, re-apply that DAG context
        # so fetch_record_data can actually read the record's data.
        if dag_used:
            try:
                project.switch_dag(dag_used)
            except Exception:
                pass

        rows = fetch_record_data(project, record_id)
        summary, repeat_rows, base_row = summarise_instruments(rows, instruments)
        with_data = sum(1 for v in summary.values() if v != "No data")
        print(f"  {sample_id:<15} -> record_id={record_id:<6} | {with_data} instrument(s) with data")

        results.append({
            "sample_id":          sample_id,
            "record_id":          record_id,
            "instrument_summary": summary,
            "repeat_rows":        repeat_rows,
            "base_row":           base_row,
        })

    # ── Output ────────────────────────────────────────────────────────────────
    if args.detail:
        print_detail_table(results, instruments, field_info)
    else:
        print_summary_table(results, instruments)

    if args.output:
        base = Path(args.output)
        stem = base.stem
        parent = base.parent
        print("\nSaving files...")
        save_summary_csv(results, instruments, base)

        if args.detail:
            detail_data = build_detail_structure(results, instruments, field_info)
            save_detail_json(detail_data,   parent / f"{stem}_detail.json")
            save_instances_csv(detail_data, parent / f"{stem}_instances.csv")

    print("\nDone.")


if __name__ == "__main__":
    main()
