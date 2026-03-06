#!/usr/bin/env python3
"""
Export all sample IDs with their REDCap TRU-IDs (record_id).

Performs a bulk export across ALL accessible DAGs to ensure every record
is captured regardless of which data access group it belongs to.
Results are deduplicated by record_id and sorted alphabetically/numerically.

Output CSV columns:  tissue_bank_id, record_id

Usage:
    python3 export_redcap_sample_list.py
    python3 export_redcap_sample_list.py --output /path/sample_list.csv
    python3 export_redcap_sample_list.py --prefixes GL GX JJ PRISM
"""

import sys
import csv
import re
import urllib3
import argparse
from datetime import datetime
from redcap import Project, RedcapError

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

API_URL    = "https://redcap.wehi.edu.au/api/"
TOKEN_FILE = "/vast/projects/GLIMMER/scripts/redcap_api_token.txt"
OUTPUT_DIR = "/vast/projects/GLIMMER/scripts/logs"

# All DAGs the API token has read access to — mirrors report_instruments.py
ACCESSIBLE_DAGS = ["bpop_all_samples", "bcrl_only", "anheart", "glimmer"]


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


def fetch_all_samples(project, prefixes=None):
    """
    Bulk-export record_id + tissue_bank_id across ALL accessible DAGs.

    The default export context may only return records from one DAG.
    To capture every sample, we also switch into each entry of ACCESSIBLE_DAGS
    and re-export, merging results and deduplicating by record_id.

    If prefixes is given, only tissue_bank_ids starting with one of those
    prefixes are kept; otherwise every non-empty tissue_bank_id is included.
    Returns a sorted list of (tissue_bank_id, record_id) tuples.
    """
    prefix_tuple = tuple(prefixes) if prefixes else None
    seen = {}   # record_id -> tissue_bank_id  (dedup across DAGs)

    def collect(rows, label):
        new = 0
        for r in rows:
            if r.get("redcap_repeat_instrument"):
                continue
            tb_id = r.get("tissue_bank_id", "").strip()
            if not tb_id:
                continue
            if prefix_tuple and not tb_id.startswith(prefix_tuple):
                continue
            rec_id = str(r["record_id"])
            if rec_id not in seen:
                seen[rec_id] = tb_id
                new += 1
        print(f"    {label:<22}  {new:>4} new  (running total: {len(seen)})")

    print("Fetching records from REDCap across all DAGs...")

    # Default context
    try:
        collect(
            project.export_records(fields=["record_id", "tissue_bank_id"]),
            "default context",
        )
    except RedcapError as e:
        print(f"    default context:       WARNING — {e}")

    # Switch into each known DAG and re-export
    for dag in ACCESSIBLE_DAGS:
        try:
            project.switch_dag(dag)
        except Exception as e:
            print(f"    {dag:<22}  could not switch — {e}")
            continue
        try:
            collect(
                project.export_records(fields=["record_id", "tissue_bank_id"]),
                dag,
            )
        except RedcapError as e:
            print(f"    {dag:<22}  WARNING — {e}")

    # Flip to (tissue_bank_id, record_id) and sort
    samples = [(tb, rid) for rid, tb in seen.items()]

    def sort_key(item):
        tb = item[0]
        m = re.match(r"([A-Za-z_]+)(\d+)", tb)
        if m:
            return (m.group(1), int(m.group(2)))
        return (tb, 0)

    samples.sort(key=sort_key)
    sep = ', '
    scope = f"prefixes: {sep.join(prefixes)}" if prefixes else "all prefixes"
    print(f"\n  {len(samples)} total sample(s) found ({scope}).")
    return samples


def save_csv(samples, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tissue_bank_id", "record_id"])
        w.writerows(samples)
    print(f"Saved {len(samples)} row(s) to: {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Export all sample IDs and TRU-IDs from REDCap (all DAGs)."
    )
    parser.add_argument(
        "--output", "-o", metavar="FILE",
        default=None,
        help="Output CSV path (default: logs/sample_list_YYYY-MM-DD.csv)",
    )
    parser.add_argument(
        "--prefixes", nargs="+", metavar="PREFIX",
        default=None,
        help="Restrict to tissue_bank_ids starting with these prefixes "
             "(e.g. GL GX JJ). Default: include every non-empty tissue_bank_id.",
    )
    parser.add_argument("--url",        default=API_URL)
    parser.add_argument("--token-file", default=TOKEN_FILE)
    args = parser.parse_args()

    try:
        with open(args.token_file) as f:
            token = f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: Token file not found: {args.token_file}")
        sys.exit(1)

    project = connect(args.url, token)
    samples = fetch_all_samples(project, args.prefixes)

    if not samples:
        print("No matching samples found. Nothing written.")
        return

    date    = datetime.today().strftime("%Y-%m-%d")
    outpath = args.output or f"{OUTPUT_DIR}/sample_list_{date}.csv"
    save_csv(samples, outpath)
    print("Done.")


if __name__ == "__main__":
    main()
