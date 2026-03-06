#!/usr/bin/env python3
"""
Quick GL number lookup — returns BCRL ID (record_id / TRU-ID) from REDCap.

Usage:
    python3 lookup_gl.py GL0229
    python3 lookup_gl.py GL0229 GL0201 GX0004
"""

import sys
import urllib3
from redcap import Project, RedcapError

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

API_URL    = "https://your-redcap-instance.org/api/"
TOKEN_FILE = "/path/to/redcap_api_token.txt"


def lookup_gl(project, gl_id):
    """Look up a GL/GX ID via filter_logic on tissue_bank_id."""
    # Try exact match first
    for query in [gl_id, gl_id.split("-")[0]]:
        try:
            records = project.export_records(
                fields=["record_id", "tissue_bank_id"],
                filter_logic=f'[tissue_bank_id] = "{query}"'
            )
            # Return first non-repeat row
            for r in records:
                if r.get("redcap_repeat_instrument") == "":
                    return r["record_id"], r["tissue_bank_id"]
        except RedcapError as e:
            print(f"  REDCap error for {query}: {e}")
    return None, None


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 lookup_gl.py <GL_NUMBER> [GL_NUMBER ...]")
        sys.exit(1)

    queries = [q.strip() for q in sys.argv[1:]]

    with open(TOKEN_FILE) as f:
        token = f.read().strip()

    print("Connecting to REDCap...")
    try:
        project = Project(API_URL, token, verify_ssl=False)
    except TypeError:
        project = Project(API_URL, token)

    print(f"\n{'GL Number':<15} {'BCRL ID (TRU-ID)':<20} {'tissue_bank_id in REDCap'}")
    print("-" * 60)

    for gl in queries:
        record_id, matched_id = lookup_gl(project, gl)
        if record_id:
            print(f"{gl:<15} {record_id:<20} {matched_id}")
        else:
            print(f"{gl:<15} {'NOT FOUND':<20} -")


if __name__ == "__main__":
    main()
