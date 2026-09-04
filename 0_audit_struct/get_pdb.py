#!/usr/bin/env python3
"""
Download a predicted structure (PDB file) from the AlphaFold Protein
Structure Database (https://alphafold.ebi.ac.uk) given a UniProt
accession ID.
"""

import argparse
import sys
from pathlib import Path

import requests

AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
AFDB_FILE_URL = (
    "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.pdb"
)


def get_latest_pdb_url(accession: str) -> str:
    """Query the AlphaFold API to get metadata (including the current
    model version) for a given UniProt accession."""
    api_url = AFDB_API_URL.format(accession=accession)
    resp = requests.get(api_url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    if not data:
        raise ValueError(f"No AlphaFold entry found for accession '{accession}'")

    entry = data[0]
    pdb_url = entry.get("pdbUrl")
    if not pdb_url:
        raise ValueError(f"No PDB file URL returned for accession '{accession}'")

    return pdb_url


def download_pdb(accession: str, out_dir: Path, version: int | None = None) -> Path:
    """Download the AlphaFold PDB file for a given accession.

    If version is None, queries the API for the current/latest model
    version and URL. Otherwise, builds the URL directly.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    if version is None:
        pdb_url = get_latest_pdb_url(accession)
    else:
        pdb_url = AFDB_FILE_URL.format(accession=accession, version=version)

    out_path = out_dir / "alphafold_raw.pdb"

    resp = requests.get(pdb_url, timeout=60)
    if resp.status_code == 404:
        raise FileNotFoundError(
            f"PDB file not found at {pdb_url} "
            f"(check accession '{accession}' and/or model version)"
        )
    resp.raise_for_status()

    out_path.write_bytes(resp.content)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Download a PDB structure from the AlphaFold Database."
    )
    parser.add_argument(
        "accession",
        help="UniProt accession ID (e.g. P12345)",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=Path("."),
        help="Output directory for the downloaded PDB file (default: current dir)",
    )
    parser.add_argument(
        "--version",
        type=int,
        default=None,
        help="Specific AlphaFold model version to fetch (e.g. 4). "
        "If omitted, the current version is looked up via the API.",
    )
    args = parser.parse_args()

    try:
        out_path = download_pdb(args.accession, args.out_dir, args.version)
    except (requests.HTTPError, FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Downloaded: {out_path}")


if __name__ == "__main__":
    main()
