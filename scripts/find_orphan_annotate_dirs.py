#!/usr/bin/env python3
"""
Print annotate/ subdirectories that won't be produced by the current samples.csv,
i.e. leftovers from a previous version of the manifest.

Mirrors the folder-name logic from pipeline/nextflow/funannotate.nf:
    strain = strain.replaceAll(/;.*$/, '').trim()
    out    = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
"""

import csv
import os
import re
import sys

SAMPLES_CSV = os.path.join(os.path.dirname(__file__), '..', 'samples.csv')
ANNOTATE_DIR = os.path.join(os.path.dirname(__file__), '..', 'annotate')


def strip_markdown_italics(s: str) -> str:
    # Remove underscores used as markdown italics markers, e.g. "_Candida_ aaseri" → "Candida aaseri"
    return re.sub(r'(^|(?<=\s))_|_($|(?=\s))', '', s)


def compute_out(species: str, strain: str) -> str:
    species = strip_markdown_italics((species or '').strip())
    strain = (strain or '').strip()
    strain = re.sub(r';.*$', '', strain).strip()
    parts = [p for p in [species, strain] if p]
    return re.sub(r'\s+', '_', '_'.join(parts))


def main():
    samples_csv = os.path.realpath(SAMPLES_CSV)
    annotate_dir = os.path.realpath(ANNOTATE_DIR)

    # Build the set of expected output folder names from samples.csv
    expected = set()
    with open(samples_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = (row.get('ASMID') or '').strip()
            species = row.get('SPECIES') or ''
            strain = row.get('STRAIN') or ''
            out = compute_out(species, strain)
            if out and asmid:
                expected.add(out)

    # List actual subdirectories in annotate/
    try:
        actual = {
            d for d in os.listdir(annotate_dir)
            if os.path.isdir(os.path.join(annotate_dir, d))
        }
    except FileNotFoundError:
        print(f"ERROR: annotate dir not found: {annotate_dir}", file=sys.stderr)
        sys.exit(1)

    orphans = sorted(actual - expected)

    if not orphans:
        print("No orphan directories found.")
        return

    print(f"# {len(orphans)} annotate/ dirs not in current samples.csv:")
    for d in orphans:
        print(d)


if __name__ == '__main__':
    main()
