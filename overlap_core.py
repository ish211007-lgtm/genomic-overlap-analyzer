"""Core logic for the Genomic Feature Overlap Analyzer (BT 2113).

Only topics from the course are used: functions, dictionaries, loops,
conditionals, file I/O, the re module, Biopython (SeqIO / Seq), NumPy, Pandas.

Coordinates follow the BED convention: 0-based, half-open [start, end).
Two features overlap only if they share at least one base.
"""
import re

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# chromosome names like chr1, chrX, chrM, 1, X, MT
CHROM_PATTERN = re.compile(r"^(chr)?([0-9]{1,2}|[XYM]|MT)$", re.IGNORECASE)


# ---------- 1. Representing features (dictionaries) ----------
def make_feature(chrom, start, end, name):
    """A genomic feature is stored as a dictionary."""
    return {"chrom": chrom, "start": start, "end": end, "name": name}


# ---------- 2. Reading / parsing files (file I/O + re) ----------
def parse_bed_text(text):
    """Parse BED-style text (tab or space separated: chrom start end name).

    Returns (features, errors). Bad lines are reported, not fatal.
    """
    features = []
    errors = []
    for line_no, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if line == "" or line.startswith(("#", "track", "browser")):
            continue
        parts = re.split(r"\s+", line)
        if len(parts) < 3:
            errors.append(f"Line {line_no}: needs at least chrom, start, end")
            continue
        chrom = parts[0]
        if not CHROM_PATTERN.match(chrom):
            errors.append(f"Line {line_no}: unrecognised chromosome '{chrom}'")
            continue
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            errors.append(f"Line {line_no}: start/end must be integers")
            continue
        if start < 0 or end <= start:
            errors.append(f"Line {line_no}: need 0 <= start < end")
            continue
        name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
        features.append(make_feature(chrom, start, end, name))
    return features, errors


def read_feature_file(filepath):
    """Read a BED file from disk and parse it."""
    with open(filepath, "r") as f:
        return parse_bed_text(f.read())


# ---------- 3. Overlap logic (functions, conditionals) ----------
def features_overlap(a, b):
    if a["chrom"] != b["chrom"]:
        return False
    return a["start"] < b["end"] and b["start"] < a["end"]


def overlap_length(a, b):
    return min(a["end"], b["end"]) - max(a["start"], b["start"])


def group_by_chrom(features):
    """Dictionary: chromosome -> list of features (avoids useless comparisons)."""
    groups = {}
    for feat in features:
        groups.setdefault(feat["chrom"], []).append(feat)
    return groups


def find_all_overlaps(features_a, features_b, min_overlap=1):
    """Return a list of dictionaries, one per overlapping pair."""
    results = []
    groups_b = group_by_chrom(features_b)
    for a in features_a:
        for b in groups_b.get(a["chrom"], []):
            if features_overlap(a, b):
                length = overlap_length(a, b)
                if length >= min_overlap:
                    results.append({
                        "chrom": a["chrom"],
                        "feature_A": a["name"],
                        "feature_B": b["name"],
                        "overlap_start": max(a["start"], b["start"]),
                        "overlap_end": min(a["end"], b["end"]),
                        "overlap_bp": length,
                        "pct_of_A": round(100 * length / (a["end"] - a["start"]), 2),
                        "pct_of_B": round(100 * length / (b["end"] - b["start"]), 2),
                        "relationship": relationship(a, b),
                    })
    return results


def relationship(a, b):
    """Describe how A and B overlap."""
    if a["start"] == b["start"] and a["end"] == b["end"]:
        return "identical"
    if a["start"] <= b["start"] and b["end"] <= a["end"]:
        return "B inside A"
    if b["start"] <= a["start"] and a["end"] <= b["end"]:
        return "A inside B"
    return "partial"


def features_without_overlap(features_a, overlaps):
    """Names of features in A that have no partner in B."""
    hit = {row["feature_A"] for row in overlaps}
    return [f["name"] for f in features_a if f["name"] not in hit]


# ---------- 4. Summary with NumPy / Pandas ----------
def overlaps_to_dataframe(overlaps):
    columns = ["chrom", "feature_A", "feature_B", "overlap_start", "overlap_end",
               "overlap_bp", "pct_of_A", "pct_of_B", "relationship"]
    return pd.DataFrame(overlaps, columns=columns)


def summarise(overlaps_df, features_a, features_b):
    lengths = overlaps_df["overlap_bp"].to_numpy()
    return {
        "Features in A": len(features_a),
        "Features in B": len(features_b),
        "Overlapping pairs": len(overlaps_df),
        "Total overlap (bp)": int(np.sum(lengths)) if len(lengths) else 0,
        "Mean overlap (bp)": round(float(np.mean(lengths)), 1) if len(lengths) else 0,
        "Max overlap (bp)": int(np.max(lengths)) if len(lengths) else 0,
    }


# ---------- 5. Biopython: sequence of the overlapping region ----------
def load_reference_fasta(handle):
    """Read a (multi-)FASTA into a dictionary {record id: Seq} using SeqIO."""
    return {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}


def add_sequence_info(overlaps_df, reference):
    """Add the overlap sequence and its GC content from a reference genome."""
    seqs, gcs = [], []
    for _, row in overlaps_df.iterrows():
        seq = reference.get(row["chrom"])
        if seq is None or row["overlap_end"] > len(seq):
            seqs.append("n/a")
            gcs.append(np.nan)
        else:
            piece = seq[row["overlap_start"]:row["overlap_end"]]
            seqs.append(str(piece) if len(piece) <= 50 else str(piece[:50]) + "...")
            gcs.append(round(100 * gc_fraction(piece), 1))
    out = overlaps_df.copy()
    out["overlap_sequence"] = seqs
    out["GC_percent"] = gcs
    return out
