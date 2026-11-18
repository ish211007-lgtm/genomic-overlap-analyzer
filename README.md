# 🧬 Genomic Feature Overlap Analyzer

BT 2113 – Python & Biopython (Fall 2026) class project.

Upload two BED files (e.g. genes and SNPs) and the app reports every pair of
features that overlap, how many bases they share, what fraction of each feature
that is, and how they relate (identical / inside / partial). Optionally add a
reference FASTA to see the overlap sequence and its GC content.

**Live app:** _paste your Streamlit link here_

## Run locally
```
pip install -r requirements.txt
streamlit run app.py
pytest tests
```

## Files
| File | Purpose |
|---|---|
| `overlap_core.py` | All logic: parsing, overlap detection, summary, Biopython step |
| `app.py` | Streamlit web interface |
| `data/` | Sample genes.bed, snps.bed, reference.fa |
| `tests/` | Unit tests |

## Course topics used
| Syllabus topic | Where |
|---|---|
| Variables, loops, conditionals, functions | `overlap_core.py` throughout |
| Dictionaries | features stored as dicts; `group_by_chrom` |
| File I/O, parsing delimited files | `read_feature_file`, `parse_bed_text` |
| Regular expressions (`re`) | chromosome-name check, whitespace splitting |
| Biopython `SeqIO`, `Seq` | `load_reference_fasta`, `add_sequence_info` |
| NumPy | `summarise` (sum / mean / max) |
| Pandas | results table, CSV export |

## Method
Features use BED coordinates (0-based, half-open). Features A and B overlap if
they are on the same chromosome and `A.start < B.end` and `B.start < A.end`.
Overlap length = `min(A.end, B.end) - max(A.start, B.start)`. Features that only
touch (A.end == B.start) do not overlap.

## Team
- Member 1: _name – role_
- Member 2: _name – role_
