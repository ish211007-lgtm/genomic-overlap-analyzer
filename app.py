"""Streamlit web app: Genomic Feature Overlap Analyzer (BT 2113)."""
import streamlit as st

import overlap_core as oc

st.set_page_config(page_title="Genomic Feature Overlap Analyzer", page_icon="🧬", layout="wide")
st.title("🧬 Genomic Feature Overlap Analyzer")
st.caption("Find which features in file A (e.g. genes) overlap features in file B "
           "(e.g. SNPs). Coordinates are BED style: 0-based, half-open.")


def load_text(uploaded, sample_path):
    if uploaded is not None:
        return uploaded.getvalue().decode("utf-8")
    with open(sample_path) as f:
        return f.read()


with st.sidebar:
    st.header("Input")
    up_a = st.file_uploader("Feature file A (BED)", type=["bed", "txt", "tsv"])
    up_b = st.file_uploader("Feature file B (BED)", type=["bed", "txt", "tsv"])
    st.caption("No upload? The sample genes/SNPs files are used.")
    min_ov = st.number_input("Minimum overlap (bp)", min_value=1, value=1, step=1)
    up_ref = st.file_uploader("Optional reference FASTA", type=["fa", "fasta", "fna"])
    st.caption("Add a FASTA to get the overlap sequence and GC content (Biopython).")

text_a = load_text(up_a, "data/genes.bed")
text_b = load_text(up_b, "data/snps.bed")
features_a, errors_a = oc.parse_bed_text(text_a)
features_b, errors_b = oc.parse_bed_text(text_b)

for label, errs in (("File A", errors_a), ("File B", errors_b)):
    if errs:
        with st.expander(f"⚠️ {label}: {len(errs)} line(s) skipped"):
            for e in errs:
                st.write(e)

if not features_a or not features_b:
    st.error("Both files need at least one valid feature.")
    st.stop()

overlaps = oc.find_all_overlaps(features_a, features_b, min_overlap=int(min_ov))
df = oc.overlaps_to_dataframe(overlaps)

if up_ref is not None:
    import io
    reference = oc.load_reference_fasta(io.StringIO(up_ref.getvalue().decode("utf-8")))
    df = oc.add_sequence_info(df, reference)

stats = oc.summarise(df, features_a, features_b)
cols = st.columns(len(stats))
for col, (k, v) in zip(cols, stats.items()):
    col.metric(k, v)

st.subheader("Overlapping pairs")
if df.empty:
    st.info("No overlaps found.")
else:
    st.dataframe(df, width="stretch", hide_index=True)
    st.download_button("Download results (CSV)", df.to_csv(index=False),
                       "overlaps.csv", "text/csv")

missing = oc.features_without_overlap(features_a, overlaps)
st.subheader("Features in A with no overlap")
st.write(", ".join(missing) if missing else "None — every feature in A overlaps something in B.")

with st.expander("Show parsed input"):
    c1, c2 = st.columns(2)
    c1.write("File A"); c1.dataframe(features_a, hide_index=True)
    c2.write("File B"); c2.dataframe(features_b, hide_index=True)
