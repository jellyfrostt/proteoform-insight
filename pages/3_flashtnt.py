"""Page 3: FLASHTnT Results -- Proteoform/Tag Tables + Sequence Coverage.

FLASHTnT integration and proteoform-level results.
"""

import shutil

import polars as pl
import streamlit as st
from openms_insight import LinePlot, SequenceView, Table

from flashtnt_parser import (
    compute_tag_coverage,
    has_flashtnt_pyopenms,
    parse_flashtnt_output,
    run_flashtnt_cli,
    run_flashtnt_pyopenms,
)
from init import init_app
init_app()

sm = st.session_state["state_manager"]

_src = st.session_state.get("_data_source", "demo")
_prot = st.session_state.get("_selected_protein", "default")
_cache_suffix = f"{_src}_{_prot}"

st.header("FLASHTnT Results")

# -- Sidebar: data source ---------------------------------------------
with st.sidebar:
    st.subheader("Data Source")
    source = st.radio(
        "Source",
        ["Demo data", "Upload TSV files", "Run FLASHTnT"],
        index=0,
    )

if source == "Upload TSV files":
    col1, col2 = st.columns(2)
    protein_file = col1.file_uploader("protein.tsv", type=["tsv"])
    tag_file = col2.file_uploader("tags.tsv", type=["tsv"])

    if protein_file and tag_file:
        try:
            protein_df, tag_df = parse_flashtnt_output(protein_file, tag_file)
            st.session_state["protein_df"] = protein_df
            st.session_state["tag_df"] = tag_df
            st.success("Parsed FLASHTnT output successfully.")
        except Exception as e:
            st.error(f"Error parsing files: {e}")
            st.stop()
    else:
        st.info("Upload both protein.tsv and tags.tsv from FLASHTnT output.")
        st.stop()

elif source == "Run FLASHTnT":
    # -- On-demand FLASHTnT targeted search -----------------------------
    st.subheader("On-Demand Targeted Search")

    has_pyoms_tnt = has_flashtnt_pyopenms()
    has_cli_tnt = shutil.which("FLASHTnT") is not None
    has_any_tnt = has_pyoms_tnt or has_cli_tnt

    # Show backend status in sidebar
    with st.sidebar:
        if has_pyoms_tnt:
            st.caption("FLASHTnT: Python bindings (pyOpenMS)")
        elif has_cli_tnt:
            st.caption("FLASHTnT: CLI backend")
        else:
            st.caption("FLASHTnT: not available")

    if not has_any_tnt:
        st.warning(
            "FLASHTnT not found (neither Python bindings nor CLI). "
            "Install OpenMS with FLASHTnT or use Demo/Upload mode."
        )

    col1, col2 = st.columns(2)
    mzml_file = col1.file_uploader("Deconvolved mzML", type=["mzML"])
    fasta_file = col2.file_uploader("FASTA database", type=["fasta", "fa"])

    selected_spectrum = sm.get_selection("spectrum")
    if selected_spectrum is not None:
        st.info(f"Selected spectrum: Scan {selected_spectrum}")
    else:
        st.caption("Select a spectrum from the LC-MS Overview page first.")

    if st.button("Run FLASHTnT", disabled=not has_any_tnt):
        if mzml_file and fasta_file:
            import tempfile
            with tempfile.TemporaryDirectory() as tmpdir:
                mzml_path = f"{tmpdir}/input.mzML"
                fasta_path = f"{tmpdir}/input.fasta"
                with open(mzml_path, "wb") as f:
                    f.write(mzml_file.read())
                with open(fasta_path, "wb") as f:
                    f.write(fasta_file.read())

                result = None
                # Try Python bindings first, fall back to CLI
                if has_pyoms_tnt:
                    with st.spinner("Running FLASHTnT via Python bindings..."):
                        result = run_flashtnt_pyopenms(mzml_path, fasta_path, tmpdir)

                if result is None and has_cli_tnt:
                    with st.spinner("Running FLASHTnT via CLI..."):
                        result = run_flashtnt_cli(mzml_path, fasta_path, tmpdir)

                if result:
                    protein_df, tag_df = parse_flashtnt_output(
                        result["protein"], result["tag"]
                    )
                    st.session_state["protein_df"] = protein_df
                    st.session_state["tag_df"] = tag_df
                    st.success("FLASHTnT search completed.")
                    st.rerun()
                else:
                    st.error("FLASHTnT search failed. Check input files.")
        else:
            st.error("Upload both mzML and FASTA files.")

    if not has_any_tnt:
        st.divider()
        st.caption("Falling back to demo data for preview.")

protein_df = st.session_state["protein_df"]
tag_df = st.session_state["tag_df"]

# -- Proteoform Table -------------------------------------------------
st.subheader("Proteoforms")

proteoform_table = Table(
    cache_id=f"flashtnt_proteoforms_{_cache_suffix}",
    data=protein_df,
    interactivity={"proteoform": "ProteoformIndex"},
    column_definitions=[
        {"field": "ProteoformIndex", "title": "Index", "sorter": "number"},
        {"field": "ProteinAccession", "title": "Accession", "sorter": "string"},
        {"field": "Score", "title": "Score", "sorter": "number"},
        {"field": "Coverage(%)", "title": "Coverage (%)", "sorter": "number"},
        {"field": "ModCount", "title": "Mods", "sorter": "number"},
        {"field": "ProteoformMass", "title": "Mass (Da)", "sorter": "number"},
        {"field": "ProteoformLevelQvalue", "title": "Q-value", "sorter": "number"},
    ],
    title="FLASHTnT Proteoforms",
    index_field="ProteoformIndex",
    default_row=0,
    initial_sort=[{"column": "Score", "dir": "desc"}],
    cache_path=".cache",
)

proteoform_table(key="pf_table", state_manager=sm)

# -- Tag Table (filtered by selected proteoform) ---------------------
st.subheader("Sequence Tags")

tag_table = Table(
    cache_id=f"flashtnt_tags_{_cache_suffix}",
    data=tag_df,
    filters={"proteoform": "ProteoformIndex"},
    column_definitions=[
        {"field": "TagIndex", "title": "Tag ID", "sorter": "number"},
        {"field": "TagSequence", "title": "Sequence", "sorter": "string"},
        {"field": "StartPosition", "title": "Start", "sorter": "number"},
        {"field": "Length", "title": "Length", "sorter": "number"},
        {"field": "DeNovoScore", "title": "Score", "sorter": "number"},
    ],
    title="Sequence Tags",
    index_field="TagIndex",
    initial_sort=[{"column": "DeNovoScore", "dir": "desc"}],
    cache_path=".cache",
)

tag_table(key="tag_table", state_manager=sm)

# -- Proteoform detail (conditional on selection) ---------------------
selected_pf = sm.get_selection("proteoform")

if selected_pf is not None:
    pf_data = protein_df.filter(pl.col("ProteoformIndex") == selected_pf).collect()

    if len(pf_data) > 0:
        row = pf_data.row(0, named=True)
        pf_seq = row.get("ProteinSequence", st.session_state["cytc_sequence"])

        st.divider()
        st.subheader(f"Proteoform {selected_pf} Detail")

        # Modification info
        mod_count = row.get("ModCount", 0)
        if mod_count and mod_count > 0:
            mod_id = row.get("ModID", "")
            mod_mass = row.get("ModMass", "")
            mod_start = row.get("ModStart", "")
            st.markdown(
                f"**Modifications:** {mod_count} -- {mod_id} "
                f"(+{mod_mass} Da at position {mod_start})"
            )

        col1, col2, col3 = st.columns(3)
        col1.metric("Score", f"{row.get('Score', 0):.1f}")
        col2.metric("Coverage", f"{row.get('Coverage(%)', 0):.1f}%")
        col3.metric("Mass", f"{row.get('ProteoformMass', 0):.2f} Da")

        # -- SequenceView for selected proteoform ---------------------
        # Use best MS2 scan for fragment annotation
        ms2_ids = st.session_state["ms2_scan_ids"]
        best_scan = ms2_ids[len(ms2_ids) // 2] if ms2_ids else None

        if best_scan:
            sm.set_selection("spectrum", best_scan)

            sv = SequenceView(
                cache_id=f"flashtnt_sv_{selected_pf}",
                sequence_data=(pf_seq, 14),
                peaks_data=st.session_state["peaks_table"],
                filters={"spectrum": "scan_id"},
                deconvolved=False,
                annotation_config={
                    "ion_types": ["c", "z"],
                    "tolerance": 20.0,
                    "tolerance_ppm": True,
                    "colors": {"c": "#E67E22", "z": "#2ECC71"},
                },
                title=f"Proteoform {selected_pf} -- Scan {best_scan}",
                cache_path=".cache",
            )
            sv(key="flashtnt_sv", state_manager=sm)

            # Annotated spectrum linked to FLASHTnT SequenceView
            st.subheader("Annotated Spectrum")
            flashtnt_annotated = LinePlot.from_sequence_view(
                sequence_view=sv,
                cache_id=f"flashtnt_annotated_{selected_pf}_{best_scan}",
                cache_path=".cache",
                title=f"Annotated Spectrum -- Proteoform {selected_pf}",
                x_label="m/z",
                y_label="Intensity",
                styling={
                    "highlightColor": "#E4572E",
                    "unhighlightedColor": "#94A3B8",
                    "selectedColor": "#F59E0B",
                },
            )
            flashtnt_annotated(
                key="flashtnt_annotated_spec",
                state_manager=sm,
                sequence_view_key="flashtnt_sv",
            )

        # -- Tag Coverage Bar Plot ------------------------------------
        st.subheader("Per-Residue Tag Coverage")

        coverage = compute_tag_coverage(pf_seq, tag_df, selected_pf)
        coverage_rows = [
            {
                "peak_id": i,
                "scan_id": 0,
                "mass": float(i),
                "intensity": cov,
                "is_covered": 1 if cov > 0 else 0,
                "label": pf_seq[i] if cov > 0.5 else "",
            }
            for i, cov in enumerate(coverage)
        ]

        coverage_df = pl.LazyFrame(coverage_rows)

        coverage_plot = LinePlot(
            cache_id=f"tag_coverage_{selected_pf}",
            data=coverage_df,
            x_column="mass",
            y_column="intensity",
            highlight_column="is_covered",
            annotation_column="label",
            title=f"Tag Coverage -- Proteoform {selected_pf}",
            x_label="Residue Position",
            y_label="Normalized Coverage",
            styling={
                "highlightColor": "#3B82F6",
                "unhighlightedColor": "#E2E8F0",
            },
            cache_path=".cache",
        )

        coverage_plot(key="coverage_plot", state_manager=sm)
