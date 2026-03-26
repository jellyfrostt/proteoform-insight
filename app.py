"""Proteoform Insight -- Interactive top-down proteomics visualization.

Built on OpenMS-Insight + Streamlit for top-down proteomics visualization.
Demonstrates proteoform-centric visualization and targeted spectrum analysis.
"""

import shutil
from pathlib import Path

import streamlit as st

from demo_data import PROTEINS
from pyopenms_integration import has_pyopenms, get_pyopenms_version


def _clear_cache():
    """Remove stale component cache on data source switch."""
    cache_dir = Path(__file__).parent / ".cache"
    if cache_dir.exists():
        shutil.rmtree(cache_dir, ignore_errors=True)


def _load_sample_data():
    """Load bundled Cytochrome C sample data.

    Uses pyOpenMS to parse mzML if available, otherwise loads
    pre-processed parquet files.
    """
    import json
    import polars as pl
    from demo_data import CYTC_SEQUENCE

    data_dir = Path(__file__).parent / "data"

    if has_pyopenms():
        mzml_path = data_dir / "FLASHDeconv_sample_input1.mzML"
        if mzml_path.exists():
            from pyopenms_integration import load_mzml
            with st.spinner("Loading sample mzML (Cytochrome C, ETD)..."):
                result = load_mzml(mzml_path.read_bytes())
            for key, value in result.items():
                st.session_state[key] = value
    else:
        # Load pre-processed parquet files
        st.session_state["spectra_table"] = pl.scan_parquet(str(data_dir / "sample_spectra.parquet"))
        st.session_state["peaks_table"] = pl.scan_parquet(str(data_dir / "sample_peaks.parquet"))
        st.session_state["heatmap_data"] = pl.scan_parquet(str(data_dir / "sample_heatmap.parquet"))
        meta = json.loads((data_dir / "sample_meta.json").read_text())
        st.session_state["ms2_scan_ids"] = meta["ms2_scan_ids"]

    st.session_state["cytc_sequence"] = CYTC_SEQUENCE
    st.session_state["_selected_protein"] = "Cytochrome C"
    st.session_state["_data_source"] = "sample_mzml"
    st.session_state["_demo_data_ready"] = True
    st.rerun()

st.set_page_config(
    page_title="Proteoform Insight",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

from init import init_app

init_app()

# Navigation
pages = st.navigation([
    st.Page("pages/1_overview.py", title="LC-MS Overview"),
    st.Page("pages/2_proteoform.py", title="Proteoform Analysis"),
    st.Page("pages/3_flashtnt.py", title="FLASHTnT Results"),
    st.Page("pages/4_feature_trace.py", title="Feature Traces"),
])

# Sidebar branding
with st.sidebar:
    st.markdown("### Proteoform Insight")
    st.caption("Top-Down Proteomics Visualization")
    st.divider()

    # pyOpenMS status
    pyoms_version = get_pyopenms_version()
    if pyoms_version:
        st.success(f"pyOpenMS {pyoms_version}")
    else:
        st.warning(
            "pyOpenMS not found — running in demo mode.  \n"
            "`pip install pyopenms>=3.0`"
        )

    # Data source selection
    if has_pyopenms():
        source_options = ["Sample mzML (Cytochrome C)", "Demo data (synthetic)", "Upload mzML"]
    else:
        source_options = ["Sample mzML (Cytochrome C)", "Demo data (synthetic)"]
    data_source = st.radio("Data source", source_options, horizontal=True)

    if data_source == "Demo data (synthetic)":
        if st.session_state.get("_data_source") == "sample_mzml":
            st.session_state["_data_source"] = None
            st.session_state["_demo_data_ready"] = False
            _clear_cache()
            st.rerun()
        protein_choice = st.selectbox(
            "Demo Protein",
            list(PROTEINS.keys()),
            index=list(PROTEINS.keys()).index(
                st.session_state.get("_selected_protein", "Cytochrome C")
            ),
        )
        if protein_choice != st.session_state.get("_selected_protein"):
            st.session_state["_selected_protein"] = protein_choice
            st.session_state["_demo_data_ready"] = False
            _clear_cache()
            st.rerun()

    if data_source == "Sample mzML (Cytochrome C)":
        if st.session_state.get("_data_source") != "sample_mzml":
            _clear_cache()
            _load_sample_data()
        else:
            st.caption("Loaded: Cytochrome C (4 MS1 + 6 MS2, Orbitrap ETD)")

    if data_source == "Upload mzML":
        uploaded = st.file_uploader("Upload mzML", type=["mzML"])
        use_flashdeconv = st.checkbox("Run FLASHDeconv", value=False,
                                      help="Deconvolve MS1 scans to neutral masses")
        if use_flashdeconv:
            min_charge = st.slider("Min charge", 1, 10, 1)
            max_charge = st.slider("Max charge", 10, 100, 50)
        else:
            min_charge, max_charge = 1, 50

        if uploaded and st.button("Process", type="primary"):
            from pyopenms_integration import load_mzml, run_flashdeconv
            file_bytes = uploaded.read()
            with st.spinner("Processing mzML..."):
                if use_flashdeconv:
                    result = run_flashdeconv(file_bytes, min_charge, max_charge)
                else:
                    result = load_mzml(file_bytes)
            for key, value in result.items():
                st.session_state[key] = value
            st.session_state["_demo_data_ready"] = True
            st.rerun()

pages.run()
