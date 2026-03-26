"""Page 1: LC-MS Overview -- Heatmap (m/z x RT) + Spectra Table.

Interactive exploration and feature detection overview.
"""

import streamlit as st
from openms_insight import Heatmap, Table

from init import init_app
init_app()

sm = st.session_state["state_manager"]

st.header("LC-MS Overview")
st.caption("Peak map heatmap and scan table -- click to select a scan.")

# Data-source-aware cache key to avoid stale renders
_src = st.session_state.get("_data_source", "demo")
_prot = st.session_state.get("_selected_protein", "default")
_cache_suffix = f"{_src}_{_prot}"

# -- Heatmap: m/z x RT peak map --------------------------------------
heatmap = Heatmap(
    cache_id=f"overview_heatmap_{_cache_suffix}",
    x_column="retention_time",
    y_column="mass",
    intensity_column="intensity",
    data=st.session_state["heatmap_data"],
    interactivity={"spectrum": "scan_id"},
    colorscale="Portland",
    log_scale=True,
    title="Peak Map (m/z x RT)",
    x_label="Retention Time (min)",
    y_label="m/z",
    intensity_label="Intensity",
    cache_path=".cache",
)

heatmap(key="overview_hm", state_manager=sm)

# -- Spectra Table ----------------------------------------------------
st.subheader("Spectra")

spectra_table = Table(
    cache_id=f"spectra_table_{_cache_suffix}",
    data=st.session_state["spectra_table"],
    interactivity={"spectrum": "scan_id"},
    column_definitions=[
        {"field": "scan_id", "title": "Scan", "sorter": "number"},
        {"field": "rt", "title": "RT (min)", "sorter": "number"},
        {"field": "ms_level", "title": "MS Level", "sorter": "number"},
        {"field": "precursor_mz", "title": "Precursor m/z", "sorter": "number"},
        {"field": "precursor_charge", "title": "Charge", "sorter": "number"},
        {"field": "num_peaks", "title": "Peaks", "sorter": "number"},
    ],
    title="Scan List",
    index_field="scan_id",
    default_row=0,
    initial_sort=[{"column": "scan_id", "dir": "asc"}],
    cache_path=".cache",
)

spectra_table(key="spectra_tbl", state_manager=sm)

# -- Selected scan info -----------------------------------------------
selected = sm.get_selection("spectrum")
if selected is not None:
    scan_data = (
        st.session_state["spectra_table"]
        .filter(__import__("polars").col("scan_id") == selected)
        .collect()
    )
    if len(scan_data) > 0:
        row = scan_data.row(0, named=True)
        level = row["ms_level"]
        col1, col2, col3 = st.columns(3)
        col1.metric("Selected Scan", row["scan_id"])
        col2.metric("MS Level", level)
        col3.metric("Peaks", row["num_peaks"])
        if level == 2:
            st.info(
                "This is an MS2 scan. Navigate to **Proteoform Analysis** "
                "to see fragment annotations and mirror plot."
            )
