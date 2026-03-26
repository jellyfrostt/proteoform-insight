"""Page 4: Feature Traces -- Theoretical vs Observed isotope traces.

Feature detection and isotope trace visualization.
"""

import polars as pl
import streamlit as st
from openms_insight import Heatmap, LinePlot

from init import init_app
init_app()

sm = st.session_state["state_manager"]

_src = st.session_state.get("_data_source", "demo")
_prot = st.session_state.get("_selected_protein", "default")
_cache_suffix = f"{_src}_{_prot}"

st.header("Feature Traces")

# -- Sidebar: view mode -----------------------------------------------
with st.sidebar:
    st.subheader("View Mode")
    view_mode = st.radio(
        "Y-axis",
        ["m/z (charge-specific)", "Neutral mass (deconvolved)"],
        index=0,
        help="m/z shows charge-resolved traces; mass shows deconvolved neutral masses.",
    )

use_mass = view_mode.startswith("Neutral")
y_col = "neutral_mass" if use_mass else "mass"
y_label = "Neutral Mass (Da)" if use_mass else "m/z"

if use_mass:
    st.caption(
        "Deconvolved neutral mass x RT -- all charge states collapse to the same mass."
    )
else:
    st.caption(
        "Theoretical isotope envelopes (red) overlaid on observed peaks (blue) "
        "in m/z x RT space."
    )

# -- Feature Trace Heatmap --------------------------------------------
trace_data = st.session_state["feature_traces"]

trace_heatmap = Heatmap(
    cache_id=f"feature_trace_heatmap_{y_col}_{_cache_suffix}",
    x_column="retention_time",
    y_column=y_col,
    intensity_column="intensity",
    data=trace_data,
    interactivity={"peak": "peak_id"},
    category_column="trace_type",
    category_colors={
        "observed": "#3B82F6",
        "theoretical": "#EF4444",
    },
    title=f"Feature Traces -- {y_label} vs RT",
    x_label="Retention Time (min)",
    y_label=y_label,
    log_scale=True,
    cache_path=".cache",
)

trace_heatmap(key="trace_hm", state_manager=sm)

# -- XIC: extracted ion chromatogram ----------------------------------
st.subheader("Extracted Ion Chromatogram (XIC)")

selected_peak = sm.get_selection("peak")

if selected_peak is not None:
    peak_row = trace_data.filter(pl.col("peak_id") == selected_peak).collect()

    if len(peak_row) > 0:
        target_val = peak_row[y_col][0]
        target_charge = peak_row["charge"][0]

        if use_mass:
            tol = 0.5  # Da window for neutral mass
            tol_label = f"+/- {tol} Da"
        else:
            tol = 0.05  # Da window for m/z
            tol_label = f"+/- {tol} Da"

        st.caption(
            f"XIC at {y_label} = {target_val:.4f} (z={target_charge}, "
            f"window = {tol_label})"
        )

        xic_data = (
            trace_data
            .filter(
                (pl.col(y_col) > target_val - tol)
                & (pl.col(y_col) < target_val + tol)
            )
            .collect()
            .sort("retention_time")
        )

        if len(xic_data) > 0:
            xic_plot = LinePlot(
                cache_id=f"xic_{selected_peak}_{y_col}",
                data=xic_data.lazy(),
                x_column="retention_time",
                y_column="intensity",
                highlight_column="trace_type",
                title=f"XIC -- {y_label} {target_val:.4f}",
                x_label="Retention Time (min)",
                y_label="Intensity",
                styling={
                    "highlightColor": "#EF4444",
                    "unhighlightedColor": "#3B82F6",
                },
                cache_path=".cache",
            )

            xic_plot(key="xic_plot", state_manager=sm)
        else:
            st.info("No data points found in the XIC window.")
    else:
        st.info("Click a point on the heatmap to extract its ion chromatogram.")
else:
    st.info("Click a point on the heatmap to extract its ion chromatogram.")

# -- Info -------------------------------------------------------------
with st.expander("About Feature Traces"):
    st.markdown("""
**Feature traces** show the isotope envelopes of detected proteoforms across
retention time.

- **Blue (observed):** Peaks detected in the raw data
- **Red (theoretical):** Expected isotope positions based on the proteoform mass
  and charge state

**View modes:**
- **m/z view:** Shows charge-resolved traces. Each charge state (z=12, 13, 14)
  appears at different m/z positions. Useful for validating individual charge
  states and detecting interference.
- **Neutral mass view:** Deconvolved view where all charge states collapse to
  the same neutral mass. Useful for quantitative validation -- if deconvolution
  is correct, theoretical and observed traces should overlap precisely.

Charge states shown: z = 12, 13, 14 for Cytochrome C (~12360 Da).

Good overlap between theoretical and observed traces indicates correct
deconvolution and feature detection. Deviations may indicate:
- Mass calibration errors
- Co-eluting species
- Charge state misassignment
""")
