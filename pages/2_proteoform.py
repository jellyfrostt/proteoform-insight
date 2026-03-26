"""Page 2: Proteoform Analysis -- SequenceView + Mirror Plot + Peak Table.

Sequence coverage, spectrum annotation, and interactive fragment analysis.
"""

import numpy as np
import polars as pl
import streamlit as st
from openms_insight import LinePlot, SequenceView, Table
from openms_insight.components.sequenceview import (
    _calculate_fragment_masses_simple,
    get_theoretical_mass,
)

from demo_data import PROTON
from init import init_app
from pyopenms_integration import has_pyopenms
init_app()

sm = st.session_state["state_manager"]

st.header("Proteoform Analysis")

# -- Sidebar controls -------------------------------------------------
with st.sidebar:
    st.subheader("Analysis Settings")

    sequence = st.text_area(
        "Protein Sequence",
        value=st.session_state["cytc_sequence"],
        height=100,
    )

    # Theoretical mass display
    theo_mass = get_theoretical_mass(sequence)
    st.metric("Theoretical Mass", f"{theo_mass:.4f} Da")

    ion_types = st.multiselect(
        "Ion Types", ["a", "b", "c", "x", "y", "z"],
        default=["c", "z"],
        help="c/z -- ETD/ECD (top-down); b/y -- CID/HCD (bottom-up)",
    )

    tolerance = st.slider("Tolerance (ppm)", 5, 50, 20)

    ms2_ids = st.session_state["ms2_scan_ids"]
    preselected = sm.get_selection("spectrum")
    default_idx = 0
    if preselected in ms2_ids:
        default_idx = ms2_ids.index(preselected)

    selected_scan = st.selectbox("MS2 Scan", ms2_ids, index=default_idx)
    sm.set_selection("spectrum", selected_scan)

    # Precursor info for selected MS2 scan
    _prec_info = (
        st.session_state["spectra_table"]
        .filter(pl.col("scan_id") == selected_scan)
        .collect()
    )
    if len(_prec_info) > 0 and _prec_info["precursor_mz"][0] > 0:
        _pc1, _pc2 = st.columns(2)
        _pc1.metric("Precursor m/z", f"{_prec_info['precursor_mz'][0]:.4f}")
        _pc2.metric("Charge", f"+{_prec_info['precursor_charge'][0]}")

    st.divider()
    st.subheader("Custom Masses")
    custom_masses_raw = st.text_area(
        "User-defined m/z values",
        placeholder="Enter one m/z per line, e.g.:\n687.345\n1031.52\n824.41",
        height=100,
        help="Paste observed m/z values to match against theoretical fragments.",
    )
    use_custom = st.checkbox("Use custom masses in mirror plot")

    st.divider()
    if has_pyopenms():
        st.caption("Fragments: pyOpenMS TheoreticalSpectrumGenerator")
    else:
        st.caption("Fragments: built-in calculator (install pyOpenMS for accuracy)")

# Ion colors
ION_COLORS = {
    "a": "#9B59B6", "b": "#E74C3C", "c": "#E67E22",
    "x": "#1ABC9C", "y": "#3498DB", "z": "#2ECC71",
}


# =====================================================================
# Helpers
# =====================================================================

def match_with_charge_states(obs_mz_list, fragment_data, ion_types_list,
                             max_z, tol_ppm):
    """Match observed m/z against theoretical fragments across charge states.

    Returns (n_term_positions, c_term_positions, match_details_list).
    """
    n_term = set()
    c_term = set()
    details = []

    for ion in ion_types_list:
        key = f"fragment_masses_{ion}"
        if key not in fragment_data:
            continue
        is_nterm = ion in ("a", "b", "c")

        for pos_idx, masses in enumerate(fragment_data[key]):
            ion_number = pos_idx + 1
            for neutral_mass in masses:
                if neutral_mass <= 0:
                    continue
                found = False
                for z in range(1, max_z + 1):
                    theo_mz = (neutral_mass + z * PROTON) / z
                    for obs_mz in obs_mz_list:
                        ppm = (obs_mz - theo_mz) / theo_mz * 1e6
                        if abs(ppm) < tol_ppm:
                            if is_nterm:
                                n_term.add(pos_idx)
                            else:
                                c_term.add(pos_idx)
                            details.append({
                                "Ion": f"{ion}{ion_number}",
                                "Type": "N-term" if is_nterm else "C-term",
                                "z": z,
                                "Theo m/z": round(theo_mz, 4),
                                "Obs m/z": round(obs_mz, 4),
                                "Delta ppm": round(ppm, 2),
                            })
                            found = True
                            break
                    if found:
                        break

    return n_term, c_term, details


def render_coverage_html(residues, n_matched, c_matched):
    """Render ProSight Lite-style sequence coverage HTML.

    Red = N-terminal (a/b/c), Blue = C-terminal (x/y/z), Green = both.
    """
    n = len(residues)
    html = '<div style="font-family:monospace;line-height:2.0;margin:10px 0;">'

    # N-terminal markers (above)
    html += '<div style="height:18px;">'
    for i in range(n):
        if i in n_matched:
            html += (
                '<span style="display:inline-block;width:18px;text-align:center;'
                'color:#E74C3C;font-size:10px;">&#9660;</span>'
            )
        else:
            html += '<span style="display:inline-block;width:18px;">&nbsp;</span>'
        if (i + 1) % 10 == 0:
            html += '<span style="display:inline-block;width:24px;">&nbsp;</span>'
    html += '</div>'

    # Residues
    html += '<div>'
    for i, res in enumerate(residues):
        both = i in n_matched and i in c_matched
        nterm = i in n_matched
        cterm = i in c_matched
        if both:
            bg, fg = "#27ae60", "white"
        elif nterm:
            bg, fg = "#E74C3C", "white"
        elif cterm:
            bg, fg = "#3498DB", "white"
        else:
            bg, fg = "#f0f0f0", "#aaa"
        html += (
            f'<span style="display:inline-block;width:18px;text-align:center;'
            f'background:{bg};color:{fg};font-size:13px;font-weight:bold;'
            f'border-radius:2px;">{res}</span>'
        )
        if (i + 1) % 10 == 0:
            html += (
                f'<span style="display:inline-block;width:24px;text-align:center;'
                f'color:#bbb;font-size:9px;">{i+1}</span>'
            )
    html += '</div>'

    # C-terminal markers (below)
    html += '<div style="height:18px;">'
    for i in range(n):
        if i in c_matched:
            html += (
                '<span style="display:inline-block;width:18px;text-align:center;'
                'color:#3498DB;font-size:10px;">&#9650;</span>'
            )
        else:
            html += '<span style="display:inline-block;width:18px;">&nbsp;</span>'
        if (i + 1) % 10 == 0:
            html += '<span style="display:inline-block;width:24px;">&nbsp;</span>'
    html += '</div></div>'

    return html


# =====================================================================
# Data
# =====================================================================

scan_meta = (
    st.session_state["spectra_table"]
    .filter(pl.col("scan_id") == selected_scan)
    .collect()
)
max_charge = scan_meta["precursor_charge"][0] if len(scan_meta) > 0 else 14

scan_peaks = (
    st.session_state["peaks_table"]
    .filter(pl.col("scan_id") == selected_scan)
    .collect()
)

# Parse custom masses
custom_mz_list = []
if use_custom and custom_masses_raw.strip():
    for line in custom_masses_raw.strip().splitlines():
        line = line.strip()
        if line:
            try:
                custom_mz_list.append(float(line))
            except ValueError:
                pass

# Choose mass list
if custom_mz_list:
    obs_masses = custom_mz_list
    obs_intensities = [1.0] * len(custom_mz_list)
else:
    obs_masses = scan_peaks["mass"].to_list() if len(scan_peaks) > 0 else []
    obs_intensities = scan_peaks["intensity"].to_list() if len(scan_peaks) > 0 else []

# Compute fragment data and matching
if has_pyopenms() and ion_types:
    from pyopenms_integration import compute_fragments_pyopenms
    frag_data = compute_fragments_pyopenms(sequence, ion_types)
    _frag_source = "pyOpenMS TheoreticalSpectrumGenerator"
else:
    frag_data = _calculate_fragment_masses_simple(sequence)
    _frag_source = "built-in calculator"
n_term_matched, c_term_matched, match_details = match_with_charge_states(
    obs_masses, frag_data, ion_types, max_charge, tolerance,
)


# =====================================================================
# 1. ProSight Lite Coverage Map + Metrics
# =====================================================================

st.subheader("Sequence Coverage Map")

residues = list(sequence)
total_sites = len(residues) - 1
covered = n_term_matched | c_term_matched
pct = len(covered) / total_sites * 100 if total_sites > 0 else 0

# Metrics row
c1, c2, c3, c4 = st.columns(4)
c1.metric("Coverage", f"{pct:.0f}%", help=f"{len(covered)}/{total_sites} cleavage sites")
c2.metric("N-terminal", f"{len(n_term_matched)} sites")
c3.metric("C-terminal", f"{len(c_term_matched)} sites")
c4.metric("Matched Ions", f"{len(match_details)}")

if custom_mz_list:
    st.caption(f"Using {len(custom_mz_list)} user-defined m/z values.")

# ProSight Lite HTML
html = render_coverage_html(residues, n_term_matched, c_term_matched)
st.markdown(html, unsafe_allow_html=True)
st.caption(
    "&#9660; N-terminal (red) | &#9650; C-terminal (blue) | "
    "Green: both | Gray: no match"
)


# =====================================================================
# 2. Interactive SequenceView (OpenMS-Insight)
# =====================================================================

st.subheader("Interactive Fragment View")

sequence_view = SequenceView(
    cache_id=f"proteoform_sv_{selected_scan}_{'-'.join(ion_types)}_{tolerance}",
    sequence_data=(sequence, max_charge),
    peaks_data=st.session_state["peaks_table"],
    filters={"spectrum": "scan_id"},
    interactivity={"peak": "peak_id"},
    deconvolved=False,
    annotation_config={
        "ion_types": ion_types,
        "neutral_losses": False,
        "tolerance": float(tolerance),
        "tolerance_ppm": True,
        "colors": {k: v for k, v in ION_COLORS.items() if k in ion_types},
    },
    title=f"{st.session_state.get('_selected_protein', 'Protein')} -- Scan {selected_scan}",
    cache_path=".cache",
)

sv_result = sequence_view(key="proteoform_sv", state_manager=sm)


# =====================================================================
# 2b. Annotated Spectrum (linked to SequenceView)
# =====================================================================

st.subheader("Annotated Spectrum")

annotated_plot = LinePlot.from_sequence_view(
    sequence_view=sequence_view,
    cache_id=f"annotated_spec_{selected_scan}_{'-'.join(ion_types)}_{tolerance}",
    cache_path=".cache",
    title=f"Annotated Spectrum -- Scan {selected_scan}",
    x_label="m/z",
    y_label="Intensity",
    styling={
        "highlightColor": "#E4572E",
        "unhighlightedColor": "#94A3B8",
        "selectedColor": "#F59E0B",
    },
)

annotated_plot(key="annotated_spec", state_manager=sm, sequence_view_key="proteoform_sv")


# =====================================================================
# 3. Mirror Plot
# =====================================================================

st.subheader("Mirror Plot (Observed vs Theoretical)")

if obs_masses:
    # Collect theoretical masses with labels
    theo_masses = []
    theo_labels = []
    for ion in ion_types:
        key = f"fragment_masses_{ion}"
        if key in frag_data:
            for i, mass_list in enumerate(frag_data[key]):
                for m in mass_list:
                    theo_masses.append(m)
                    theo_labels.append(f"{ion}{i+1}")

    max_obs = max(obs_intensities) if obs_intensities else 1.0
    obs_norm = [i / max_obs * 100 for i in obs_intensities]

    mirror_rows = []
    peak_counter = 0

    # Top: observed peaks
    for mz, intensity in zip(obs_masses, obs_norm):
        mirror_rows.append({
            "peak_id": peak_counter, "scan_id": selected_scan,
            "mass": mz, "intensity": intensity,
            "is_matched": 0, "label": "",
        })
        peak_counter += 1

    # Bottom: theoretical ions
    for theo_m, label in zip(theo_masses, theo_labels):
        matched = False
        for z in range(1, max_charge + 1):
            theo_mz = (theo_m + z * PROTON) / z
            for obs_mz in obs_masses:
                if abs(obs_mz - theo_mz) / theo_mz * 1e6 < tolerance:
                    matched = True
                    break
            if matched:
                break

        display_mz = theo_m + PROTON
        mirror_rows.append({
            "peak_id": peak_counter, "scan_id": selected_scan,
            "mass": display_mz,
            "intensity": -70.0 if matched else -30.0,
            "is_matched": 1 if matched else 0,
            "label": label if matched else "",
        })
        peak_counter += 1

    mirror_df = pl.LazyFrame(mirror_rows)

    mirror_plot = LinePlot(
        cache_id=f"mirror_{selected_scan}_{'-'.join(ion_types)}_{tolerance}_{'custom' if custom_mz_list else 'scan'}",
        data=mirror_df,
        x_column="mass",
        y_column="intensity",
        filters={"spectrum": "scan_id"},
        interactivity={"peak": "peak_id"},
        highlight_column="is_matched",
        annotation_column="label",
        title=f"Mirror Plot -- Scan {selected_scan}",
        x_label="m/z",
        y_label="Relative Intensity",
        styling={
            "highlightColor": "#E74C3C",
            "selectedColor": "#F3A712",
            "unhighlightedColor": "#94A3B8",
        },
        cache_path=".cache",
    )

    mirror_plot(key="mirror_plot", state_manager=sm)

    n_matched_theo = sum(1 for r in mirror_rows if r["is_matched"] and r["intensity"] < 0)
    n_total_theo = sum(1 for r in mirror_rows if r["intensity"] < 0)
    if n_total_theo > 0:
        st.caption(
            f"Matched {n_matched_theo}/{n_total_theo} theoretical ions "
            f"({n_matched_theo/n_total_theo*100:.1f}%) at {tolerance} ppm tolerance"
        )
else:
    st.warning("No peaks found for the selected scan.")


# =====================================================================
# 4. Fragment Match Details + Mass Accuracy
# =====================================================================

if match_details:
    st.subheader("Fragment Match Details")

    match_df = pl.DataFrame(match_details)

    # Mass accuracy scatter plot (delta ppm vs m/z)
    scatter_rows = [
        {
            "peak_id": i,
            "scan_id": 0,
            "mass": d["Obs m/z"],
            "intensity": d["Delta ppm"],
            "is_matched": 1,
            "label": d["Ion"],
        }
        for i, d in enumerate(match_details)
    ]
    scatter_df = pl.LazyFrame(scatter_rows)

    mass_accuracy_plot = LinePlot(
        cache_id=f"mass_accuracy_{selected_scan}_{'-'.join(ion_types)}_{tolerance}",
        data=scatter_df,
        x_column="mass",
        y_column="intensity",
        highlight_column="is_matched",
        annotation_column="label",
        title="Mass Accuracy (Delta ppm vs m/z)",
        x_label="Observed m/z",
        y_label="Delta (ppm)",
        styling={
            "highlightColor": "#3B82F6",
            "unhighlightedColor": "#3B82F6",
        },
        cache_path=".cache",
    )

    mass_accuracy_plot(key="mass_accuracy", state_manager=sm)

    # Summary stats
    ppm_values = [d["Delta ppm"] for d in match_details]
    mc1, mc2, mc3 = st.columns(3)
    mc1.metric("Mean error", f"{np.mean(ppm_values):.2f} ppm")
    mc2.metric("Std dev", f"{np.std(ppm_values):.2f} ppm")
    mc3.metric("Max |error|", f"{max(abs(p) for p in ppm_values):.2f} ppm")

    # Detail table
    with st.expander(f"All matched ions ({len(match_details)})", expanded=False):
        st.dataframe(match_df.to_pandas(), height=300, use_container_width=True)

    # CSV export
    csv_data = match_df.write_csv()
    st.download_button(
        "Download Fragment Table (CSV)",
        data=csv_data,
        file_name=f"fragments_scan{selected_scan}.csv",
        mime="text/csv",
    )


# =====================================================================
# 5. Peak Table
# =====================================================================

st.subheader("Peak List")

peak_table = Table(
    cache_id=f"peaks_table_{selected_scan}",
    data=st.session_state["peaks_table"],
    filters={"spectrum": "scan_id"},
    interactivity={"peak": "peak_id"},
    column_definitions=[
        {"field": "peak_id", "title": "Peak ID", "sorter": "number"},
        {"field": "mass", "title": "m/z", "sorter": "number"},
        {"field": "intensity", "title": "Intensity", "sorter": "number"},
    ],
    title=f"Peaks -- Scan {selected_scan}",
    index_field="peak_id",
    initial_sort=[{"column": "intensity", "dir": "desc"}],
    cache_path=".cache",
)

peak_table(key="peaks_tbl", state_manager=sm)
