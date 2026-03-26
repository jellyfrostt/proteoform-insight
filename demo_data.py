"""Synthetic data generation for Proteoform Insight demo mode.

Generates realistic top-down proteomics data for horse Cytochrome C
(104 AA, ~12360 Da) without requiring pyOpenMS.
"""

import numpy as np
import polars as pl
import streamlit as st

# Horse Cytochrome C sequence (104 residues)
CYTC_SEQUENCE = (
    "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKE"
    "ETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
)

# Human Ubiquitin sequence (76 residues, ~8565 Da, UniProt P62988)
UBIQ_SEQUENCE = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
    "KESTLHLVLRLRGG"
)

# Protein registry for demo mode
PROTEINS = {
    "Cytochrome C": {
        "sequence": CYTC_SEQUENCE,
        "accession": "P00004",
        "ms1_charges": list(range(8, 26)),
        "ms2_charges": [12, 13, 14, 15, 16],
        "ms2_rts": [4.5, 4.8, 5.0, 5.2, 5.5],
        "trace_charges": [12, 13, 14],
        "mod_specs": [
            {"name": "Acetyl", "mass": 42.0106, "pos": "0"},
            {"name": "Oxidation", "mass": 15.9949, "pos": "64"},
        ],
    },
    "Ubiquitin": {
        "sequence": UBIQ_SEQUENCE,
        "accession": "P62988",
        "ms1_charges": list(range(5, 16)),
        "ms2_charges": [7, 8, 9, 10, 11],
        "ms2_rts": [3.8, 4.2, 4.6, 5.0, 5.4],
        "trace_charges": [7, 8, 9],
        "mod_specs": [
            {"name": "Acetyl", "mass": 42.0106, "pos": "0"},
            {"name": "Phospho", "mass": 79.9663, "pos": "20"},
        ],
    },
}

# Amino acid monoisotopic masses
AA_MASSES = {
    "A": 71.037114, "R": 156.101111, "N": 114.042927, "D": 115.026943,
    "C": 103.009185, "E": 129.042593, "Q": 128.058578, "G": 57.021464,
    "H": 137.058912, "I": 113.084064, "L": 113.084064, "K": 128.094963,
    "M": 131.040485, "F": 147.068414, "P": 97.052764, "S": 87.032028,
    "T": 101.047679, "U": 150.953633, "W": 186.079313, "Y": 163.063329,
    "V": 99.068414,
}

H2O = 18.010565
PROTON = 1.007276

ION_OFFSETS = {
    "a": -27.994915, "b": 0.0, "c": 17.026549,
    "x": 43.989829, "y": 18.010565, "z": 1.991841,
}


def _protein_mass(sequence: str) -> float:
    """Calculate monoisotopic neutral mass of a protein."""
    return sum(AA_MASSES.get(aa, 0.0) for aa in sequence) + H2O


def _mz(mass: float, charge: int) -> float:
    """Convert neutral mass to m/z."""
    return (mass + charge * PROTON) / charge


def _fragment_masses(sequence: str, ion_types: list[str]) -> dict[str, list[float]]:
    """Calculate fragment neutral masses for given ion types."""
    n = len(sequence)
    prefix = []
    m = 0.0
    for aa in sequence:
        m += AA_MASSES.get(aa, 0.0)
        prefix.append(m)

    suffix = []
    m = 0.0
    for aa in reversed(sequence):
        m += AA_MASSES.get(aa, 0.0)
        suffix.append(m)
    suffix = list(reversed(suffix))

    result = {}
    for ion in ion_types:
        masses = []
        if ion in ("a", "b", "c"):
            for i in range(1, n):  # skip full-length
                masses.append(prefix[i - 1] + ION_OFFSETS[ion])
        else:
            for i in range(1, n):
                masses.append(suffix[n - i] + ION_OFFSETS[ion])
        result[ion] = masses
    return result


def ensure_demo_data():
    """Generate all demo data once and store in session_state. Idempotent.

    Regenerates when the selected protein changes.
    """
    selected = st.session_state.get("_selected_protein", "Cytochrome C")
    current = st.session_state.get("_current_demo_protein")

    if st.session_state.get("_demo_data_ready") and current == selected:
        return

    protein_info = PROTEINS.get(selected, PROTEINS["Cytochrome C"])
    seq = protein_info["sequence"]
    accession = protein_info["accession"]

    rng = np.random.default_rng(42)
    protein_mass = _protein_mass(seq)

    # -- MS1 scans: charge envelope across RT -------------------------
    n_ms1 = 20
    rts_ms1 = np.linspace(4.0, 6.0, n_ms1)
    charges = protein_info["ms1_charges"]

    spectra_rows = []
    peaks_rows = []
    heatmap_rows = []
    peak_id = 0

    # Gaussian elution profile centered at RT=5.0
    elution = np.exp(-0.5 * ((rts_ms1 - 5.0) / 0.4) ** 2)

    for scan_idx, rt in enumerate(rts_ms1):
        scan_id = scan_idx + 1
        scan_peaks = 0
        base_intensity = elution[scan_idx] * 1e6

        for z in charges:
            mz_val = _mz(protein_mass, z)
            # Add isotope envelope (5 peaks per charge state)
            for iso in range(5):
                iso_mz = mz_val + iso * (1.003355 / z)
                iso_int = base_intensity * (0.8 ** iso) * rng.uniform(0.8, 1.2)
                # Add sub-5 ppm error
                ppm_err = rng.normal(0, 1.5)
                iso_mz *= (1 + ppm_err * 1e-6)

                peaks_rows.append({
                    "peak_id": peak_id, "scan_id": scan_id,
                    "mass": round(iso_mz, 6), "intensity": round(iso_int, 1),
                })
                heatmap_rows.append({
                    "peak_id": peak_id, "scan_id": scan_id,
                    "retention_time": round(rt, 4), "mass": round(iso_mz, 6),
                    "intensity": round(iso_int, 1),
                })
                peak_id += 1
                scan_peaks += 1

        # Noise peaks
        n_noise = rng.integers(20, 50)
        for _ in range(n_noise):
            noise_mz = rng.uniform(400, 2000)
            noise_int = rng.exponential(5000)
            peaks_rows.append({
                "peak_id": peak_id, "scan_id": scan_id,
                "mass": round(noise_mz, 6), "intensity": round(noise_int, 1),
            })
            heatmap_rows.append({
                "peak_id": peak_id, "scan_id": scan_id,
                "retention_time": round(rt, 4), "mass": round(noise_mz, 6),
                "intensity": round(noise_int, 1),
            })
            peak_id += 1
            scan_peaks += 1

        spectra_rows.append({
            "scan_id": scan_id, "rt": round(rt, 4), "ms_level": 1,
            "precursor_mz": 0.0, "precursor_charge": 0, "num_peaks": scan_peaks,
        })

    # -- MS2 scans: ETD fragmentation ---------------------------------
    ms2_rts = protein_info["ms2_rts"]
    ms2_charges = protein_info["ms2_charges"]
    frags = _fragment_masses(seq, ["c", "z"])

    for i, (rt, prec_z) in enumerate(zip(ms2_rts, ms2_charges)):
        scan_id = n_ms1 + i + 1
        scan_peaks = 0
        prec_mz = _mz(protein_mass, prec_z)

        # ~75% of theoretical c/z ions
        for ion_type, ion_masses in frags.items():
            for frag_mass in ion_masses:
                if rng.random() > 0.75:
                    continue
                # Random charge states for each fragment
                max_frag_z = min(prec_z, max(1, int(frag_mass / 500)))
                for z in range(1, max_frag_z + 1):
                    if rng.random() > 0.5 and z > 1:
                        continue
                    frag_mz = _mz(frag_mass, z)
                    frag_int = rng.exponential(50000) * (1.0 / z)
                    ppm_err = rng.normal(0, 1.5)
                    frag_mz *= (1 + ppm_err * 1e-6)

                    peaks_rows.append({
                        "peak_id": peak_id, "scan_id": scan_id,
                        "mass": round(frag_mz, 6), "intensity": round(frag_int, 1),
                    })
                    heatmap_rows.append({
                        "peak_id": peak_id, "scan_id": scan_id,
                        "retention_time": round(rt, 4), "mass": round(frag_mz, 6),
                        "intensity": round(frag_int, 1),
                    })
                    peak_id += 1
                    scan_peaks += 1

        # MS2 noise
        n_noise = rng.integers(30, 80)
        for _ in range(n_noise):
            noise_mz = rng.uniform(100, 2000)
            noise_int = rng.exponential(2000)
            peaks_rows.append({
                "peak_id": peak_id, "scan_id": scan_id,
                "mass": round(noise_mz, 6), "intensity": round(noise_int, 1),
            })
            heatmap_rows.append({
                "peak_id": peak_id, "scan_id": scan_id,
                "retention_time": round(rt, 4), "mass": round(noise_mz, 6),
                "intensity": round(noise_int, 1),
            })
            peak_id += 1
            scan_peaks += 1

        spectra_rows.append({
            "scan_id": scan_id, "rt": round(rt, 4), "ms_level": 2,
            "precursor_mz": round(prec_mz, 6), "precursor_charge": prec_z,
            "num_peaks": scan_peaks,
        })

    # -- FLASHTnT mock data: 3 proteoforms ----------------------------
    end_pos = len(seq) - 1
    mod_specs = protein_info["mod_specs"]

    proteoforms = [
        {
            "ProteoformIndex": 0,
            "ProteinAccession": accession,
            "DatabaseSequence": seq,
            "ProteinSequence": seq,
            "StartPosition": 0, "EndPosition": end_pos,
            "ModCount": 0, "ModMass": "", "ModStart": "", "ModEnd": "", "ModID": "",
            "Coverage(%)": 78.5, "Score": 245.3,
            "ProteoformMass": round(protein_mass, 4),
            "ProteoformLevelQvalue": 0.001,
        },
    ]
    # Add modified proteoforms from mod_specs
    for idx, mod in enumerate(mod_specs, start=1):
        proteoforms.append({
            "ProteoformIndex": idx,
            "ProteinAccession": accession,
            "DatabaseSequence": seq,
            "ProteinSequence": seq,
            "StartPosition": 0, "EndPosition": end_pos,
            "ModCount": 1, "ModMass": str(mod["mass"]),
            "ModStart": mod["pos"], "ModEnd": mod["pos"],
            "ModID": mod["name"],
            "Coverage(%)": round(78.5 - idx * 13.3, 1),
            "Score": round(245.3 - idx * 55.6, 1),
            "ProteoformMass": round(protein_mass + mod["mass"], 4),
            "ProteoformLevelQvalue": round(0.001 * (idx * 4 + 1), 3),
        })

    # Sequence tags: 15-20 per proteoform
    tag_rows = []
    tag_idx = 0
    for pf in proteoforms:
        n_tags = rng.integers(15, 21)
        pf_seq = pf["ProteinSequence"]
        for _ in range(n_tags):
            length = rng.integers(4, 9)
            start = rng.integers(0, len(pf_seq) - length)
            tag_rows.append({
                "TagIndex": tag_idx,
                "ProteoformIndex": pf["ProteoformIndex"],
                "TagSequence": pf_seq[start:start + length],
                "StartPosition": int(start),
                "Length": int(length),
                "DeNovoScore": round(float(rng.uniform(50, 200)), 2),
            })
            tag_idx += 1

    # -- Feature traces: theoretical isotope envelopes ----------------
    trace_rows = []
    trace_charges = protein_info["trace_charges"]
    trace_peak_id = 0

    for z in trace_charges:
        base_mz = _mz(protein_mass, z)
        for iso in range(6):
            theo_mz = base_mz + iso * (1.003355 / z)
            for rt_idx, rt in enumerate(rts_ms1):
                intensity = elution[rt_idx] * 1e6 * (0.75 ** iso)
                # Neutral mass for this isotope peak
                theo_neutral = protein_mass + iso * 1.003355
                # Theoretical trace
                trace_rows.append({
                    "peak_id": trace_peak_id, "scan_id": rt_idx + 1,
                    "retention_time": round(rt, 4),
                    "mass": round(theo_mz, 6),
                    "neutral_mass": round(theo_neutral, 6),
                    "intensity": round(intensity, 1),
                    "trace_type": "theoretical", "charge": z,
                })
                trace_peak_id += 1
                # Observed trace (with noise and ppm error)
                obs_mz = theo_mz * (1 + rng.normal(0, 2.0) * 1e-6)
                obs_neutral = obs_mz * z - z * PROTON
                obs_int = intensity * rng.uniform(0.6, 1.1)
                trace_rows.append({
                    "peak_id": trace_peak_id, "scan_id": rt_idx + 1,
                    "retention_time": round(rt, 4),
                    "mass": round(obs_mz, 6),
                    "neutral_mass": round(obs_neutral, 6),
                    "intensity": round(obs_int, 1),
                    "trace_type": "observed", "charge": z,
                })
                trace_peak_id += 1

    # -- Store in session_state ----------------------------------------
    st.session_state["spectra_table"] = pl.LazyFrame(spectra_rows)
    st.session_state["peaks_table"] = pl.LazyFrame(peaks_rows)
    st.session_state["heatmap_data"] = pl.LazyFrame(heatmap_rows)
    st.session_state["protein_df"] = pl.LazyFrame(proteoforms)
    st.session_state["tag_df"] = pl.LazyFrame(tag_rows)
    st.session_state["feature_traces"] = pl.LazyFrame(trace_rows)
    st.session_state["cytc_sequence"] = seq
    st.session_state["ms2_scan_ids"] = [n_ms1 + i + 1 for i in range(len(ms2_rts))]
    st.session_state["_current_demo_protein"] = selected
    st.session_state["_demo_data_ready"] = True
