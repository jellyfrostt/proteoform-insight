"""pyOpenMS integration for Proteoform Insight.

Core data processing layer using pyOpenMS APIs:
- mzML file loading via MzMLFile + MSExperiment
- FLASHDeconv deconvolution on MS1 scans (FLASHDeconvAlgorithm)
- Fragment mass calculation via TheoreticalSpectrumGenerator + AASequence

All functions return data in the same schema as demo_data.py so the
visualization layer works identically with real or synthetic data.
Falls back gracefully if pyOpenMS is not installed (demo mode only).
"""

import os
import tempfile

import polars as pl


def has_pyopenms() -> bool:
    """Check if pyOpenMS is available."""
    try:
        import pyopenms  # noqa: F401
        return True
    except ImportError:
        return False


def get_pyopenms_version() -> str | None:
    """Return pyOpenMS version string, or None if not installed."""
    try:
        import pyopenms
        return pyopenms.__version__
    except (ImportError, AttributeError):
        return None


# =====================================================================
# mzML Loading
# =====================================================================

def load_mzml(file_bytes: bytes) -> dict:
    """Load an mzML file and extract spectra/peaks in session_state format.

    Parameters
    ----------
    file_bytes : bytes
        Raw bytes of the uploaded mzML file.

    Returns
    -------
    dict with keys: spectra_table, peaks_table, heatmap_data, ms2_scan_ids,
                    cytc_sequence (empty string -- user provides in page 2).
    """
    import pyopenms as oms

    exp = _load_experiment(file_bytes)

    spectra_rows = []
    peaks_rows = []
    heatmap_rows = []
    ms2_scan_ids = []
    peak_id = 0

    for spec_idx in range(exp.size()):
        spectrum = exp[spec_idx]
        rt = spectrum.getRT() / 60.0
        ms_level = spectrum.getMSLevel()

        precursor_mz = 0.0
        precursor_charge = 0
        if ms_level > 1 and spectrum.getPrecursors():
            prec = spectrum.getPrecursors()[0]
            precursor_mz = prec.getMZ()
            precursor_charge = prec.getCharge()
            ms2_scan_ids.append(spec_idx)

        mz_array, int_array = spectrum.get_peaks()

        spectra_rows.append({
            "scan_id": spec_idx, "rt": round(rt, 3), "ms_level": ms_level,
            "precursor_mz": round(precursor_mz, 4),
            "precursor_charge": precursor_charge,
            "num_peaks": len(mz_array),
        })

        for i in range(len(mz_array)):
            peaks_rows.append({
                "peak_id": peak_id, "scan_id": spec_idx,
                "mass": round(float(mz_array[i]), 4),
                "intensity": float(int_array[i]),
            })
            heatmap_rows.append({
                "peak_id": peak_id, "scan_id": spec_idx,
                "retention_time": round(rt, 4),
                "mass": round(float(mz_array[i]), 4),
                "intensity": float(int_array[i]),
            })
            peak_id += 1

    return {
        "spectra_table": pl.LazyFrame(spectra_rows),
        "peaks_table": pl.LazyFrame(peaks_rows),
        "heatmap_data": pl.LazyFrame(heatmap_rows),
        "ms2_scan_ids": ms2_scan_ids,
        "cytc_sequence": "",
    }


# =====================================================================
# FLASHDeconv
# =====================================================================

PROTON_MASS = 1.00727646677


def run_flashdeconv(file_bytes: bytes, min_charge: int = 1,
                    max_charge: int = 50) -> dict:
    """Load mzML and run FLASHDeconv on MS1 scans.

    MS1 scans are deconvolved to neutral masses. MS2 scans keep raw m/z.
    Also generates feature_traces from deconvolved data.

    Returns same schema as load_mzml, plus feature_traces LazyFrame.
    """
    import pyopenms as oms

    exp = _load_experiment(file_bytes)

    # Run FLASHDeconv
    fd = oms.FLASHDeconvAlgorithm()
    params = fd.getDefaults()
    params.setValue("SD:min_charge", min_charge)
    params.setValue("SD:max_charge", max_charge)
    fd.setParameters(params)

    deconvolved_spectra = []
    deconvolved_features = []
    fd.run(exp, deconvolved_spectra, deconvolved_features)

    deconv_by_scan = {}
    for ds in deconvolved_spectra:
        if ds.size() > 0:
            deconv_by_scan[ds.getScanNumber()] = ds

    spectra_rows = []
    peaks_rows = []
    heatmap_rows = []
    ms2_scan_ids = []
    trace_rows = []
    peak_id = 0
    trace_peak_id = 0

    for spec_idx in range(exp.size()):
        spectrum = exp[spec_idx]
        rt = spectrum.getRT() / 60.0
        ms_level = spectrum.getMSLevel()

        precursor_mz = 0.0
        precursor_charge = 0
        if ms_level > 1 and spectrum.getPrecursors():
            prec = spectrum.getPrecursors()[0]
            precursor_mz = prec.getMZ()
            precursor_charge = prec.getCharge()
            ms2_scan_ids.append(spec_idx)

        mz_array, int_array = spectrum.get_peaks()

        spectra_rows.append({
            "scan_id": spec_idx, "rt": round(rt, 3), "ms_level": ms_level,
            "precursor_mz": round(precursor_mz, 4),
            "precursor_charge": precursor_charge,
            "num_peaks": len(mz_array),
        })

        if ms_level == 1 and spec_idx in deconv_by_scan:
            # MS1: use FLASHDeconv deconvolved neutral masses
            ds = deconv_by_scan[spec_idx]
            ms_spec = ds.toSpectrum(1, 10.0, False)
            d_mz, d_int = ms_spec.get_peaks()
            for i in range(len(d_mz)):
                neutral_mass = float(d_mz[i]) - PROTON_MASS
                peaks_rows.append({
                    "peak_id": peak_id, "scan_id": spec_idx,
                    "mass": round(neutral_mass, 4),
                    "intensity": float(d_int[i]),
                })
                heatmap_rows.append({
                    "peak_id": peak_id, "scan_id": spec_idx,
                    "retention_time": round(rt, 4),
                    "mass": round(neutral_mass, 4),
                    "intensity": float(d_int[i]),
                })
                # Feature trace entry for deconvolved peaks
                trace_rows.append({
                    "peak_id": trace_peak_id, "scan_id": spec_idx,
                    "retention_time": round(rt, 4),
                    "mass": round(float(d_mz[i]), 6),
                    "neutral_mass": round(neutral_mass, 6),
                    "intensity": float(d_int[i]),
                    "trace_type": "observed", "charge": 0,
                })
                trace_peak_id += 1
                peak_id += 1
        else:
            # MS2 + non-deconvolved MS1: raw m/z
            for i in range(len(mz_array)):
                peaks_rows.append({
                    "peak_id": peak_id, "scan_id": spec_idx,
                    "mass": round(float(mz_array[i]), 4),
                    "intensity": float(int_array[i]),
                })
                heatmap_rows.append({
                    "peak_id": peak_id, "scan_id": spec_idx,
                    "retention_time": round(rt, 4),
                    "mass": round(float(mz_array[i]), 4),
                    "intensity": float(int_array[i]),
                })
                peak_id += 1

    return {
        "spectra_table": pl.LazyFrame(spectra_rows),
        "peaks_table": pl.LazyFrame(peaks_rows),
        "heatmap_data": pl.LazyFrame(heatmap_rows),
        "ms2_scan_ids": ms2_scan_ids,
        "cytc_sequence": "",
        "feature_traces": pl.LazyFrame(trace_rows),
    }


# =====================================================================
# Fragment Mass Calculation
# =====================================================================

def compute_fragments_pyopenms(sequence: str,
                               ion_types: list[str]) -> dict:
    """Compute theoretical fragment masses using pyOpenMS.

    Uses AASequence + TheoreticalSpectrumGenerator for accurate masses.

    Returns dict matching _calculate_fragment_masses_simple() format:
    {"fragment_masses_c": [[mass], [mass], ...], ...}
    """
    import pyopenms as oms

    aa_seq = oms.AASequence.fromString(sequence)
    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getDefaults()
    # Enable only the requested ion types
    for ion in ("a", "b", "c", "x", "y", "z"):
        params.setValue(f"add_{ion}_ions", "true" if ion in ion_types else "false")
    params.setValue("add_metainfo", "true")
    params.setValue("add_precursor_peaks", "false")
    params.setValue("add_losses", "false")
    tsg.setParameters(params)

    spec = oms.MSSpectrum()
    tsg.getSpectrum(spec, aa_seq, 1, 1)  # charge 1 only -> neutral + proton

    # Parse generated peaks back into per-position lists
    n = len(sequence)
    result = {}
    for ion in ion_types:
        result[f"fragment_masses_{ion}"] = [[] for _ in range(n - 1)]

    for i in range(spec.size()):
        peak = spec[i]
        mz = peak.getMZ()
        # Recover neutral mass: mz at z=1 means mz = (M + H)/1, so M = mz - H
        neutral_mass = mz - PROTON_MASS

        # Get ion annotation from metadata
        if spec.getStringDataArrays():
            annotations = spec.getStringDataArrays()
            if len(annotations) > 0 and i < annotations[0].size():
                label = annotations[0][i].decode() if isinstance(
                    annotations[0][i], bytes) else str(annotations[0][i])
                # Label format: e.g. "c3" or "z5"
                if len(label) >= 2 and label[0] in ("a", "b", "c", "x", "y", "z"):
                    ion_type = label[0]
                    try:
                        ion_number = int(label[1:].split("+")[0].split("-")[0].split("^")[0])
                    except ValueError:
                        continue
                    key = f"fragment_masses_{ion_type}"
                    if key in result and 1 <= ion_number <= n - 1:
                        pos_idx = ion_number - 1
                        result[key][pos_idx].append(neutral_mass)

    return result


# =====================================================================
# Internal helpers
# =====================================================================

def _load_experiment(file_bytes: bytes):
    """Write bytes to temp file and load as MSExperiment.

    Uses close-before-open pattern for Windows file locking compatibility.
    """
    import pyopenms as oms

    tmp = tempfile.NamedTemporaryFile(suffix=".mzML", delete=False)
    tmp_path = tmp.name
    try:
        tmp.write(file_bytes)
        tmp.close()
        exp = oms.MSExperiment()
        oms.MzMLFile().load(tmp_path, exp)
        return exp
    finally:
        os.unlink(tmp_path)
