# Proteoform Insight

Interactive proteoform-centric visualization and targeted spectrum analysis for top-down proteomics. Built on [OpenMS](https://openms.de/) and [OpenMS-Insight](https://github.com/t0mdavid-m/OpenMS-Insight).

## Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

Requires **pyOpenMS >= 3.0** for the full pipeline:

- **mzML loading** — `MzMLFile` + `MSExperiment` for real mass spectrometry data
- **FLASHDeconv** — `FLASHDeconvAlgorithm` to deconvolve MS1 scans to neutral masses
- **TheoreticalSpectrumGenerator** — accurate fragment mass calculation via `AASequence`

A built-in **demo mode** with synthetic Cytochrome C data is available as a fallback when pyOpenMS is not installed.

## Pages

| Page | Description |
|------|-------------|
| **LC-MS Overview** | Heatmap (m/z × RT peak map) + spectra table. Click to select scans. |
| **Proteoform Analysis** | SequenceView with fragment annotations, mirror plot (observed vs theoretical), peak table. |
| **FLASHTnT Results** | Proteoform/tag tables from FLASHTnT output, per-residue tag coverage, sequence view. Supports upload of real TSV files. |
| **Feature Traces** | Theoretical vs observed isotope envelope heatmap with XIC extraction. |

### Features

- **Proteoform-centric sequence coverage** — SequenceView with fragment ion annotation
- **Targeted spectrum analysis** — Mirror plot comparing observed and theoretical spectra
- **FLASHTnT integration** — Parser for protein.tsv/tags.tsv with tag coverage visualization
- **Interactive exploration** — Cross-component linking via StateManager (click heatmap → select scan → view fragments)
- **Feature detection overview** — Isotope trace heatmap with theoretical overlay and XIC extraction

## Architecture

Built on **OpenMS-Insight** (Streamlit component library for mass spectrometry).

```
StateManager("proteoform_insight")
    │
    ├── "spectrum" (scan_id)
    │     ← Heatmap click, Table click, selectbox
    │     → filters: LinePlot, peaks Table, SequenceView
    │
    ├── "peak" (peak_id)
    │     ← SequenceView click, LinePlot click, Table click
    │     → highlights across proteoform page components
    │
    └── "proteoform" (ProteoformIndex)
          ← FLASHTnT proteoform Table click
          → filters: tag Table, coverage LinePlot, SequenceView
```

## Demo Data

Synthetic data for **horse Cytochrome C** (P00004, 104 AA, ~12,360 Da):

- 20 MS1 scans (RT 4.0–6.0 min) with charge envelopes z=8–25
- 5 MS2 scans (ETD fragmentation) with ~75% theoretical c/z ions + sub-5 ppm error
- 3 mock proteoforms: wild-type, N-terminal acetylation, oxidized Met65
- 15–20 sequence tags per proteoform
- Feature traces at z=12, 13, 14 with theoretical isotope envelopes

## Key Dependencies

- [pyOpenMS](https://pyopenms.readthedocs.io/) — Python bindings for OpenMS (mzML I/O, FLASHDeconv, TheoreticalSpectrumGenerator, AASequence)
- [OpenMS-Insight](https://github.com/t0mdavid-m/OpenMS-Insight) — Streamlit components for mass spectrometry (Heatmap, Table, LinePlot, SequenceView)
- [Polars](https://pola.rs/) — Data processing
- [Streamlit](https://streamlit.io/) — Web framework

## References

- [FLASHDeconv](https://doi.org/10.1016/j.celrep.2020.107961) — Ultrafast mass deconvolution
- [FLASHTnT](https://doi.org/10.1021/acs.analchem.3c04283) — Top-down proteoform identification
- [FLASHApp](https://github.com/OpenMS/FLASHApp) — Streamlit app for FLASH tools
- [OpenMS](https://openms.de/) — Open-source mass spectrometry framework
