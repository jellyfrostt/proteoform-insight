"""FLASHTnT output parser and coverage computation.

Adapted from FLASHApp's src/parse/tnt.py for parsing FLASHTnT TSV output files.
"""

import shutil
import subprocess
from pathlib import Path

import polars as pl


def parse_flashtnt_output(
    protein_tsv: str | Path,
    tag_tsv: str | Path,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Parse FLASHTnT protein.tsv and tags.tsv output files.

    Args:
        protein_tsv: Path to *_protein.tsv from FLASHTnT.
        tag_tsv: Path to *_tags.tsv from FLASHTnT.

    Returns:
        (protein_df, tag_df) as polars LazyFrames.
    """
    protein_df = pl.scan_csv(str(protein_tsv), separator="\t")
    tag_df = pl.scan_csv(str(tag_tsv), separator="\t")

    # Standardize column names (FLASHTnT uses spaces/caps inconsistently)
    rename_protein = {
        "Proteoform Index": "ProteoformIndex",
        "Protein Accession": "ProteinAccession",
        "Database Sequence": "DatabaseSequence",
        "Protein Sequence": "ProteinSequence",
        "Start Position": "StartPosition",
        "End Position": "EndPosition",
        "Mod Count": "ModCount",
        "Mod Mass": "ModMass",
        "Mod Start": "ModStart",
        "Mod End": "ModEnd",
        "Mod ID": "ModID",
        "Coverage (%)": "Coverage(%)",
        "Proteoform Mass": "ProteoformMass",
        "Proteoform Level Qvalue": "ProteoformLevelQvalue",
    }

    rename_tag = {
        "Tag Index": "TagIndex",
        "Proteoform Index": "ProteoformIndex",
        "Tag Sequence": "TagSequence",
        "Start Position": "StartPosition",
        "DeNovo Score": "DeNovoScore",
    }

    # Only rename columns that exist
    protein_cols = protein_df.collect_schema().names()
    tag_cols = tag_df.collect_schema().names()

    protein_renames = {k: v for k, v in rename_protein.items() if k in protein_cols}
    tag_renames = {k: v for k, v in rename_tag.items() if k in tag_cols}

    if protein_renames:
        protein_df = protein_df.rename(protein_renames)
    if tag_renames:
        tag_df = tag_df.rename(tag_renames)

    # Expand semicolon-delimited multi-proteoform tags
    if "ProteoformIndex" in tag_df.collect_schema().names():
        tag_collected = tag_df.collect()
        # Check if any ProteoformIndex values contain semicolons
        pf_col = tag_collected["ProteoformIndex"]
        if pf_col.dtype == pl.Utf8:
            # Split semicolon-delimited values and explode
            tag_collected = tag_collected.with_columns(
                pl.col("ProteoformIndex").str.split(";")
            ).explode("ProteoformIndex").with_columns(
                pl.col("ProteoformIndex").cast(pl.Int64)
            )
        tag_df = tag_collected.lazy()

    # Adjust to 0-based positions if they appear 1-based
    # (FLASHTnT uses 1-based in some versions)
    # We detect by checking if StartPosition min is 1
    protein_schema = protein_df.collect_schema()
    if "StartPosition" in protein_schema.names():
        sample = protein_df.select("StartPosition").head(1).collect()
        if len(sample) > 0 and sample["StartPosition"][0] == 1:
            protein_df = protein_df.with_columns(
                (pl.col("StartPosition") - 1).alias("StartPosition"),
                (pl.col("EndPosition") - 1).alias("EndPosition"),
            )

    tag_schema = tag_df.collect_schema()
    if "StartPosition" in tag_schema.names():
        sample = tag_df.select("StartPosition").head(1).collect()
        if len(sample) > 0 and sample["StartPosition"][0] >= 1:
            tag_df = tag_df.with_columns(
                (pl.col("StartPosition") - 1).alias("StartPosition"),
            )

    return protein_df, tag_df


def run_flashtnt_cli(
    deconv_mzml: str | Path,
    fasta: str | Path,
    output_dir: str | Path,
) -> dict[str, Path] | None:
    """Run FLASHTnT CLI tool if available.

    Args:
        deconv_mzml: Path to deconvolved mzML file.
        fasta: Path to FASTA database file.
        output_dir: Directory for output files.

    Returns:
        Dict with 'protein' and 'tag' paths, or None if FLASHTnT not found.
    """
    if not shutil.which("FLASHTnT"):
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "FLASHTnT",
        "-in", str(deconv_mzml),
        "-fasta", str(fasta),
        "-out_dir", str(output_dir),
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return None

    protein_path = output_dir / "protein.tsv"
    tag_path = output_dir / "tags.tsv"

    if protein_path.exists() and tag_path.exists():
        return {"protein": protein_path, "tag": tag_path}

    # Try alternative naming (some FLASHTnT versions use prefix)
    for f in output_dir.iterdir():
        if f.name.endswith("_protein.tsv"):
            protein_path = f
        elif f.name.endswith("_tags.tsv"):
            tag_path = f

    if protein_path.exists() and tag_path.exists():
        return {"protein": protein_path, "tag": tag_path}

    return None


def run_flashtnt_pyopenms(
    deconv_mzml: str | Path,
    fasta: str | Path,
    output_dir: str | Path,
) -> dict[str, Path] | None:
    """Run FLASHTnT via pyOpenMS Python bindings if available.

    Checks for oms.FLASHTnTAlgorithm (newer pyOpenMS builds).
    Returns same format as run_flashtnt_cli, or None if bindings unavailable.
    """
    try:
        import pyopenms as oms
    except ImportError:
        return None

    if not hasattr(oms, "FLASHTnTAlgorithm"):
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        algo = oms.FLASHTnTAlgorithm()
        params = algo.getDefaults()
        algo.setParameters(params)

        exp = oms.MSExperiment()
        oms.MzMLFile().load(str(deconv_mzml), exp)

        fasta_entries = []
        oms.FASTAFile().load(str(fasta), fasta_entries)

        algo.run(exp, fasta_entries)

        protein_path = output_dir / "protein.tsv"
        tag_path = output_dir / "tags.tsv"
        algo.writeOutput(str(protein_path), str(tag_path))

        if protein_path.exists() and tag_path.exists():
            return {"protein": protein_path, "tag": tag_path}
    except Exception:
        pass

    return None


def has_flashtnt_pyopenms() -> bool:
    """Check if FLASHTnT Python bindings are available."""
    try:
        import pyopenms as oms
        return hasattr(oms, "FLASHTnTAlgorithm")
    except ImportError:
        return False


def compute_tag_coverage(
    sequence: str,
    tag_df: pl.LazyFrame | pl.DataFrame,
    proteoform_index: int | None = None,
) -> list[float]:
    """Compute per-residue coverage from sequence tags.

    Args:
        sequence: Protein sequence string.
        tag_df: Tag DataFrame with StartPosition, Length columns.
        proteoform_index: If provided, filter tags to this proteoform.

    Returns:
        List of floats (0.0-1.0) per residue, normalized by max coverage.
    """
    if isinstance(tag_df, pl.LazyFrame):
        tag_df = tag_df.collect()

    if proteoform_index is not None and "ProteoformIndex" in tag_df.columns:
        tag_df = tag_df.filter(pl.col("ProteoformIndex") == proteoform_index)

    n = len(sequence)
    coverage = [0] * n

    for row in tag_df.iter_rows(named=True):
        start = row.get("StartPosition", 0)
        length = row.get("Length", 0)
        for i in range(start, min(start + length, n)):
            coverage[i] += 1

    max_cov = max(coverage) if max(coverage) > 0 else 1
    return [c / max_cov for c in coverage]
