import subprocess
from pathlib import Path
from typing import Dict, Optional, Set

# from utils.config import CFG

# References outputs/labels assigned by by Williams et al. 2025 (phenix.barbed_wire_analysis)
# Relies on their implementation of pLDDT/packing/outlier criteria.

ALL_MODES = {
    "Predictive",
    "Unpacked high pLDDT",
    "Unpacked possible",
    "Near-predictive",
    "Pseudostructure",
    "Barbed wire",
    "Unphysical",
    "Unassigned",
} # exactly matching the phenix output labels except for "unassigned" which is custom to account for blank entries


# FILTERED_MODES: Set[str] = set(CFG.get("phenix_filtered_modes", ["barbed_wire"]))
PHENIX_FILTERED_MODES: Set[str] = {"Barbed wire"}
PHENIX_THRESHOLDS: Dict[str, float] = {
    "Barbed wire": 1.0,
    "Pseudostructure": 1.03
}

def _phenix_output_path(pdb_path: Path) -> Path:
    return pdb_path.with_suffix(".phenix")


def _run_phenix(pdb_path: Path, output_path: Path) -> bool:
    """
    Runs phenix.barbed_wire_analysis and writes per-residue mode annotations
    to output_path. Returns True on success.
    """
    if output_path.exists() and output_path.stat().st_size > 0:
        print("Cached output preview:")
        for line in output_path.read_text(errors="replace").splitlines()[:40]:
            print(line)
        return True
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["phenix.barbed_wire_analysis", str(pdb_path), "output.type=text"],
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        print("STDERROR:")
        print(result.stderr)
        return False
    # debug
    output_path.write_text(result.stdout)
    return True


def _parse_phenix_output(output_path: Path) -> Dict[str, float]:
    """
    Parses phenix.barbed_wire_analysis text output.
    Lines look like:
      A,   3 to A,  29 Unpacked possible 27
      A,  30 to A,  30 Barbed wire 1
      A,   1 to A,   2  2          <- blank label, counted as unassigned
    Format: chain, start, to, chain, end, [mode words...], count
    The range header is always 5 tokens (A, X to A, Y), mode label starts at index 5.
    """
    counts: Dict[str, int] = {mode: 0 for mode in ALL_MODES}
    total = 0
    for line in output_path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("Starting job") or line.startswith("="):
            continue
        parts = line.split()
        # minimum: A, X to A, Y count = 6 tokens
        if len(parts) < 6:
            continue
        try:
            n_residues = int(parts[-1])
        except ValueError:
            continue
        # mode label is everything between the 5-token range header and the count
        label_parts = parts[5:-1]
        if not label_parts:
            counts["Unassigned"] += n_residues
            total += n_residues
            continue
        label = " ".join(label_parts)
        if label not in ALL_MODES:
            print(f"Unrecognized phenix label: {label!r} in line: {line}")
            continue
        counts[label] += n_residues
        total += n_residues

    print(f"Residues counted from phenix output: {total}")
    for mode, count in sorted(counts.items()):
        print(f"  {mode}: {count}")
    if total == 0:
        return {mode: 0.0 for mode in ALL_MODES}
    return {mode: count / total for mode, count in counts.items()}


def get_mode_fractions(pdb_path: Path) -> Optional[Dict[str, float]]:
    """
    Returns per-mode residue fractions for a PDB file, running phenix if needed
    and caching the result alongside the PDB file.
    Returns None if phenix is unavailable or fails.
    """
    output_path = _phenix_output_path(pdb_path)
    try:
        if not _run_phenix(pdb_path, output_path):
            return None
        # debug
        # fractions = _parse_phenix_output(output_path)

        # for mode, fraction in fractions.items():
        #     print(mode, fraction)
        # debug
        return _parse_phenix_output(output_path)
    except FileNotFoundError:
        return None

def phenix_fields_for_model(prefix: str, pdb_path: Path) -> Dict[str, object]:
    fractions = get_mode_fractions(pdb_path)
    fields: Dict[str, object] = {}
    for mode in ALL_MODES:
        fields[f"{prefix}_phenix_{mode}"] = (
            f"{fractions.get(mode, 0.0):.3f}" if fractions is not None else ""
        )
    return fields


# def main():
#     pdb_dir = Path("../examples_failed_cases/E_Coli/failed/WP_000752961/cfr/blind_prediction/WP_000752961/WP_000752961_predicted_models_rand_57_max_1_ext_2/")
#     pdb_files = sorted(pdb_dir.glob("*.pdb"))
#     for pdb_path in pdb_files:
#         print(f"\nProcessing: {pdb_path.name}")

#         get_mode_fractions(pdb_path)

# if __name__ == "__main__":
#     main()