import subprocess
from pathlib import Path
from typing import List, Optional, Sequence

from utils.config import (
    HELIX_DSSP_CODES,
    BETA_DSSP_CODES,
    KAPPA_DSSP_CODES,
    DISORDERED_DSSP_CODES,
    GROUP_HELIX,
    GROUP_BETA,
    GROUP_KAPPA,
    GROUP_DISORDERED,
    NONDISORDERED_GROUPS,
)
from utils.classes import Residue


# ---------------------------------------------------------------------------
# DSSP secondary-structure classification
# ---------------------------------------------------------------------------

def classify_dssp_code(code: str) -> str:
    code = (code or "-").strip()
    if code in HELIX_DSSP_CODES:
        return GROUP_HELIX
    if code in BETA_DSSP_CODES:
        return GROUP_BETA
    if code in KAPPA_DSSP_CODES:
        return GROUP_KAPPA
    return GROUP_DISORDERED


def classify_foldswitch(dominant_group: str, alternative_group: str) -> Optional[str]:
    if dominant_group == alternative_group:
        return None
    if (
        dominant_group in NONDISORDERED_GROUPS and alternative_group in NONDISORDERED_GROUPS
    ):  # Helix, Beta, or Kappa
        return "switch"
    return "disordered"


def normalize_dssp_code(code: str) -> str:
    return (code or "-").strip() or "-"


def get_sequence(residues: Sequence[Residue]) -> str:
    return "".join(row.amino_acid for row in residues)


def get_dssp(residues: Sequence[Residue]) -> str:
    return "".join(normalize_dssp_code(residue.dssp_code) for residue in residues)


# ---------------------------------------------------------------------------
# Running DSSP and parsing its output
# ---------------------------------------------------------------------------

def dssp_path_helper(pdb_path: Path, dssp_path: Path) -> None:
    if dssp_path.exists() and dssp_path.stat().st_size > 0:
        return
    dssp_path.parent.mkdir(parents=True, exist_ok=True)
    for binary, args in [
        ("dssp", [str(pdb_path), str(dssp_path)]),
        ("mkdssp", ["-i", str(pdb_path), "-o", str(dssp_path)]),
    ]:
        result = subprocess.run([binary, *args], text=True, capture_output=True)
        if result.returncode == 0:
            return
    raise RuntimeError(
        f"Both dssp and mkdssp failed on {pdb_path}: {result.stdout}\n{result.stderr}"
    )


def parse_dssp(dssp_path: Path) -> List[Residue]:
    # Every relevant line in a DSSP looks like this: https://pdb-redo.eu/dssp/about
    #   XX   YY A B  C...
    # XX is entry number [3-4], identical to YY [8-9] which is residue number
    # A [11] is chain ID, B [13] is amino acid, C [16] is DSSP code (xn VW VLnk)

    residues: List[Residue] = []
    in_residue_table = False
    with dssp_path.open("r", errors="replace") as handle:
        for line in handle:
            if line.lstrip().startswith("#") and "RESIDUE" in line and "AA" in line:
                in_residue_table = True
                continue
            if not in_residue_table:
                continue
            if len(line) < 17:
                continue
            amino_acid = line[13].strip()
            if not amino_acid:
                continue
            elif amino_acid in {"!", "*"}:
                print(f"Chain Break '{amino_acid}' in {dssp_path}")
                continue
            dssp_code = line[16].strip() or "-"
            residue_number = line[5:10].strip()
            chain_id = line[11].strip()
            residues.append(
                Residue(
                    amino_acid=amino_acid.upper(),
                    dssp_code=dssp_code,
                    group=classify_dssp_code(dssp_code),
                    chain_id=chain_id,
                    residue_number=residue_number,
                )
            )
    if not residues:
        raise ValueError(f"No DSSP residue rows were parsed from {dssp_path}")
    return residues


def run_dssp(model) -> List[Residue]:
    dssp_path_helper(model.pdb_path, model.dssp_path)
    return parse_dssp(model.dssp_path)


# # ---------------------------------------------------------------------------
# # Disorder metrics
# # ---------------------------------------------------------------------------

# def disordered_region_lengths(
#     residues: Sequence[Residue], disorder_codes: Sequence[str]
# ) -> List[int]:
#     res_to_analyze = {normalize_dssp_code(code) for code in disorder_codes}
#     lengths: List[int] = []
#     current_length = 0
#     for residue in residues:
#         if normalize_dssp_code(residue.dssp_code) in res_to_analyze:
#             current_length += 1
#         elif current_length:
#             lengths.append(current_length)
#             current_length = 0
#     if current_length:
#         lengths.append(current_length)
#     return lengths


# def percent_disordered(
#     residues: Sequence[Residue], disorder_codes: Sequence[str]
# ) -> float:
#     if not residues:
#         return 0.0
#     res_to_analyze = {normalize_dssp_code(code) for code in disorder_codes}
#     n_disordered = sum(
#         1
#         for residue in residues
#         if normalize_dssp_code(residue.dssp_code) in res_to_analyze
#     )
#     return 100.0 * n_disordered / len(residues)


# def number_of_disordered_regions(
#     residues: Sequence[Residue], disorder_codes: Sequence[str]
# ) -> int:
#     return len(disordered_region_lengths(residues, disorder_codes))


# def average_disordered_region_length(
#     residues: Sequence[Residue], disorder_codes: Sequence[str]
# ) -> float:
#     lengths = disordered_region_lengths(residues, disorder_codes)
#     if not lengths:
#         return 0.0
#     return sum(lengths) / len(lengths)


# def model_disorder_metric_text(residues: Sequence[Residue]) -> str:
#     return (
#         "all="
#         f"{percent_disordered(residues, DISORDERED_DSSP_CODES):.1f}%/"
#         f"{number_of_disordered_regions(residues, DISORDERED_DSSP_CODES)}/"
#         f"{average_disordered_region_length(residues, DISORDERED_DSSP_CODES):.1f}"
#     )


# def disorder_metrics_for_residues(
#     prefix: str, residues: Sequence[Residue]
# ) -> dict:
#     return {
#         f"{prefix}_all_dssp_disorder_percent": f"{percent_disordered(residues, DISORDERED_DSSP_CODES):.2f}",
#         f"{prefix}_all_dssp_disorder_n_regions": number_of_disordered_regions(
#             residues, DISORDERED_DSSP_CODES
#         ),
#         f"{prefix}_all_dssp_disorder_avg_region_length": f"{average_disordered_region_length(residues, DISORDERED_DSSP_CODES):.2f}",
#     }


# def disorder_summary_fields(
#     dominant_residues: Sequence[Residue],
#     alternative_residues: Sequence[Residue],
# ) -> dict:
#     fields: dict = {}
#     fields.update(disorder_metrics_for_residues("dominant", dominant_residues))
#     fields.update(disorder_metrics_for_residues("putative_tier1", alternative_residues))
#     return fields


# def disorder_filter_reason(disorder_fields: dict) -> str:
#     return ""
#     # dominant_percent = float(disorder_fields["dominant_all_dssp_disorder_percent"])
#     # alternative_percent = float(disorder_fields["putative_tier1_all_dssp_disorder_percent"])
#     # maximum_percent = dominant_percent * (DISORDER_PERCENT_MAX_PERCENT_OF_DOMINANT / 100.0)
#     # if alternative_percent > maximum_percent:
#     #     return f"all_dssp_disorder_percent_above_dominant_relative_threshold:{alternative_percent:.2f}>{maximum_percent:.2f}"

#     # dominant_regions = int(disorder_fields["dominant_all_dssp_disorder_n_regions"])
#     # alternative_regions = int(disorder_fields["putative_tier1_all_dssp_disorder_n_regions"])
#     # minimum_regions = dominant_regions * (DISORDER_REGION_COUNT_MIN_PERCENT_OF_DOMINANT / 100.0)
#     # if alternative_regions < minimum_regions:
#     #     return f"all_dssp_disorder_region_count_below_dominant_relative_threshold:{alternative_regions}>{minimum_regions:.2f}"

#     # dominant_average = float(disorder_fields["dominant_all_dssp_disorder_avg_region_length"])
#     # alternative_average = float(disorder_fields["putative_tier1_all_dssp_disorder_avg_region_length"])
#     # maximum_average = dominant_average * (DISORDER_AVERAGE_REGION_LENGTH_MAX_PERCENT_OF_DOMINANT / 100.0)
#     # if alternative_average > maximum_average:
#     #     return f"all_dssp_disorder_average_region_length_above_dominant_relative_threshold:{alternative_average:.2f}>{maximum_average:.2f}"
