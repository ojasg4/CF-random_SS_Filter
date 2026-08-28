import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence
import csv
import fcntl

from utils.classes import ExtractedModel
from utils.dssp import normalize_dssp_code
from utils.utils import relative_path_for_summary
from .recovery_classes import RecoveryResult, ResidueSS

# DSSP Helpers for input parsing

def parse_dssp_with_bridges(dssp_path: Path) -> List[ResidueSS]:
    rows: List[ResidueSS] = []
    in_table = False
    with dssp_path.open("r", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#") and "RESIDUE" in line and "AA" in line:
                in_table = True
                continue
            if not in_table or len(line) < 17:
                continue
            aa = line[13].strip()
            if not aa or aa in {"!", "*"}:
                continue
            bp1 = int(line[25:29].strip()) if len(line) >= 29 and line[25:29].strip().lstrip("-").isdigit() else 0
            bp2 = int(line[29:33].strip()) if len(line) >= 33 and line[29:33].strip().lstrip("-").isdigit() else 0
            rows.append(ResidueSS(
                position_index=len(rows) + 1,
                amino_acid=aa.upper(),
                dssp_code=line[16].strip() or "-",
                chain_id=line[11].strip(),
                residue_number=line[5:10].strip(),
                bp1=bp1, bp2=bp2,
            ))
    if not rows:
        raise ValueError(f"No DSSP rows parsed from {dssp_path}")
    return rows


def load_residues(model: ExtractedModel) -> List[ResidueSS]:
    if not model.dssp_path.exists() or model.dssp_path.stat().st_size == 0:
        raise ValueError(
            f"Missing cached DSSP for {model.model_name} at {model.dssp_path}. "
            f"Run the main pipeline with --keep-work first."
        )
    return parse_dssp_with_bridges(model.dssp_path)


def dssp_string(residues: Sequence[ResidueSS]) -> str:
    return "".join(normalize_dssp_code(r.dssp_code) for r in residues)


def sequence_from_residues(residues: Sequence[ResidueSS]) -> str:
    return "".join(r.amino_acid for r in residues)


def write_recovery_log(
    log_path: Path,
    job_name: str,
    input_pse: Path,
    dominant_name: str,
    assignments: Dict[str, List[ResidueSS]],
    results_by_model: Dict[str, RecoveryResult],
    selected: Optional[RecoveryResult],
    args,
    base_dir: Path,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as fh:
        fh.write(f"# job={job_name}\n")
        fh.write(f"# dominant={dominant_name}\n")
        fh.write(f"# selected={selected.alternative_model if selected else 'none'}\n")
        fh.write(f"# selected_detector={selected.detector_hit if selected else 'none'}\n")
        fh.write(f"# bridge_reshuffle_enabled={args.beta_bridges}\n")
        fh.write(f"# min_run_length={args.min_run_length}\n")
        ordered = list(assignments.keys())
        if not ordered:
            fh.write("# no_successful_dssp_assignments\n")
            return
        fh.write(sequence_from_residues(assignments[ordered[0]]) + "\n")
        for name in ordered:
            structure = dssp_string(assignments[name])
            if name == dominant_name:
                fh.write(f"{structure}\tmodel={name} role=dominant\n")
                continue
            result = results_by_model.get(name)
            if result is None:
                fh.write(f"{structure}\tmodel={name} not_evaluated\n")
                continue
            fh.write(
                f"{structure}\tmodel={name} "
                f"detector={result.detector_hit or 'none'} "
                f"bridge_reshuffled={result.bridge_reshuffled_residue_count} "
                f"disorder_changed={result.disorder_change_residue_count} "
                f"rigid_body_regions={len(result.rigid_body_regions)} "
                f"reason={result.invalid_reason or '-'}\n"
            )

# CSV

OUTPUT_FILES = {
    "recovered": "recovered.csv",
    "failed":    "recovery_failed.csv",
}

OUTPUT_COLUMNS = [
    "input_pse", "alternative_model", "source", "classification",
    "recovery_detector",
    "bridge_reshuffled_residue_count", "bridge_reshuffled_regions",
    "disorder_change_regions", "disorder_change_residue_count",
    "rigid_body_regions", "rigid_body_max_window_distance",
    "sequences_identical", "n_candidates_analyzed",
    "log_file", "raw_dssp_dir", "invalid_reason", "dom_chainsaw", "alt_chainsaw"
]

def result_to_row(
    result: Optional[RecoveryResult],
    input_pse: Path,
    classification: str,
    n_analyzed: int,
    log_file: Path,
    raw_dssp_dir: Path,
    base_dir: Path,
    invalid_reason: str = "",
    chainsaw_fields: Optional[Dict] = None
) -> Dict[str, object]:
    base = {
        "input_pse":             relative_path_for_summary(input_pse, base_dir),
        "classification":        classification,
        "n_candidates_analyzed": n_analyzed,
        "log_file":              relative_path_for_summary(log_file, base_dir),
        "raw_dssp_dir":          relative_path_for_summary(raw_dssp_dir, base_dir),
        "invalid_reason":        invalid_reason,
    }
    if result is None:
        return {
            **base,
            "alternative_model": "none",
            "source": "",
            "recovery_detector": "no_recovery",
        }
    base.update({
        "alternative_model":   result.alternative_model,
        "source":              result.source,
        "recovery_detector":   result.detector_hit or "no_recovery",
        "sequences_identical": (
            "" if result.sequences_identical is None
            else str(result.sequences_identical).lower()
        ),
        "invalid_reason":      result.invalid_reason or invalid_reason,
    })
    if result.bridge_reshuffled_regions:
        base["bridge_reshuffled_residue_count"] = result.bridge_reshuffled_residue_count
        base["bridge_reshuffled_regions"] = ",".join(
            r.annotation for r in result.bridge_reshuffled_regions
        )
    if result.disorder_change_regions:
        base["disorder_change_regions"] = ",".join(
            r.annotation for r in result.disorder_change_regions
        )
        base["disorder_change_residue_count"] = result.disorder_change_residue_count
    if result.rigid_body_regions:
        base["rigid_body_regions"] = ",".join(
            r.annotation for r in result.rigid_body_regions
        )
        base["rigid_body_max_window_distance"] = (
            f"{result.rigid_body_max_window_distance:.3f}"
        )
    if chainsaw_fields:
        base.update(chainsaw_fields)
    return base


def write_output_csvs(
    summary_dir: Path,
    recovered_rows: List[Dict],
    failed_rows: List[Dict],
) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    lock_path = summary_dir / ".recovery.lock"
    with lock_path.open("w") as lh:
        fcntl.flock(lh, fcntl.LOCK_EX)
        try:
            for filename, rows in [
                (OUTPUT_FILES["recovered"], recovered_rows),
                (OUTPUT_FILES["failed"],    failed_rows),
            ]:
                tmp = (summary_dir / filename).with_suffix(".csv.tmp")
                with tmp.open("w", newline="") as fh:
                    writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS)
                    writer.writeheader()
                    for row in rows:
                        writer.writerow({
                            col: row.get(col, "") for col in OUTPUT_COLUMNS
                        })
                tmp.replace(summary_dir / filename)
        finally:
            fcntl.flock(lh, fcntl.LOCK_UN)
