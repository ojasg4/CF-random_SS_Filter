import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

CHAINSAW_ROOT = Path(__file__).resolve().parents[2] / "chainsaw"

TSV_COLNAMES = ["chain_id", "sequence_md5", "nres", "ndom",
                "chopping", "confidence", "time_sec"]
REPORTED_COLNAMES = TSV_COLNAMES[2:]   # drop chain_id and sequence_md5

CHAINSAW_MIN_NDOM = 2
CHAINSAW_MIN_REGION_COVERAGE_PERCENT = 50.0

# One model per process. Never pass this across a multiprocessing.Pool boundary.
_MODEL = None


def _load_chainsaw_model(model_dir: Optional[Path] = None):
    global _MODEL
    if _MODEL is not None:
        return _MODEL
    # Insert at position 0 so chainsaw's own `src` package wins over ours
    if str(CHAINSAW_ROOT) not in sys.path:
        sys.path.insert(0, str(CHAINSAW_ROOT))
    from get_predictions import load_model
    if model_dir is None:
        model_dir = CHAINSAW_ROOT / "saved_models" / "model_v3"
    _MODEL = load_model(
        model_dir=str(model_dir),
        remove_disordered_domain_threshold=0.35,
        min_ss_components=2,
        min_domain_length=30,
        post_process_domains=True,
    )
    return _MODEL


def parse_chainsaw_tsv(tsv_path: Path) -> Dict[str, dict]:
    """Returns {chain_id: row_dict}. chain_id is the PDB filename stem."""
    if not tsv_path.exists():
        return {}
    with tsv_path.open(newline="") as fh:
        return {row["chain_id"]: row for row in csv.DictReader(fh, delimiter="\t")}


def score_work_directory(
    work_dir: Path,
    logs_dir: Path,
    model_dir: Optional[Path] = None,
) -> Dict[str, dict]:
    """
    Runs Chainsaw once over every PDB in work_dir (non-recursive, matching
    Chainsaw's own directory mode). Writes chainsaw_results.tsv into work_dir.
    Returns {chain_id: row}.
    """
    tsv_path = work_dir / "chainsaw_results.tsv"
    if tsv_path.exists() and tsv_path.stat().st_size > 0:
        return parse_chainsaw_tsv(tsv_path)
    

    model = _load_chainsaw_model(model_dir)
    from get_predictions import predict, get_csv_writer, write_csv_results

    logs_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for fname in sorted(os.listdir(work_dir)):
        if Path(fname).suffix not in (".pdb", ".cif"):
            continue
        pdb_path = work_dir / fname
        try:
            results.append(
                predict(model, str(pdb_path), renumber_pdbs=False, pdbchain=None)
            )
        except Exception as exc:
            print(f"Chainsaw failed on {fname}: {exc}")

    with tsv_path.open("w", newline="") as fh:
        writer = get_csv_writer(fh)
        writer.writeheader()
        write_csv_results(writer, results)

    return parse_chainsaw_tsv(tsv_path)


def chopping_coverage_percent(row: dict) -> float:
    """
    Percent of residues covered by accepted domains.
    chopping format: domains joined by ',', discontinuous segments by '_'.
    """
    chopping = (row.get("chopping") or "").strip()
    try:
        nres = int(row.get("nres") or 0)
    except ValueError:
        return 0.0
    if not chopping or chopping == "NULL" or nres <= 0:
        return 0.0
    covered = 0
    for domain in chopping.split(","):
        for seg in domain.split("_"):
            if "-" not in seg:
                continue
            start, end = seg.split("-", 1)
            try:
                covered += int(end) - int(start) + 1
            except ValueError:
                continue
    return 100.0 * covered / nres


def passes_chainsaw_filter(row: Optional[dict]) -> bool:
    """
    Domain-count and coverage filter. DISABLED for now — always passes.
    Uncomment the body to enable.
    """
    return True
    # if row is None:
    #     return False
    # try:
    #     ndom = int(row.get("ndom") or 0)
    # except ValueError:
    #     return False
    # if ndom < CHAINSAW_MIN_NDOM:
    #     return False
    # if chopping_coverage_percent(row) < CHAINSAW_MIN_REGION_COVERAGE_PERCENT:
    #     return False
    # return True


def chainsaw_summary_fields(prefix: str, row: Optional[dict]) -> Dict[str, object]:
    """One cell per model holding nres/ndom/chopping/confidence/time_sec."""
    key = f"{prefix}_chainsaw"
    if row is None:
        return {key: ""}
    return {key: ";".join(f"{c}={row.get(c, '')}" for c in REPORTED_COLNAMES)}