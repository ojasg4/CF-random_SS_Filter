#!/usr/bin/env python3
"""
input_only.py

Standalone pre-screen that flags PSE files whose Dominant structure exceeds
the unassigned-DSSP threshold. Writes each flagged job to too_disordered.txt
in the summary directory (job_name<TAB>percent), matching the main pipeline.

Lightweight by design: extracts only the Dominant model from each PSE, runs
DSSP on it, checks the unassigned percent. No comparison, no phenix, no
alternatives — and it does NOT import the pipeline, so startup stays fast.

If --keep-work is set, extracted Dominant PDBs and their .dssp files are
written into <work-root>/<job_name>/, so a later main.py run with the same
--work-root finds the cached DSSP and skips recomputing it.

Usage:
    python src/input_only.py <dir_or_pse> ... [options]

    python src/input_only.py \
        ../examples_failed_cases/E_Coli/failed \
        --pse-glob "structures_of_interest.pse" \
        --summary-dir fold_switch_dssp_summary \
        --work-root fold_switch_dssp_work \
        --keep-work \
        --verbose
"""

import argparse
import sys
import tempfile
from pathlib import Path
 
sys.path.insert(0, str(Path(__file__).resolve().parent))
 
from utils.cli_helpers import derive_job_name, collect_pse_files
from utils.pdb_helpers import (
    pymol_session,
    extract_models_from_pse,
    find_dominant_model,
)
from utils.dssp import run_dssp
from utils.comparison import input_is_too_disordered
 
 
def check_one(pse_file: Path, work_dir: Path, job_name: str, args, cmd) -> bool:
    """
    Extract the Dominant model into work_dir, run DSSP, and delegate the
    threshold check + logging to input_is_too_disordered.
    Returns True if the dominant was flagged as too disordered.
    """
    models = extract_models_from_pse(pse_file, work_dir, work_dir, cmd)
    dominant_model = find_dominant_model(models, args.dominant_label_patterns)
    dominant_residues = run_dssp(dominant_model)
    return input_is_too_disordered(job_name, args.summary_dir, dominant_residues)
 
 
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Flag PSE files whose Dominant structure is too disordered."
    )
    parser.add_argument("inputs", type=Path, nargs="+",
                        help="PSE files or directories to search recursively.")
    parser.add_argument("--pse-glob", default="structures_of_interest.pse",
                        help="Glob for PSE files inside input directories.")
    parser.add_argument("--summary-dir", type=Path,
                        default=Path("fold_switch_dssp_summary"),
                        help="Directory where too_disordered.txt is written.")
    parser.add_argument("--work-root", type=Path,
                        default=Path("fold_switch_dssp_work"),
                        help="Root for per-job work directories (used only with --keep-work).")
    parser.add_argument("--keep-work", action="store_true",
                        help="Persist extracted Dominant PDBs and .dssp files under "
                             "<work-root>/<job_name>/ so a later main.py run reuses the DSSP cache.")
    parser.add_argument("--dominant-label-patterns", nargs="+",
                        default=["Dominant", "Pred_Dominant"])
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()
 
 
def main() -> int:
    args = parse_args()
    try:
        pse_files = collect_pse_files(args.inputs, args.pse_glob)
        if not pse_files:
            raise ValueError(
                f"No PSE files found in {args.inputs} matching '{args.pse_glob}'."
            )
 
        n_flagged = 0
        with pymol_session() as cmd:
            for idx, pse_file in enumerate(pse_files, start=1):
                job_name = derive_job_name(pse_file)
 
                # Resolve work_dir: persistent if --keep-work, temporary otherwise
                if args.keep_work:
                    work_dir = args.work_root / job_name
                    work_dir.mkdir(parents=True, exist_ok=True)
                    tmp_ctx = None
                else:
                    tmp_ctx = tempfile.TemporaryDirectory(prefix=f"{job_name}_too_disordered_")
                    work_dir = Path(tmp_ctx.name)
 
                try:
                    flagged = check_one(pse_file, work_dir, job_name, args, cmd)
                except Exception as exc:
                    if args.verbose:
                        print(f"[{idx}/{len(pse_files)}] {job_name}: ERROR {exc}")
                    continue
                finally:
                    if tmp_ctx is not None:
                        tmp_ctx.cleanup()
 
                if flagged:
                    n_flagged += 1
                elif args.verbose:
                    print(f"[{idx}/{len(pse_files)}] {job_name}: ok")
 
        print(f"\nFlagged {n_flagged}/{len(pse_files)} as too disordered.")
        print(f"Written to {args.summary_dir / 'too_disordered.txt'}")
        return 0
 
    except Exception as exc:
        if args.verbose:
            print(f"ERROR: {exc}", file=sys.stderr)
            raise
        return 1
 
 
if __name__ == "__main__":
    raise SystemExit(main())