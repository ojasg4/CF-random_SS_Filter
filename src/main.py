#!/usr/bin/env python3
"""
Entry point for the fold-switch DSSP screen.

Usage:
    python src/main.py <pse_file> --job-name <name> [options]
    python src/main.py <dir_or_pse> ... [--processes N] [options]

Run with --help for the full option list.
"""

import argparse
import os
import re
import sys
import tempfile
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from multiprocessing import Pool
from functools import partial

from utils.config import MIN_CONSECUTIVE_CHANGED_RESIDUES, RMSD_THRESHOLD
from utils.pipelines import run_pipeline_default, run_pipeline_keep_going
from utils.cli_helpers import derive_job_name, collect_pse_files
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def _init_worker():
    from utils.chainsaw_filter import _load_chainsaw_model
    _load_chainsaw_model()

def run_for_one(args, job_name: str) -> None:
    """
    Runs the pipeline for a single protein with work-dir handling.

    Important:
      - If --keep-work is set, work_dir stays:
            <work_root>/<job_name>/
      - If --keep-work is not set, work_dir is temporary and deleted after the job.
        That means any DSSP/Phenix caches written into work_dir are also deleted.

    Quiet-mode stdout/stderr are written to:
        <logs_root>/<job_name>/pipeline.log
    """
    protein_log_dir = args.logs_root / job_name
    protein_log_dir.mkdir(parents=True, exist_ok=True)

    pipeline_log_path = protein_log_dir / "pipeline.log"

    pipeline = run_pipeline_keep_going if args.keep_going else run_pipeline_default

    temp_dir = None
    if args.keep_work:
        work_dir = args.work_root / job_name
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        temp_dir = tempfile.TemporaryDirectory(prefix=f"{job_name}_dssp_work_")
        work_dir = Path(temp_dir.name)

    try:
        if args.verbose:
            print(f"[{job_name}] work_dir={work_dir}")
            if not args.keep_work:
                print(
                    f"[{job_name}] WARNING: --keep-work is not set; "
                    "loaded pdbs/DSSP/Phenix caches in work_dir will be deleted after this job."
                )
            pipeline(args, job_name, protein_log_dir, work_dir)
        else:
            with pipeline_log_path.open("w") as log_handle:
                with redirect_stdout(log_handle), redirect_stderr(log_handle):
                    print(f"[{job_name}] work_dir={work_dir}")
                    if not args.keep_work:
                        print(
                            f"[{job_name}] WARNING: --keep-work is not set; "
                            "loaded pdbs/DSSP/Phenix caches in work_dir will be deleted after this job."
                        )
                    pipeline(args, job_name, protein_log_dir, work_dir)
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()


def process_one_protein(pse_file: Path, args_dict: dict) -> tuple:
    """
    extra handler/wrapper of run_for_one when more than pse is found in file directory.
    """
    args = argparse.Namespace(**args_dict)
    args.pse_file = pse_file

    job_name = derive_job_name(pse_file)
    protein_log_dir = args.logs_root / job_name
    protein_log_dir.mkdir(parents=True, exist_ok=True)
    pipeline_log_path = protein_log_dir / "pipeline.log"

    try:
        run_for_one(args, job_name)
        return (job_name, True, "")
    except Exception as exc:
        with pipeline_log_path.open("a") as log_handle:
            log_handle.write("\n\nERROR TRACEBACK\n")
            log_handle.write(traceback.format_exc())
            log_handle.write("\n")
        return (job_name, False, str(exc))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract structures from a PSE, skip alternatives below the minimum RMSD to Dominant, compare DSSP secondary structure, and write likely/unlikely summaries."
    )
    parser.add_argument(
        "inputs",
        type=Path,
        nargs="+",
        help="One or more .pse files, or directories to search recursively for PSE files (see --pse-glob).",
    )
    parser.add_argument(
        "--job-name",
        default=None,
        help="Name used for per-run output files and the per-run log subdirectory. Only valid with a single PSE file; otherwise job names are derived from each path.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of parallel worker processes when multiple PSE files are found.",
    )
    parser.add_argument(
        "--pse-glob",
        default="*structures_of_interest.pse",
        help="Glob pattern used to find PSE files inside input directories (recursive).",
    )
    parser.add_argument(
        "--logs-root",
        type=Path,
        default=Path("logs"),
        help="Root directory for per-run DSSP logs.",
    )
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=Path("fold_switch_dssp_summary"),
        help="Directory containing likely/unlikely summary CSVs and selected colored PSE outputs.",
    )
    parser.add_argument(
        "--min-run-length",
        type=int,
        default=MIN_CONSECUTIVE_CHANGED_RESIDUES,
        help="Minimum consecutive changed residues required for a qualifying region.",
    )
    parser.add_argument(
        "--dominant-label-patterns",
        nargs="+",
        default=["Dominant", "Pred_Dominant"],
        help="Labels used to identify the dominant/reference model.",
    )
    parser.add_argument(
        "--rmsd-threshold",
        type=float,
        default=RMSD_THRESHOLD,
        help="Minimum RMSD to Dominant required before running DSSP on an alternative structure. Alternatives below this value are skipped as too similar to Dominant.",
    )
    parser.add_argument(
        "--save-folds",
        action="store_true",
        help="Save conf1.fasta, conf2.fasta, conf1.pdb, conf2.pdb, and metadata.tsv for the selected pair. sequence.fasta is also written when the two sequences are identical.",
    )
    parser.add_argument(
        "--folds-root",
        type=Path,
        default=Path("saved_folds"),
        help="Root directory for --save-folds outputs.",
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep temporary extracted PDB files for debugging.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Use the Dominant model from the input PSE and compare it against sibling [protein]_predicted_models_rand_* PDBs, starting with the largest max_* directory and stopping after the first likely hit.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress messages and errors to the terminal.",
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path("fold_switch_dssp_work"),
        help="Root directory used only when --keep-work is set.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to config.json (defaults to config/config.json from the repo root).",
    )
    parser.add_argument(
        "--chainsaw",
        action="store_true",
        help="Run Chainsaw domain analysis on the work directory after selection. "
             "Default mode only; requires --keep-work.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print("DEBUG: script started", flush=True)
    print(f"DEBUG: inputs={args.inputs}", flush=True)
    print(f"DEBUG: pse_glob={args.pse_glob}", flush=True)

    try:
        if args.min_run_length < 1:
            raise ValueError("--min-run-length must be >= 1")

        if args.rmsd_threshold < 0:
            raise ValueError("--rmsd-threshold must be >= 0")

        if args.processes < 1:
            raise ValueError("--processes must be >= 1")

        pse_files = collect_pse_files(args.inputs, args.pse_glob)
        if not pse_files:
            raise ValueError(
                f"No PSE files found. Searched inputs {args.inputs} "
                f"for pattern '{args.pse_glob}'."
            )

        if len(pse_files) > 1 and args.job_name:
            raise ValueError(
                "--job-name cannot be used with multiple PSE files; "
                "job names are derived per file."
            )

        if not args.keep_work:
            print(
                "WARNING: --keep-work is not set. DSSP/Phenix caches written "
                "under work_dir will be temporary and will not be available "
                "for later visualization.",
                file=sys.stderr,
            )

        # Single file, single process: run directly.
        if len(pse_files) == 1 and args.processes == 1:
            args.pse_file = pse_files[0]
            job_name = args.job_name if args.job_name else derive_job_name(pse_files[0])
            run_for_one(args, job_name)
            return 0

        # Multiple files or parallel processing.
        args_dict = {
            key: value
            for key, value in vars(args).items()
            if key not in {"inputs", "pse_file"}
        }

        worker = partial(process_one_protein, args_dict=args_dict)

        results = []
        initializer = _init_worker if args.chainsaw else None
        
        with Pool(processes=args.processes) as pool:
            for job_name, success, error in pool.imap_unordered(worker, pse_files):
                if success:
                    print(f"[{job_name}] round 1 completed")
                else:
                    print(
                        f"[{job_name}] FAILED: {error} "
                        f"(see {args.logs_root / job_name / 'pipeline.log'})"
                    )

                results.append((job_name, success, error))

        n_ok = sum(1 for _, success, _ in results if success)
        print(f"\nCompleted: {n_ok}/{len(results)} succeeded")

        return 0 if n_ok == len(results) else 1

    except Exception as exc:
        if args.verbose:
            print(f"ERROR: {exc}", file=sys.stderr)
            raise

        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
