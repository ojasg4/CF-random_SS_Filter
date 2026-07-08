#!/usr/bin/env python3
"""
Entry point for the fold-switch DSSP screen.

Usage:
    python src/main.py <pse_file> --job-name <name> [options]

Run with --help for the full option list.
"""

import argparse
import os
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from utils.config import MIN_CONSECUTIVE_CHANGED_RESIDUES, RMSD_THRESHOLD
from utils.pipelines import run_pipeline_default, run_pipeline_keep_going


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract structures from a PSE, skip alternatives below the minimum RMSD to Dominant, compare DSSP secondary structure, and write likely/unlikely summaries."
    )
    parser.add_argument(
        "pse_file",
        type=Path,
        help="Input .pse file containing the exact case-sensitive Dominant structure.",
    )
    parser.add_argument(
        "--job-name",
        required=True,
        help="Name used for per-run output files and the per-run log subdirectory.",
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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.min_run_length < 1:
            raise ValueError("--min-run-length must be >= 1")
        if args.rmsd_threshold < 0:
            raise ValueError("--rmsd-threshold must be >= 0")

        job_name = args.job_name
        protein_log_dir = args.logs_root / job_name
        protein_log_dir.mkdir(parents=True, exist_ok=True)

        pipeline = run_pipeline_keep_going if args.keep_going else run_pipeline_default

        temp_dir = None
        if args.keep_work:
            work_dir = args.work_root / job_name
            work_dir.mkdir(parents=True, exist_ok=True)
        else:
            temp_dir = tempfile.TemporaryDirectory(prefix=f"{job_name}_dssp_work_")
            work_dir = Path(temp_dir.name)

        # Try and finally block ensures the temp_dir is cleaned even if code fails
        try:
            if args.verbose:
                pipeline(args, job_name, protein_log_dir, work_dir)
            else:
                with open(os.devnull, "w") as devnull:
                    with redirect_stdout(devnull), redirect_stderr(devnull):
                        pipeline(args, job_name, protein_log_dir, work_dir)
        finally:
            if temp_dir is not None:
                temp_dir.cleanup() # Method of TemporaryDirectory
        return 0

    except Exception as exc:
        if args.verbose:
            print(f"ERROR: {exc}", file=sys.stderr)
            raise
        return 1

if __name__ == "__main__":
    raise SystemExit(main())

